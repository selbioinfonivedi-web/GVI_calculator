package org.gvi.algorithms.dnds;

import org.gvi.algorithms.common.CodonUtil;
import org.gvi.core.exception.GviComputationException;
import org.gvi.core.exception.GviInputException;
import org.gvi.core.util.GeneticCode;
import org.gvi.core.util.MathUtil;

import java.util.ArrayList;
import java.util.List;

/**
 * Index 5 - dN/dS Ratio (omega), Section 5.5. Implements the Nei-Gojobori
 * (1986) counting method: per-codon synonymous/nonsynonymous site counts
 * via the fractional-site method, observed differences counted by
 * enumerating every mutational pathway between two differing codons and
 * averaging synonymous/nonsynonymous steps across pathways, then a
 * Jukes-Cantor correction applied to the pooled proportions pS/pN.
 * <p>
 * {@link #compute} pools this over an entire gene into one ratio.
 * {@link #slidingWindow} applies the same math over a moving window of
 * codons, so localized positive selection isn't diluted/hidden by
 * surrounding purifying selection elsewhere in the gene -- a tractable,
 * documented stand-in for full ML site-models (e.g. MEME/FEL), not
 * equivalent to them.
 */
public final class DnDsCalculator {

    private static final char[] BASES = {'A', 'C', 'G', 'T'};
    private final NeiGojoboriSelectionTest selectionTest = new NeiGojoboriSelectionTest();

    /**
     * Fraction of the 3 codon positions' single-nucleotide substitutions
     * that are synonymous, i.e. the classic Nei-Gojobori per-codon
     * synonymous site count (range 0..3). Nonsynonymous sites = 3 - this.
     */
    public double synonymousSites(String codon) {
        char aa = GeneticCode.translate(codon);
        double sSites = 0.0;
        char[] chars = codon.toCharArray();
        for (int pos = 0; pos < 3; pos++) {
            char original = chars[pos];
            int synCount = 0;
            for (char alt : BASES) {
                if (alt == original) continue;
                char[] mutant = chars.clone();
                mutant[pos] = alt;
                char mutantAa = GeneticCode.translate(new String(mutant));
                if (mutantAa == aa) synCount++;
            }
            sSites += synCount / 3.0;
        }
        return sSites;
    }

    private record PathwayCounts(double syn, double nonsyn, boolean passesThroughInternalStop) {
    }

    /** Averages synonymous/nonsynonymous step counts across every mutational pathway between two codons. */
    PathwayResult observedDifferences(String codonA, String codonB) {
        List<Integer> diffPositions = new ArrayList<>();
        for (int i = 0; i < 3; i++) {
            if (codonA.charAt(i) != codonB.charAt(i)) diffPositions.add(i);
        }
        if (diffPositions.isEmpty()) {
            return new PathwayResult(0.0, 0.0, 1);
        }

        List<PathwayCounts> pathways = new ArrayList<>();
        for (List<Integer> order : permutations(diffPositions)) {
            char[] current = codonA.toCharArray();
            double syn = 0, nonsyn = 0;
            boolean internalStop = false;
            for (int step = 0; step < order.size(); step++) {
                int pos = order.get(step);
                char[] next = current.clone();
                next[pos] = codonB.charAt(pos);
                char aaCurrent = GeneticCode.translate(new String(current));
                char aaNext = GeneticCode.translate(new String(next));
                if (aaNext == GeneticCode.STOP && step < order.size() - 1) {
                    internalStop = true;
                }
                if (aaCurrent == aaNext) syn++;
                else nonsyn++;
                current = next;
            }
            pathways.add(new PathwayCounts(syn, nonsyn, internalStop));
        }

        List<PathwayCounts> usable = pathways.stream().filter(p -> !p.passesThroughInternalStop()).toList();
        if (usable.isEmpty()) usable = pathways; // every pathway hits a premature stop; use them anyway, nothing better available

        double avgSyn = usable.stream().mapToDouble(PathwayCounts::syn).average().orElse(0.0);
        double avgNonsyn = usable.stream().mapToDouble(PathwayCounts::nonsyn).average().orElse(0.0);
        return new PathwayResult(avgSyn, avgNonsyn, usable.size());
    }

    record PathwayResult(double synonymousDiffs, double nonsynonymousDiffs, int pathwaysUsed) {
    }

    private List<List<Integer>> permutations(List<Integer> items) {
        if (items.size() == 1) {
            List<List<Integer>> single = new ArrayList<>();
            single.add(items);
            return single;
        }
        List<List<Integer>> result = new ArrayList<>();
        for (int i = 0; i < items.size(); i++) {
            List<Integer> rest = new ArrayList<>(items);
            Integer picked = rest.remove(i);
            for (List<Integer> perm : permutations(rest)) {
                List<Integer> withPicked = new ArrayList<>();
                withPicked.add(picked);
                withPicked.addAll(perm);
                result.add(withPicked);
            }
        }
        return result;
    }

    private record Accumulation(double totalS, double totalN, double totalSd, double totalNd, int compared, int excluded) {
    }

    /** Pools site counts and observed differences over codons [startIdx, endIdxExclusive) of the two codon lists. */
    private Accumulation accumulate(List<String> refCodons, List<String> qryCodons, int startIdx, int endIdxExclusive) {
        double totalS = 0, totalN = 0, totalSd = 0, totalNd = 0;
        int compared = 0, excluded = 0;
        for (int i = startIdx; i < endIdxExclusive; i++) {
            String rc = refCodons.get(i);
            String qc = qryCodons.get(i);
            if (!GeneticCode.isStandardCodon(rc) || !GeneticCode.isStandardCodon(qc)) {
                excluded++;
                continue;
            }
            if (GeneticCode.translate(rc) == GeneticCode.STOP || GeneticCode.translate(qc) == GeneticCode.STOP) {
                excluded++;
                continue;
            }
            double sRef = synonymousSites(rc);
            double sQry = synonymousSites(qc);
            double avgS = (sRef + sQry) / 2.0;
            double avgN = 3.0 - avgS;
            totalS += avgS;
            totalN += avgN;

            PathwayResult diffs = observedDifferences(rc, qc);
            totalSd += diffs.synonymousDiffs();
            totalNd += diffs.nonsynonymousDiffs();
            compared++;
        }
        return new Accumulation(totalS, totalN, totalSd, totalNd, compared, excluded);
    }

    private record Omega(double dN, double dS, double omega, boolean jcApplied, String note) {
    }

    /** Jukes-Cantor-corrects the pooled pS/pN from an {@link Accumulation} into dN, dS, and their ratio. */
    private Omega omegaFrom(Accumulation acc) {
        String note = null;

        // Zero *observed* synonymous differences makes pS=0, dS=0, and omega=dN/dS mathematically undefined.
        // A Haldane-Anscombe-style continuity correction (add half a pseudo-difference) is applied so dS
        // becomes a small positive number and omega stays finite and downstream-safe.
        //
        // What that corrected omega is NOT is evidence of selection. This code previously asserted that
        // "zero observed synonymous divergence is itself a strong positive-selection signal" -- that is
        // backwards. Synonymous sites are roughly a quarter of a codon's sites, so on sequences with very
        // little total divergence the *expected* number of synonymous differences is near zero by chance
        // alone. Observing zero of them says the sequences are barely diverged, not that selection stripped
        // them out. Real diversifying selection is inferred from an excess of nonsynonymous change against a
        // dS estimated with actual precision -- which is exactly what is missing here.
        //
        // The magnitude of the correction is also arbitrary rather than measured: omega scales as
        // 1/(0.5/(totalS+1)), so it is set by the site count and the choice of pseudo-count, not by data.
        // It is reported as an unbounded lower bound and must not be read as a point estimate.
        // NeiGojoboriSelectionTest returns no verdict in this case, which keeps QualityGate from letting
        // this value into the composite as though it were a confirmed finding.
        double pS;
        if (acc.totalS() > 0 && acc.totalSd() == 0) {
            pS = 0.5 / (acc.totalS() + 1);
            note = "No synonymous differences were observed in " + acc.totalS() + " synonymous sites, so dS is 0 and "
                    + "dN/dS is mathematically UNDEFINED. A continuity correction (add-0.5) was applied to keep the "
                    + "ratio finite, but the resulting omega is an arbitrary lower bound set by the site count and the "
                    + "pseudo-count -- NOT a point estimate, and NOT evidence of positive selection. Zero synonymous "
                    + "change is the expected result of very low overall divergence (synonymous sites are ~1/4 of a "
                    + "codon), not a signature of selection. Treat this gene as having insufficient divergence to "
                    + "estimate omega; add more divergent sequences or a longer coding region.";
        } else {
            pS = acc.totalS() > 0 ? acc.totalSd() / acc.totalS() : 0.0;
        }
        double pN = acc.totalN() > 0 ? acc.totalNd() / acc.totalN() : 0.0;

        double dS, dN;
        boolean jcApplied;
        double jcArgS = 1.0 - (4.0 / 3.0) * pS;
        double jcArgN = 1.0 - (4.0 / 3.0) * pN;
        if (jcArgS <= 0.0 || jcArgN <= 0.0) {
            note = "Jukes-Cantor correction saturated (pS or pN too high); reporting uncorrected pN/pS ratio instead";
            dS = pS;
            dN = pN;
            jcApplied = false;
        } else {
            dS = MathUtil.stripNegativeZero(-0.75 * Math.log(jcArgS));
            dN = MathUtil.stripNegativeZero(-0.75 * Math.log(jcArgN));
            jcApplied = true;
        }

        if (dS == 0.0) {
            // totalS itself was 0 (no synonymous sites existed at all in this window/gene) -- nothing to correct.
            return new Omega(dN, dS, dN, jcApplied,
                    "No synonymous sites available in this region (totalS=0); dN/dS is undefined, reporting dN as a lower-bound proxy");
        }
        return new Omega(dN, dS, dN / dS, jcApplied, note);
    }

    public DnDsResult compute(String sequenceId, String geneName, String refCds, String queryCds) {
        if (refCds.length() != queryCds.length()) {
            throw new GviComputationException("dN/dS for '" + sequenceId + "'/'" + geneName
                    + "': reference and query CDS are different lengths (" + refCds.length() + " vs " + queryCds.length() + ")");
        }
        CodonUtil.SplitResult refSplit = CodonUtil.splitIntoCodons(refCds);
        CodonUtil.SplitResult qrySplit = CodonUtil.splitIntoCodons(queryCds);
        List<String> refCodons = refSplit.codons();
        List<String> qryCodons = qrySplit.codons();

        Accumulation acc = accumulate(refCodons, qryCodons, 0, refCodons.size());
        if (acc.compared() == 0) {
            throw new GviComputationException("dN/dS for '" + sequenceId + "'/'" + geneName + "': no comparable codons after excluding stops/invalid");
        }

        List<String> diagnostics = new ArrayList<>();
        if (acc.excluded() > 0) {
            diagnostics.add(acc.excluded() + " codon(s) excluded (stop codon or gap/ambiguous bases)");
        }
        if (refSplit.droppedTrailingBases() > 0) {
            diagnostics.add(refSplit.droppedTrailingBases() + " trailing base(s) dropped (CDS length not a multiple of 3)");
        }

        Omega omega = omegaFrom(acc);
        if (omega.note() != null) diagnostics.add(omega.note());
        diagnostics.add(selectionTest.test(acc.totalSd(), acc.totalS(), acc.totalNd(), acc.totalN()).note());

        String category = DnDsGeneReferenceTable.match(geneName)
                .map(b -> b.geneCategory() + ": " + b.selectionType() + " (typical range " + b.dnDsRange()
                        + "); " + b.biologicalRole() + "; pandemic risk: " + b.pandemicRisk())
                .orElse(DnDsGeneReferenceTable.genericClassification(omega.omega()) + " (no gene-specific reference band matched '" + geneName + "')");

        return new DnDsResult(sequenceId, geneName, omega.dN(), omega.dS(), omega.omega(), omega.jcApplied(),
                acc.totalS(), acc.totalN(), acc.totalSd(), acc.totalNd(), acc.compared(), acc.excluded(), category, diagnostics);
    }

    /**
     * Pools several {@link DnDsResult}s (e.g. one per sequence in a
     * dataset, each already computed by {@link #compute}) into ONE
     * dataset-wide dN/dS by summing their raw Nei-Gojobori site and
     * difference counts and re-deriving a single ratio from the totals --
     * <b>not</b> by averaging the individual omega ratios.
     * <p>
     * This distinction matters: averaging ratios is statistically biased
     * (Jensen's-inequality-type distortion, and especially misleading when
     * some inputs have zero synonymous differences, an edge case whose
     * continuity-corrected omega can be a large number that would dominate
     * a naive mean). Pooling counts first and computing one ratio from the
     * totals is the standard, correct way population-genetics tools
     * combine multiple pairwise comparisons into one summary statistic.
     */
    public DnDsResult pool(String label, String geneLabel, List<DnDsResult> results) {
        if (results == null || results.isEmpty()) {
            throw new GviComputationException("Cannot pool dN/dS: no results supplied");
        }
        double totalS = 0, totalN = 0, totalSd = 0, totalNd = 0;
        int compared = 0, excluded = 0;
        for (DnDsResult r : results) {
            totalS += r.synonymousSites();
            totalN += r.nonsynonymousSites();
            totalSd += r.synonymousDifferences();
            totalNd += r.nonsynonymousDifferences();
            compared += r.codonsCompared();
            excluded += r.excludedCodons();
        }
        Accumulation acc = new Accumulation(totalS, totalN, totalSd, totalNd, compared, excluded);
        Omega omega = omegaFrom(acc);

        List<String> diagnostics = new ArrayList<>();
        diagnostics.add("Pooled Nei-Gojobori counts across " + results.size() + " result(s) (not an average of individual ratios)");
        if (omega.note() != null) diagnostics.add(omega.note());
        diagnostics.add(selectionTest.test(acc.totalSd(), acc.totalS(), acc.totalNd(), acc.totalN()).note());

        String category = DnDsGeneReferenceTable.genericClassification(omega.omega())
                + " (pooled across " + results.size() + " sequence(s), no single gene-specific reference band applies)";

        return new DnDsResult(label, geneLabel, omega.dN(), omega.dS(), omega.omega(), omega.jcApplied(),
                totalS, totalN, totalSd, totalNd, compared, excluded, category, diagnostics);
    }

    /**
     * Site-resolved dN/dS: slides a {@code windowCodons}-wide window across
     * the gene in steps of {@code stepCodons}, pooling Nei-Gojobori counts
     * within each window independently. Windows with zero comparable codons
     * (e.g. entirely gaps) are silently skipped rather than reported as a
     * spurious 0/0 ratio.
     */
    public List<SlidingWindowDnDsResult> slidingWindow(String sequenceId, String geneName, String refCds, String queryCds,
                                                         int windowCodons, int stepCodons) {
        if (windowCodons <= 0 || stepCodons <= 0) {
            throw new GviInputException("windowCodons and stepCodons must both be positive (got window=" + windowCodons + ", step=" + stepCodons + ")");
        }
        if (refCds.length() != queryCds.length()) {
            throw new GviComputationException("Sliding-window dN/dS for '" + sequenceId + "'/'" + geneName
                    + "': reference and query CDS are different lengths");
        }
        List<String> refCodons = CodonUtil.splitIntoCodons(refCds).codons();
        List<String> qryCodons = CodonUtil.splitIntoCodons(queryCds).codons();
        int total = refCodons.size();

        List<SlidingWindowDnDsResult> windows = new ArrayList<>();
        for (int start = 0; start + windowCodons <= total; start += stepCodons) {
            int end = start + windowCodons;
            Accumulation acc = accumulate(refCodons, qryCodons, start, end);
            if (acc.compared() == 0) continue;
            Omega omega = omegaFrom(acc);
            NeiGojoboriSelectionTest.Result significance = selectionTest.test(acc.totalSd(), acc.totalS(), acc.totalNd(), acc.totalN());
            windows.add(new SlidingWindowDnDsResult(sequenceId, geneName, start + 1, end, omega.dN(), omega.dS(), omega.omega(),
                    acc.totalS(), acc.totalN(), acc.totalSd(), acc.totalNd(),
                    significance.zScore(), significance.pValueTwoTailed(), significance.verdict()));
        }
        return windows;
    }
}
