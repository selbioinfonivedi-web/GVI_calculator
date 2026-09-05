package org.gvi.algorithms.cai;

import org.apache.commons.math3.distribution.HypergeometricDistribution;
import org.gvi.core.exception.GviComputationException;
import org.gvi.core.model.NucleotideSequence;
import org.gvi.core.util.IupacUtil;

import java.util.ArrayList;
import java.util.List;
import java.util.function.Predicate;

/**
 * Detects whether the substitutions between a reference and query genome
 * are consistent with known host-immune RNA-editing mutational pressure,
 * rather than being an unexplained compositional anomaly (GC_Deviation's
 * default "HGT suspected" reading, Section 5.8/8.4):
 * <p>
 * - <b>APOBEC3-family cytidine deaminases</b> hyper-edit C-&gt;U (C-&gt;T in
 * genome notation) at a 5'-U-C ("TCW") dinucleotide hotspot -- the same
 * mechanism behind COSMIC's SBS2/SBS13 mutational signatures in human
 * cancer genomics, and documented directly in coronavirus genomes
 * (Simmonds 2020, "Rampant C-&gt;U hypermutation in the genomes of
 * SARS-CoV-2 and other coronaviruses").
 * <p>
 * - <b>ADAR1 adenosine deaminases</b> hyper-edit A-&gt;I (read as A-&gt;G)
 * preferring a 5'-U or 5'-A neighbor -- documented in measles SSPE genomes
 * and other persistent RNA virus infections.
 * <p>
 * Method: for each substitution type, split the reference sites of the
 * source base into "hotspot context" (5' neighbor matches the editor's
 * known preference) vs "other context", and run a one-sided Fisher exact
 * test (via the hypergeometric distribution) for enrichment of the
 * substitution in the hotspot context, plus a Haldane-Anscombe-corrected
 * odds ratio. This is a real, literature-standard mutational-signature
 * detection method, deliberately scoped down to the two mechanisms with
 * the strongest, best-documented single-dinucleotide-context signal in
 * viral genomes -- not a full 96-trinucleotide-channel signature
 * deconvolution (which needs far more mutations than one genome pair
 * provides to be statistically identifiable).
 */
public final class MutationalSignatureAnalyzer {

    private static final double SIGNIFICANCE_LEVEL = 0.05;
    /** Below this many sites in either context group, the Fisher test is underpowered and reported as inconclusive rather than "not significant". */
    private static final int MIN_SITES_PER_CONTEXT = 10;

    public MutationalSignatureResult analyze(NucleotideSequence reference, NucleotideSequence query) {
        String ref = reference.getSequence();
        String qry = query.getSequence();
        if (ref.length() != qry.length()) {
            throw new GviComputationException("Cannot scan for mutational signatures: '" + reference.getId()
                    + "' and '" + query.getId() + "' are not aligned to equal length");
        }

        MutationalSignatureResult.SignatureTest apobec = scan(ref, qry, 'C', 'T',
                prev -> prev == 'T', "APOBEC3-like C->T (5'-U/T hotspot context)");
        MutationalSignatureResult.SignatureTest adar = scan(ref, qry, 'A', 'G',
                prev -> prev == 'T' || prev == 'A', "ADAR-like A->G (5'-U/A hotspot context)");

        List<String> diagnostics = new ArrayList<>();
        diagnostics.add(describe(apobec));
        diagnostics.add(describe(adar));

        return new MutationalSignatureResult(apobec, adar, diagnostics);
    }

    private MutationalSignatureResult.SignatureTest scan(String ref, String qry, char fromBase, char toBase,
                                                           Predicate<Character> hotspotContext, String label) {
        int hotspotMutated = 0, hotspotTotal = 0, otherMutated = 0, otherTotal = 0;
        for (int i = 1; i < ref.length(); i++) { // start at 1: every site needs a preceding reference base for context
            char r = Character.toUpperCase(ref.charAt(i));
            char prev = Character.toUpperCase(ref.charAt(i - 1));
            if (r != fromBase || !IupacUtil.isUnambiguousBase(prev)) continue;

            char q = Character.toUpperCase(qry.charAt(i));
            if (!IupacUtil.isComparable(r) || !IupacUtil.isComparable(q)) continue; // gaps/ambiguous excluded, same convention as GD/MB

            boolean mutatedToTarget = (q == toBase);
            if (hotspotContext.test(prev)) {
                hotspotTotal++;
                if (mutatedToTarget) hotspotMutated++;
            } else {
                otherTotal++;
                if (mutatedToTarget) otherMutated++;
            }
        }

        if (hotspotTotal < MIN_SITES_PER_CONTEXT || otherTotal < MIN_SITES_PER_CONTEXT) {
            return MutationalSignatureResult.SignatureTest.inconclusive(label, hotspotMutated, hotspotTotal, otherMutated, otherTotal);
        }

        double[] oddsRatioAndPValue = fisherExactEnrichment(hotspotMutated, hotspotTotal, otherMutated, otherTotal);
        double oddsRatio = oddsRatioAndPValue[0];
        double pValue = oddsRatioAndPValue[1];
        boolean significant = pValue < SIGNIFICANCE_LEVEL && oddsRatio > 1.0;
        return new MutationalSignatureResult.SignatureTest(label, hotspotMutated, hotspotTotal, otherMutated, otherTotal,
                oddsRatio, pValue, significant);
    }

    /**
     * One-sided Fisher exact test (via the hypergeometric distribution) for
     * enrichment of "mutated" counts in the hotspot context versus the
     * other context, plus a Haldane-Anscombe-corrected odds ratio (the same
     * zero-cell correction already used for dS=0 in {@code DnDsCalculator}).
     * Returns {oddsRatio, oneSidedPValue}.
     */
    private double[] fisherExactEnrichment(int hotspotMutated, int hotspotTotal, int otherMutated, int otherTotal) {
        int a = hotspotMutated, b = hotspotTotal - hotspotMutated;
        int c = otherMutated, d = otherTotal - otherMutated;

        double aC = a, bC = b, cC = c, dC = d;
        if (a == 0 || b == 0 || c == 0 || d == 0) {
            aC += 0.5;
            bC += 0.5;
            cC += 0.5;
            dC += 0.5;
        }
        double oddsRatio = (aC * dC) / (bC * cC);

        int population = hotspotTotal + otherTotal;
        int successesInPopulation = a + c;
        int sampleSize = hotspotTotal;
        if (successesInPopulation == 0 || successesInPopulation == population) {
            return new double[]{oddsRatio, 1.0}; // no variation between contexts to test
        }
        HypergeometricDistribution dist = new HypergeometricDistribution(population, successesInPopulation, sampleSize);
        double pValue = a == 0 ? 1.0 : dist.upperCumulativeProbability(a);
        return new double[]{oddsRatio, pValue};
    }

    private String describe(MutationalSignatureResult.SignatureTest t) {
        if (Double.isNaN(t.pValue())) {
            return t.label() + ": inconclusive (too few sites -- hotspot context=" + t.hotspotContextTotal()
                    + ", other context=" + t.otherContextTotal() + ", need >=" + MIN_SITES_PER_CONTEXT + " each)";
        }
        String hotspotRate = String.format("%.1f%% (%d/%d)", 100.0 * t.hotspotContextMutated() / t.hotspotContextTotal(),
                t.hotspotContextMutated(), t.hotspotContextTotal());
        String otherRate = String.format("%.1f%% (%d/%d)", 100.0 * t.otherContextMutated() / t.otherContextTotal(),
                t.otherContextMutated(), t.otherContextTotal());
        return String.format("%s: hotspot-context rate=%s vs other-context rate=%s, odds ratio=%.2f, one-sided p=%.4f%s",
                t.label(), hotspotRate, otherRate, t.oddsRatio(), t.pValue(),
                t.significant() ? " -- SIGNIFICANT enrichment" : "");
    }
}
