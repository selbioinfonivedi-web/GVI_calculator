package org.gvi.cli;

import org.gvi.algorithms.mu.MuResult;
import org.gvi.composite.IndexKey;
import org.gvi.core.model.NucleotideSequence;
import org.gvi.core.model.SequenceAlignment;
import org.gvi.core.spi.IndexResult;

import java.util.ArrayList;
import java.util.EnumMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

/**
 * Decides which indices are unfit to contribute to the composite GVI, based on
 * the same data-quality checks the pipeline already performs.
 * <p>
 * Previously those checks produced advisory text only, and ran <em>after</em>
 * the composite had already been computed, so a dataset could emit "38% of the
 * alignment is gap characters" and "mutation burden implausibly high" and still
 * report dN/dS = 1.067 as positive selection and mu = 0.209 sub/site/yr -- a
 * rate roughly five orders of magnitude above anything observed for a
 * bacterium -- with both carrying full weight in a headline score.
 * <p>
 * The composite engine already knows how to drop an index and renormalize the
 * remaining weights; it simply was never told to. This class produces that
 * instruction. Excluded indices are still computed and still reported in full
 * under their own sections -- the exclusion only removes them from the weighted
 * score, so nothing is hidden from the user.
 */
public final class QualityGate {

    /** Gap fraction above which the input is more likely padded/concatenated than genuinely aligned. */
    public static final double GAP_FRACTION_LIMIT = 0.15;
    /**
     * Mutation burden as a fraction of alignment length, above which the sequences cannot plausibly be
     * the closely-related set the reference-based indices assume.
     */
    static final double MB_FRACTION_LIMIT = 0.15;
    /** Temporal-signal R² below which a molecular-clock rate is not interpretable (the pipeline's own threshold). */
    static final double MIN_CLOCK_R_SQUARED = 0.30;

    /** One index's exclusion, with the reason surfaced to the user. */
    public record Exclusion(IndexKey key, String reason) {
    }

    private final List<Exclusion> exclusions = new ArrayList<>();

    private QualityGate() {
    }

    /**
     * Evaluates the finished indices against the alignment they came from.
     * Must run BEFORE the composite is computed.
     */
    public static QualityGate evaluate(SequenceAlignment alignment,
                                       Map<IndexKey, IndexResult> population,
                                       Map<IndexKey, IndexResult> datasetIndices) {
        return evaluate(alignment, population, datasetIndices, null);
    }

    /**
     * @param dnDsSelection the Nei-Gojobori significance verdict for the pooled dN/dS, or null when
     *                      dN/dS was not computed. Passed separately because the report serializes
     *                      {@code datasetIndices} directly, so the verdict cannot ride inside it.
     */
    public static QualityGate evaluate(SequenceAlignment alignment,
                                       Map<IndexKey, IndexResult> population,
                                       Map<IndexKey, IndexResult> datasetIndices,
                                       org.gvi.algorithms.dnds.NeiGojoboriSelectionTest.Result dnDsSelection) {
        return evaluate(alignment, population, datasetIndices, dnDsSelection,
                org.gvi.algorithms.mu.GenomeType.UNSPECIFIED);
    }

    /**
     * @param genomeType the declared genome chemistry, needed to judge whether the fitted substitution rate is
     *                   physically credible at all. Without it no scale check is possible, so an unspecified
     *                   genome type disables that one check rather than guessing.
     */
    public static QualityGate evaluate(SequenceAlignment alignment,
                                       Map<IndexKey, IndexResult> population,
                                       Map<IndexKey, IndexResult> datasetIndices,
                                       org.gvi.algorithms.dnds.NeiGojoboriSelectionTest.Result dnDsSelection,
                                       org.gvi.algorithms.mu.GenomeType genomeType) {
        QualityGate gate = new QualityGate();
        gate.checkAlignmentIntegrity(alignment, datasetIndices);
        gate.checkClockSignal(population);
        gate.checkRatePlausibility(population, genomeType);
        gate.checkReEstimatorQuality(population);
        gate.checkSelectionSignificance(datasetIndices, dnDsSelection);
        return gate;
    }

    /**
     * A substitution rate can have impeccable temporal signal and still be biologically impossible. The
     * date-randomization test and R² both ask whether the collection dates explain the divergence; neither asks
     * whether the resulting rate is a rate this kind of genome can actually have.
     * <p>
     * Lumpy Skin Disease passed both and still returned 1.85e-4 sub/site/yr for a capripoxvirus -- roughly two
     * orders of magnitude above the ~1e-6 the reference table gives for that group. {@code MuReferenceTable}
     * classified it as "Implausible for a DNA genome (No known DNA virus replicates this fast)" and printed the
     * verdict, but nothing consumed it, so the rate still scored. Worse, mu normalizes against a 1e-4 DNA
     * ceiling, so an implausible rate saturates at 1.0 and becomes the largest single contributor -- the more
     * impossible the value, the more it inflated the headline GVI.
     */
    private void checkRatePlausibility(Map<IndexKey, IndexResult> population,
                                       org.gvi.algorithms.mu.GenomeType genomeType) {
        if (!(population.get(IndexKey.MU) instanceof MuResult mu)) return;
        if (isExcluded(IndexKey.MU)) return; // already gated by alignment integrity or clock signal
        if (!org.gvi.algorithms.mu.MuReferenceTable.isImplausible(mu.muSubPerSiteYear(), genomeType)) return;

        exclude(IndexKey.MU, String.format(Locale.ROOT,
                "the fitted rate is %.4g sub/site/yr, which the reference table classifies as implausible for a DNA "
                        + "genome -- no known DNA virus replicates faster than about %.0e, and the fastest credible band "
                        + "(small ssDNA, no proofreading) tops out there. Note this is NOT a weak-signal problem: the "
                        + "temporal-signal checks may well have passed, because the dates genuinely do explain the "
                        + "divergence. What fails is the magnitude. A rate this far above the ceiling usually means the "
                        + "alignment pools separate lineages or loci, or that the marker is under selection strong enough "
                        + "to break the clock assumption itself. Excluded from the composite GVI (still reported below) -- "
                        + "otherwise it would saturate mu's normalization ceiling and contribute maximum weight precisely "
                        + "because it is impossible. Verify the input before trusting this rate.",
                mu.muSubPerSiteYear(), org.gvi.algorithms.mu.MuReferenceTable.MAX_CREDIBLE_DNA_RATE));
    }

    /**
     * An omega above 1 is only evidence of positive selection when it is statistically
     * distinguishable from 1. Small codon counts make omega extremely noisy, and a shifted reading
     * frame inflates it mechanically by turning synonymous change into apparent nonsynonymous
     * change -- so an unqualified omega greater than 1 is at least as likely to be an artifact as a
     * finding.
     * <p>
     * The Nei-Gojobori Z-test was already computed and printed; nothing acted on it. A real case
     * from this project's corpus: a pooled omega of 2.536 was reported as positive selection and
     * carried full composite weight while its own Z-test could not separate it from neutrality.
     * <p>
     * Only omega above 1 is gated. An omega below 1 that is not significant is the unremarkable
     * "no detectable selection either way" case, and the purifying-selection reading it receives
     * does not overstate anything.
     */
    private void checkSelectionSignificance(Map<IndexKey, IndexResult> datasetIndices,
                                            org.gvi.algorithms.dnds.NeiGojoboriSelectionTest.Result selection) {
        if (selection == null) return;
        IndexResult dnds = datasetIndices.get(IndexKey.DNDS);
        if (dnds == null) return;
        if (isExcluded(IndexKey.DNDS)) return; // already gated by alignment integrity
        if (dnds.primaryValue() <= 1.0) return;

        if (selection.verdict() == org.gvi.algorithms.dnds.NeiGojoboriSelectionTest.Verdict.POSITIVE_SELECTION) {
            return; // omega > 1 AND significant -- a real finding, let it score
        }

        // The Z-test declines to issue a verdict at all (NaN statistic) when it would have to divide by a
        // zero-variance term -- most often because zero synonymous differences were observed, which makes
        // omega undefined rather than large. Printing "Z=NaN, p=NaN" there would read as a failed
        // computation; the actual situation is a well-understood boundary case with a specific remedy.
        String evidence = Double.isNaN(selection.zScore())
                ? "the Nei-Gojobori Z-test could not be evaluated at all (" + selection.note() + ")"
                : String.format(Locale.ROOT,
                        "the Nei-Gojobori Z-test cannot distinguish it from neutrality (Z=%.3f, two-tailed p=%.4g, verdict %s)",
                        selection.zScore(), selection.pValueTwoTailed(), selection.verdict());

        exclude(IndexKey.DNDS, String.format(Locale.ROOT,
                "the pooled dN/dS is %.3f but %s, so it is not evidence of positive selection. An "
                        + "unqualified ratio above 1 on few codons is at least as likely to be a frameshift or "
                        + "misalignment artifact -- a shifted reading frame turns synonymous change into apparent "
                        + "nonsynonymous change. Excluded from the composite GVI (still reported below); confirm the "
                        + "reading frame and gene boundaries (--gff) before treating this as selection.",
                dnds.primaryValue(), evidence));
    }

    /**
     * The lineages-through-time fallback fits a regression slope to a LTT curve. On this project's
     * corpus it returned Re in 1.0007-1.0262 for every pathogen regardless of biology, and barely
     * responded to the generation time it is nominally driven by (varying T from 5 to 365 days moved
     * FMD's Re only from 1.0016 to 1.0056). A quantity that cannot distinguish a fast-spreading RNA
     * virus from a soil bacterium is not carrying information, and Re holds the single largest
     * configured weight -- so when Re came from that fallback rather than a real birth-death fit or
     * genuine case-incidence data, it is reported but not scored.
     */
    private void checkReEstimatorQuality(Map<IndexKey, IndexResult> population) {
        if (!(population.get(IndexKey.RE) instanceof org.gvi.algorithms.re.ReResult re)) return;
        if (re.method() != org.gvi.algorithms.re.ReMethod.PHYLODYNAMIC_FALLBACK) return;

        exclude(IndexKey.RE, "Re came from the lineages-through-time fallback, which on real data returns ~1.0 for every "
                + "pathogen and is largely insensitive to the generation time driving it -- it cannot distinguish growth "
                + "rates well enough to carry the largest weight in a composite score. Excluded from the composite GVI "
                + "(still reported below). Supply --incidence for a Cori et al. estimate, or use a dataset where the "
                + "birth-death (BDSKY) fit succeeds.");
    }

    /**
     * A misaligned or wrongly-pooled alignment invalidates every index derived from
     * position-wise comparison: the "differences" being counted are not homologous.
     * That is all of MB, GD, dN/dS, pi and mu -- so this check gates them together
     * rather than one at a time.
     */
    private void checkAlignmentIntegrity(SequenceAlignment alignment, Map<IndexKey, IndexResult> datasetIndices) {
        double gapFraction = gapFraction(alignment);
        boolean gapFail = gapFraction >= GAP_FRACTION_LIMIT;

        Double mbFraction = null;
        IndexResult mb = datasetIndices.get(IndexKey.MB);
        if (mb != null && alignment.length() > 0) {
            mbFraction = mb.primaryValue() / alignment.length();
        }
        boolean mbFail = mbFraction != null && mbFraction >= MB_FRACTION_LIMIT;

        if (!gapFail && !mbFail) return;

        String cause;
        if (gapFail && mbFail) {
            cause = String.format(Locale.ROOT,
                    "the alignment is %.0f%% gap characters and the mean mutation burden is %.0f%% of its length",
                    gapFraction * 100, mbFraction * 100);
        } else if (gapFail) {
            cause = String.format(Locale.ROOT, "the alignment is %.0f%% gap characters", gapFraction * 100);
        } else {
            cause = String.format(Locale.ROOT,
                    "the mean mutation burden is %.0f%% of the alignment length", mbFraction * 100);
        }

        String reason = cause + ", so the sequences are not the closely-related, genuinely-aligned set these "
                + "position-wise indices assume. Excluded from the composite GVI (still reported in full below). "
                + "Re-align with e.g. \"mafft --auto input.fasta > aligned.fasta\", or split genuinely distinct "
                + "lineages into separate alignments, then re-run.";

        for (IndexKey key : new IndexKey[]{IndexKey.MB, IndexKey.GD, IndexKey.DNDS, IndexKey.PI, IndexKey.MU}) {
            exclude(key, reason);
        }
    }

    /**
     * A clock rate is only interpretable when the sampling dates actually explain the
     * divergence. Below the pipeline's own R² threshold -- or at a negative rate, which
     * is not a physically meaningful quantity at all -- mu is noise, and anything derived
     * from it inherits that.
     */
    private void checkClockSignal(Map<IndexKey, IndexResult> population) {
        if (!(population.get(IndexKey.MU) instanceof MuResult mu)) return;
        if (isExcluded(IndexKey.MU)) return; // already gated by alignment integrity

        if (!mu.temporalSignalConfirmed()) {
            exclude(IndexKey.MU, "the date-randomization test failed -- randomly reshuffled collection dates reproduce "
                    + "this rate, so it reflects the alignment's divergence structure rather than elapsed time (usually "
                    + "several independently-introduced lineages pooled together). A high R² does not rescue this: R² "
                    + "measures how tightly the points sit on a line, not whether the dates are what put them there. "
                    + "Excluded from the composite GVI (still reported below); split by lineage and re-run.");
        } else if (mu.muSubPerSiteYear() <= 0) {
            exclude(IndexKey.MU, String.format(Locale.ROOT,
                    "the fitted substitution rate is %.6g sub/site/yr -- a non-positive rate is not a physically "
                            + "meaningful quantity, it means the regression found no temporal signal. Excluded from the "
                            + "composite GVI (still reported below).", mu.muSubPerSiteYear()));
        } else if (mu.rSquared() < MIN_CLOCK_R_SQUARED) {
            exclude(IndexKey.MU, String.format(Locale.ROOT,
                    "temporal signal is too weak to interpret a molecular clock (R²=%.2f, threshold %.2f): the collection "
                            + "dates explain almost none of the divergence, usually because the sampling window is too short "
                            + "or the alignment pools independently-introduced lineages. Excluded from the composite GVI "
                            + "(still reported below).", mu.rSquared(), MIN_CLOCK_R_SQUARED));
        }
    }

    private void exclude(IndexKey key, String reason) {
        if (isExcluded(key)) return;
        exclusions.add(new Exclusion(key, reason));
    }

    public boolean isExcluded(IndexKey key) {
        return exclusions.stream().anyMatch(e -> e.key() == key);
    }

    public List<Exclusion> exclusions() {
        return List.copyOf(exclusions);
    }

    /** The indices that may enter the weighted composite, with the failures removed. */
    public Map<IndexKey, IndexResult> filter(Map<IndexKey, IndexResult> available) {
        Map<IndexKey, IndexResult> kept = new EnumMap<>(IndexKey.class);
        available.forEach((k, v) -> {
            if (!isExcluded(k)) kept.put(k, v);
        });
        return kept;
    }

    public static double gapFraction(SequenceAlignment alignment) {
        long gaps = 0;
        long total = 0;
        for (NucleotideSequence seq : alignment.getSequences()) {
            String s = seq.getSequence();
            total += s.length();
            for (int i = 0; i < s.length(); i++) {
                char c = s.charAt(i);
                if (c == '-' || c == '.' || c == '?') gaps++;
            }
        }
        return total > 0 ? (double) gaps / total : 0.0;
    }
}
