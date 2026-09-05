package org.gvi.algorithms.dnds;

import org.apache.commons.math3.distribution.NormalDistribution;

/**
 * The Nei-Gojobori codon-based Z-test of selection (Nei &amp; Gojobori 1986;
 * variance formula from Nei &amp; Kumar 2000, "Molecular Evolution and
 * Phylogenetics") -- the standard way to ask, within the Nei-Gojobori
 * counting framework, whether one specific gene, pooled dataset, or sliding
 * window's dN is <em>statistically significantly</em> different from its
 * dS, rather than just reporting a ratio and leaving significance to eyeballing.
 * This is the same test MEGA's "Codon-Based Test of Selection" menu runs.
 * <p>
 * It is deliberately <b>not</b> a maximum-likelihood codon site-model
 * (PAML's M1a/M2a, HyPhy's FEL/MEME, branch-site tests) -- those model an
 * entire alignment jointly with per-site omega mixture classes and
 * empirical-Bayes site classification, a substantially larger undertaking
 * this project has not built. This is a real, named, standard significance
 * test layered on the counting method already implemented here, scoped to
 * answer "is this window's dN vs dS difference likely real" rather than
 * "which exact class does this site belong to".
 * <p>
 * Test statistic: Z = (dN - dS) / sqrt(Var(dN) + Var(dS)), where dN/dS are
 * the Jukes-Cantor-corrected distances and Var(d) = p(1-p) / [n*(1-4p/3)^2]
 * is the standard large-sample variance of the JC-corrected distance
 * estimator (p = uncorrected proportion of differences, n = site count).
 */
public final class NeiGojoboriSelectionTest {

    private static final NormalDistribution STANDARD_NORMAL = new NormalDistribution(0, 1);
    private static final double SIGNIFICANCE_LEVEL = 0.05;

    public enum Verdict {POSITIVE_SELECTION, PURIFYING_SELECTION, NOT_SIGNIFICANT}

    public record Result(double zScore, double pValueTwoTailed, Verdict verdict, String note) {
    }

    public Result test(DnDsCounts counts) {
        return test(counts.synonymousDifferences(), counts.synonymousSites(),
                counts.nonsynonymousDifferences(), counts.nonsynonymousSites());
    }

    public Result test(double synonymousDifferences, double synonymousSites,
                        double nonsynonymousDifferences, double nonsynonymousSites) {
        if (synonymousSites <= 0 || nonsynonymousSites <= 0) {
            return inconclusive("zero synonymous or nonsynonymous sites available");
        }
        // A boundary count of zero breaks the large-sample normal approximation this test rests on.
        // Var(d) = p(1-p)/[n(1-4p/3)^2] is the binomial variance of the JC-corrected distance, and it
        // collapses to EXACTLY 0 at p=0 -- so an observed count of zero is treated as "dS is known with
        // infinite precision to be exactly 0" when it actually means "we have no information about dS".
        // The zero then vanishes from the denominator sqrt(varN + varS), and Z = dN/sqrt(varN) comes back
        // large and "significant" no matter how little divergence was actually observed.
        //
        // This was not hypothetical. On this project's Lumpy Skin Disease alignment -- 8 near-identical
        // poxvirus sequences averaging ~1 mutation each -- zero observed synonymous differences produced
        // Z=2.998, p=0.0027, "POSITIVE_SELECTION". That verdict then satisfied QualityGate's
        // significance check, so a continuity-corrected omega of 4.86 saturated the composite's 0-3
        // normalization at 1.0 and contributed 0.142 of a 0.401 headline GVI: 35% of the score, from an
        // artifact of dividing by zero variance. Zero synonymous change on near-identical sequences is
        // the expected outcome of almost no divergence (nonsynonymous sites outnumber synonymous ~3:1),
        // not evidence of diversifying selection.
        //
        // Both directions are guarded: pN=0 is the same degeneracy mirrored, and would manufacture
        // spurious *purifying* significance the same way.
        if (synonymousDifferences <= 0 || nonsynonymousDifferences <= 0) {
            return inconclusive(String.format(
                    "%s differences were observed (Sd=%.1f in %.1f synonymous sites, Nd=%.1f in %.1f nonsynonymous sites). "
                            + "The Jukes-Cantor variance of a zero proportion is exactly zero, which would give that distance "
                            + "infinite apparent precision and manufacture a significant Z from no evidence. A zero count means "
                            + "no information, not perfect information -- the normal approximation does not hold at the boundary, "
                            + "so no selection verdict is issued. Add more divergent sequences or a longer coding region.",
                    synonymousDifferences <= 0 && nonsynonymousDifferences <= 0 ? "No synonymous and no nonsynonymous"
                            : synonymousDifferences <= 0 ? "No synonymous" : "No nonsynonymous",
                    synonymousDifferences, synonymousSites, nonsynonymousDifferences, nonsynonymousSites));
        }

        double pS = synonymousDifferences / synonymousSites;
        double pN = nonsynonymousDifferences / nonsynonymousSites;

        double varS = jcVariance(pS, synonymousSites);
        double varN = jcVariance(pN, nonsynonymousSites);
        if (Double.isNaN(varS) || Double.isNaN(varN)) {
            return inconclusive("Jukes-Cantor correction saturated (pS or pN too high for a reliable distance estimate)");
        }
        if (varS + varN <= 0.0) {
            return inconclusive("no variation in observed differences to test against");
        }

        double dS = jcCorrect(pS);
        double dN = jcCorrect(pN);
        double z = (dN - dS) / Math.sqrt(varN + varS);
        double pTwoTailed = 2.0 * (1.0 - STANDARD_NORMAL.cumulativeProbability(Math.abs(z)));

        if (pTwoTailed < SIGNIFICANCE_LEVEL && z > 0) {
            return new Result(z, pTwoTailed, Verdict.POSITIVE_SELECTION, String.format(
                    "dN significantly exceeds dS (Nei-Gojobori Z-test, Z=%.3f, two-tailed p=%.4f) -- positive/diversifying selection", z, pTwoTailed));
        }
        if (pTwoTailed < SIGNIFICANCE_LEVEL) {
            return new Result(z, pTwoTailed, Verdict.PURIFYING_SELECTION, String.format(
                    "dS significantly exceeds dN (Nei-Gojobori Z-test, Z=%.3f, two-tailed p=%.4f) -- purifying/negative selection", z, pTwoTailed));
        }
        return new Result(z, pTwoTailed, Verdict.NOT_SIGNIFICANT, String.format(
                "dN and dS not significantly different (Nei-Gojobori Z-test, Z=%.3f, two-tailed p=%.4f) -- neutral evolution cannot be rejected", z, pTwoTailed));
    }

    private Result inconclusive(String reason) {
        return new Result(Double.NaN, Double.NaN, Verdict.NOT_SIGNIFICANT, "Selection test inconclusive: " + reason);
    }

    private double jcVariance(double p, double n) {
        double arg = 1.0 - (4.0 / 3.0) * p;
        if (arg <= 0.0) return Double.NaN;
        return p * (1.0 - p) / (n * arg * arg);
    }

    private double jcCorrect(double p) {
        return -0.75 * Math.log(1.0 - (4.0 / 3.0) * p);
    }
}
