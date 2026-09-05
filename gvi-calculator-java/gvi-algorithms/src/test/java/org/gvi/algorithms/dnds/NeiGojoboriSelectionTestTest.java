package org.gvi.algorithms.dnds;

import org.junit.jupiter.api.Test;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.within;

class NeiGojoboriSelectionTestTest {

    private final NeiGojoboriSelectionTest test = new NeiGojoboriSelectionTest();

    /**
     * Hand-computed: synonymousSites=20, synonymousDifferences=2 (pS=0.1),
     * nonsynonymousSites=80, nonsynonymousDifferences=40 (pN=0.5) -- strong,
     * clearly significant dN >> dS signal. Z and its variance components are
     * independently recomputed here from the textbook formulas (Nei & Kumar
     * 2000) rather than trusting the implementation's own internals.
     */
    @Test
    void detectsSignificantPositiveSelectionAndMatchesHandComputedZ() {
        double sSites = 20, sDiffs = 2, nSites = 80, nDiffs = 40;

        NeiGojoboriSelectionTest.Result result = test.test(sDiffs, sSites, nDiffs, nSites);

        double expectedZ = handComputedZ(sDiffs, sSites, nDiffs, nSites);
        assertThat(result.zScore()).isCloseTo(expectedZ, within(1e-6));
        assertThat(result.zScore()).isGreaterThan(3.0); // sanity: this is a strong, obvious signal
        assertThat(result.pValueTwoTailed()).isLessThan(0.05);
        assertThat(result.verdict()).isEqualTo(NeiGojoboriSelectionTest.Verdict.POSITIVE_SELECTION);
    }

    @Test
    void detectsSignificantPurifyingSelectionWhenRolesAreReversed() {
        // pN small, pS large: the mirror image of the positive-selection case above.
        NeiGojoboriSelectionTest.Result result = test.test(40, 80, 2, 20);

        assertThat(result.zScore()).isLessThan(-3.0);
        assertThat(result.pValueTwoTailed()).isLessThan(0.05);
        assertThat(result.verdict()).isEqualTo(NeiGojoboriSelectionTest.Verdict.PURIFYING_SELECTION);
    }

    @Test
    void reportsNotSignificantWhenDnAndDsAreEssentiallyEqual() {
        // Small counts, pS and pN nearly identical -- no real signal either direction.
        NeiGojoboriSelectionTest.Result result = test.test(5, 50, 5, 50);

        assertThat(Math.abs(result.zScore())).isLessThan(1.0);
        assertThat(result.pValueTwoTailed()).isGreaterThan(0.05);
        assertThat(result.verdict()).isEqualTo(NeiGojoboriSelectionTest.Verdict.NOT_SIGNIFICANT);
    }

    /**
     * Regression test for the Lumpy Skin Disease false positive.
     * <p>
     * These are the real pooled Nei-Gojobori counts from this project's LSD alignment: 8 near-identical
     * poxvirus sequences, ~773 synonymous sites, and ZERO observed synonymous differences. The JC variance
     * p(1-p)/[n(1-4p/3)^2] is exactly 0 at p=0, so the zero dropped out of sqrt(varN + varS) and the test
     * returned Z=2.998, p=0.0027, POSITIVE_SELECTION -- from a dataset averaging about one mutation per
     * sequence. That verdict satisfied QualityGate, letting a continuity-corrected omega of 4.86 saturate
     * the composite and supply 35% of the headline GVI.
     * <p>
     * A zero count carries no information about dS; it does not pin dS down with infinite precision.
     */
    @Test
    void refusesToTestWhenNoSynonymousDifferencesWereObserved() {
        NeiGojoboriSelectionTest.Result result = test.test(0, 773.33, 12, 2320);

        assertThat(result.verdict()).isEqualTo(NeiGojoboriSelectionTest.Verdict.NOT_SIGNIFICANT);
        assertThat(result.zScore()).isNaN();
        assertThat(result.pValueTwoTailed()).isNaN();
        assertThat(result.note()).contains("No synonymous differences were observed");
    }

    /** The same zero-variance degeneracy mirrored: pN=0 would manufacture spurious *purifying* significance. */
    @Test
    void refusesToTestWhenNoNonsynonymousDifferencesWereObserved() {
        NeiGojoboriSelectionTest.Result result = test.test(12, 773.33, 0, 2320);

        assertThat(result.verdict()).isEqualTo(NeiGojoboriSelectionTest.Verdict.NOT_SIGNIFICANT);
        assertThat(result.zScore()).isNaN();
        assertThat(result.note()).contains("No nonsynonymous differences were observed");
    }

    /**
     * The guard must be narrow: one more observed synonymous difference is still very little evidence, but
     * the variance is no longer degenerate and the test is legitimately applicable again. Guarding on
     * "small" rather than "zero" would silently discard real low-divergence comparisons.
     */
    @Test
    void stillTestsWhenASingleSynonymousDifferenceIsObserved() {
        NeiGojoboriSelectionTest.Result result = test.test(1, 773.33, 12, 2320);

        assertThat(result.zScore()).isNotNaN();
        assertThat(result.pValueTwoTailed()).isBetween(0.0, 1.0);
    }

    @Test
    void isInconclusiveWhenSiteCountsAreZero() {
        NeiGojoboriSelectionTest.Result result = test.test(0, 0, 3, 10);

        assertThat(result.zScore()).isNaN();
        assertThat(result.verdict()).isEqualTo(NeiGojoboriSelectionTest.Verdict.NOT_SIGNIFICANT);
        assertThat(result.note()).contains("inconclusive");
    }

    @Test
    void acceptsAnyDnDsCountsImplementationIncludingSlidingWindowResults() {
        SlidingWindowDnDsResult window = new SlidingWindowDnDsResult("s1", "geneX", 1, 10, 0.8, 0.1, 8.0, 20, 80, 2, 40,
                0.0, 1.0, NeiGojoboriSelectionTest.Verdict.NOT_SIGNIFICANT);

        NeiGojoboriSelectionTest.Result result = test.test(window);

        assertThat(result.verdict()).isEqualTo(NeiGojoboriSelectionTest.Verdict.POSITIVE_SELECTION);
    }

    private double handComputedZ(double sDiffs, double sSites, double nDiffs, double nSites) {
        double pS = sDiffs / sSites;
        double pN = nDiffs / nSites;
        double dS = jc(pS);
        double dN = jc(pN);
        double varS = jcVariance(pS, sSites);
        double varN = jcVariance(pN, nSites);
        return (dN - dS) / Math.sqrt(varN + varS);
    }

    private double jc(double p) {
        return -0.75 * Math.log(1.0 - (4.0 / 3.0) * p);
    }

    private double jcVariance(double p, double n) {
        double arg = 1.0 - (4.0 / 3.0) * p;
        return p * (1 - p) / (n * arg * arg);
    }
}
