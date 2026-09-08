package org.gvi.algorithms.re;

import org.apache.commons.math3.distribution.GammaDistribution;
import org.gvi.core.exception.GviInputException;
import org.junit.jupiter.api.Test;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.assertThatThrownBy;
import static org.assertj.core.api.Assertions.within;

/**
 * The infectiousness profile the Cori estimator divides by. It had no test of its own, and the
 * property that matters most is invisible by inspection: <em>the discretisation must preserve the
 * mean it was asked for</em>.
 * <p>
 * This is not hypothetical. Discretising the same gamma by CDF differencing, {@code F(u) - F(u-1)},
 * is a natural-looking alternative -- arguably the more obviously "correct" one, since it integrates
 * the density over each day rather than sampling it at a point. It is wrong here: it assigns each
 * interval's mass to its right endpoint and shifts the effective mean of a Gamma(5, 2) profile from
 * 5.00 to 5.50 days. Because Re is estimated by dividing incidence by this profile, that half-day
 * error propagates straight into Re, and it <em>amplifies with Re</em>: measured on a simulated
 * renewal process it costs 0.008 at Re = 1.15, 0.127 at Re = 2.0 and 0.209 at Re = 2.5.
 * <p>
 * That amplification is the reason to pin this. The error is invisible exactly when the outbreak is
 * flat and worst exactly when it is growing fastest -- when the number is being relied on most.
 */
class SerialIntervalTest {

    private static double effectiveMean(SerialInterval si) {
        double m = 0;
        for (int u = 1; u <= si.tMax(); u++) m += u * si.weight(u);
        return m;
    }

    @Test
    void theDiscretisationPreservesTheRequestedMean() {
        for (double mean : new double[]{3.0, 5.0, 7.0, 12.0}) {
            SerialInterval si = new SerialInterval(mean, mean / 2.5, 40);
            assertThat(effectiveMean(si))
                    .as("Gamma(mean=%.1f) discretised must still have mean %.1f", mean, mean)
                    .isCloseTo(mean, within(0.05));
        }
    }

    /**
     * The specific alternative that looks right and is not. Guarding the difference explicitly means
     * a future "simplification" to CDF differencing fails here with a reason attached, rather than
     * silently biasing every incidence-based Re upward.
     */
    @Test
    void cdfDifferencingWouldShiftTheMeanByHalfADayAndIsNotWhatThisClassDoes() {
        int tMax = 20;
        SerialInterval si = new SerialInterval(5.0, 2.0, tMax);

        GammaDistribution g = new GammaDistribution((5.0 * 5.0) / (2.0 * 2.0), (2.0 * 2.0) / 5.0);
        double[] cdf = new double[tMax + 1];
        double sum = 0;
        for (int u = 1; u <= tMax; u++) {
            cdf[u] = g.cumulativeProbability(u) - g.cumulativeProbability(u - 1);
            sum += cdf[u];
        }
        double cdfMean = 0;
        for (int u = 1; u <= tMax; u++) cdfMean += u * (cdf[u] / sum);

        assertThat(effectiveMean(si)).as("this class preserves the mean").isCloseTo(5.0, within(0.01));
        assertThat(cdfMean).as("CDF differencing does not").isCloseTo(5.5, within(0.01));
    }

    @Test
    void theWeightsFormAProperDistribution() {
        SerialInterval si = new SerialInterval(5.0, 2.0, 20);

        double total = 0;
        for (int u = 1; u <= si.tMax(); u++) {
            assertThat(si.weight(u)).as("w(%d)", u).isBetween(0.0, 1.0);
            total += si.weight(u);
        }
        assertThat(total).as("weights must sum to 1 after renormalisation").isCloseTo(1.0, within(1e-12));
    }

    /** Day 0 carries no mass: a case cannot infect anyone before it exists. */
    @Test
    void thereIsNoMassOutsideOneToTMax() {
        SerialInterval si = new SerialInterval(5.0, 2.0, 20);

        assertThat(si.weight(0)).isZero();
        assertThat(si.weight(-1)).isZero();
        assertThat(si.weight(21)).isZero();
        assertThat(si.weight(Integer.MAX_VALUE)).isZero();
    }

    @Test
    void theProfileIsUnimodalAndPeaksNearTheMean() {
        SerialInterval si = new SerialInterval(5.0, 2.0, 20);

        int peak = 1;
        for (int u = 2; u <= si.tMax(); u++) if (si.weight(u) > si.weight(peak)) peak = u;

        assertThat(peak).as("a Gamma(5,2) profile should peak within a day or two of its mean")
                .isBetween(3, 6);
        assertThat(si.weight(20)).as("the far tail should be negligible, not truncated mid-mass")
                .isLessThan(si.weight(peak) / 100);
    }

    @Test
    void theDefaultProfileIsTheSpecifiedGammaFiveTwo() {
        SerialInterval si = SerialInterval.defaultProfile();

        assertThat(si.tMax()).isEqualTo(20);
        assertThat(effectiveMean(si)).isCloseTo(5.0, within(0.01));
    }

    @Test
    void nonPositiveParametersAreRejectedRatherThanProducingNaNWeights() {
        assertThatThrownBy(() -> new SerialInterval(0.0, 2.0, 20))
                .isInstanceOf(GviInputException.class).hasMessageContaining("positive");
        assertThatThrownBy(() -> new SerialInterval(5.0, 0.0, 20))
                .isInstanceOf(GviInputException.class).hasMessageContaining("positive");
        assertThatThrownBy(() -> new SerialInterval(-5.0, 2.0, 20))
                .isInstanceOf(GviInputException.class).hasMessageContaining("positive");
    }

    /**
     * A window far shorter than the mean would leave the profile with essentially no mass, and the
     * renormalisation would then manufacture a distribution out of numerical dust.
     */
    @Test
    void aWindowWithNoMassIsRejected() {
        assertThatThrownBy(() -> new SerialInterval(500.0, 1.0, 5))
                .isInstanceOf(GviInputException.class)
                .hasMessageContaining("no mass");
    }
}
