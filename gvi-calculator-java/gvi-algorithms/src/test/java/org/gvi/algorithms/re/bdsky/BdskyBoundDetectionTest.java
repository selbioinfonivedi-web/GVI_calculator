package org.gvi.algorithms.re.bdsky;

import org.junit.jupiter.api.Test;

import static org.assertj.core.api.Assertions.assertThat;

/**
 * A birth-death fit that comes to rest against its own search bound is not an estimate of Re --
 * the likelihood was still climbing when the optimizer ran out of interval, so the number
 * reported describes the bracket rather than the data. Reporting it is worse than reporting
 * nothing: an Re of 20 reads as an explosive epidemic when it actually means "not identifiable".
 */
class BdskyBoundDetectionTest {

    @Test
    void detectsTheRealPinnedValueSeenOnActualData() {
        // Trypanosomiasis fitted to 19.99997674 against a bound of 20.0 -- printed as "20.0000".
        // Brent stops NEAR a bound rather than exactly on it, which is why the tolerance is
        // interval-relative rather than machine epsilon.
        assertThat(BdskyMlFitter.reIsAtSearchBound(19.99997673997714)).isTrue();
    }

    @Test
    void detectsBothBounds() {
        assertThat(BdskyMlFitter.reIsAtSearchBound(BdskyMlFitter.RE_MAX)).isTrue();
        assertThat(BdskyMlFitter.reIsAtSearchBound(BdskyMlFitter.RE_MIN)).isTrue();
    }

    @Test
    void acceptsGenuineInteriorOptima() {
        // The real fits from this project's corpus must all still be treated as estimates.
        for (double re : new double[]{1.0045, 1.0180, 1.0917, 1.3227, 1.5294, 2.5, 10.0}) {
            assertThat(BdskyMlFitter.reIsAtSearchBound(re)).as("Re=%s is an interior optimum", re).isFalse();
        }
    }

    @Test
    void toleranceIsNarrowEnoughToNotSwallowPlausibleValues() {
        // A tolerance loose enough to reject Re = 19.0 would be discarding real (if extreme) fits.
        assertThat(BdskyMlFitter.reIsAtSearchBound(19.0)).isFalse();
    }
}
