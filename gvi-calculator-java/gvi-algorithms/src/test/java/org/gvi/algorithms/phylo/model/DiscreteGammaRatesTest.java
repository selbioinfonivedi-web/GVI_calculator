package org.gvi.algorithms.phylo.model;

import org.junit.jupiter.api.Test;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.within;

class DiscreteGammaRatesTest {

    @Test
    void categoryRatesAverageToOneRegardlessOfShape() {
        for (double alpha : new double[]{0.1, 0.5, 1.0, 2.0, 10.0}) {
            double[] rates = DiscreteGammaRates.categories(alpha, 4);
            double mean = 0;
            for (double r : rates) mean += r;
            mean /= rates.length;
            assertThat(mean).as("alpha=%s", alpha).isCloseTo(1.0, within(1e-4));
        }
    }

    @Test
    void categoriesAreMonotonicallyIncreasing() {
        double[] rates = DiscreteGammaRates.categories(0.5, 4);
        for (int i = 1; i < rates.length; i++) {
            assertThat(rates[i]).isGreaterThan(rates[i - 1]);
        }
    }

    @Test
    void allRatesArePositive() {
        double[] rates = DiscreteGammaRates.categories(0.3, 6);
        for (double r : rates) assertThat(r).isPositive();
    }

    @Test
    void singleCategoryIsTrivialUnitRate() {
        assertThat(DiscreteGammaRates.categories(1.0, 1)).containsExactly(1.0);
    }

    @Test
    void veryLargeAlphaApproachesNoHeterogeneity() {
        // as alpha -> infinity, Gamma(alpha,1/alpha) concentrates at 1 -- all categories should collapse toward rate 1
        double[] rates = DiscreteGammaRates.categories(100_000, 4);
        for (double r : rates) {
            assertThat(r).isCloseTo(1.0, within(0.01));
        }
    }

    @Test
    void smallAlphaGivesHighVarianceAcrossCategories() {
        // low alpha = strong rate heterogeneity: the slowest category should be much slower than the fastest
        double[] rates = DiscreteGammaRates.categories(0.1, 4);
        assertThat(rates[0]).isLessThan(0.1);
        assertThat(rates[3]).isGreaterThan(2.0);
    }

    @Test
    void rejectsNonPositiveAlpha() {
        org.assertj.core.api.Assertions.assertThatThrownBy(() -> DiscreteGammaRates.categories(0.0, 4))
                .isInstanceOf(org.gvi.core.exception.GviInputException.class);
    }
}
