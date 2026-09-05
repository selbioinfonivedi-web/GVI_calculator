package org.gvi.algorithms.phylo.model;

import org.junit.jupiter.api.Test;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.within;

class GtrModelTest {

    @Test
    void atTimeZeroTransitionMatrixIsIdentity() {
        GtrModel gtr = GtrModel.of(new double[]{0.3, 0.2, 0.2, 0.3}, new double[]{1.2, 3.1, 0.8, 0.9, 2.5, 1.0});
        double[][] p = gtr.transitionProbabilities(0.0);
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                assertThat(p[i][j]).as("P[%d][%d] at t=0", i, j).isCloseTo(i == j ? 1.0 : 0.0, within(1e-9));
            }
        }
    }

    @Test
    void everyRowIsAValidProbabilityDistributionAtAnyBranchLength() {
        GtrModel gtr = GtrModel.of(new double[]{0.15, 0.35, 0.35, 0.15}, new double[]{2.0, 4.0, 0.5, 1.0, 3.0, 0.7});
        for (double t : new double[]{0.001, 0.05, 0.3, 1.5, 10.0}) {
            double[][] p = gtr.transitionProbabilities(t);
            for (int i = 0; i < 4; i++) {
                double rowSum = 0;
                for (int j = 0; j < 4; j++) {
                    assertThat(p[i][j]).as("P[%d][%d] at t=%s must be a valid probability", i, j, t).isBetween(-1e-9, 1.0 + 1e-9);
                    rowSum += p[i][j];
                }
                assertThat(rowSum).as("row %d sums to 1 at t=%s", i, t).isCloseTo(1.0, within(1e-6));
            }
        }
    }

    @Test
    void satisfiesDetailedBalanceReversibility() {
        // pi_i * P_ij(t) == pi_j * P_ji(t) is the defining property of a time-reversible model
        double[] freq = {0.1, 0.4, 0.4, 0.1};
        GtrModel gtr = GtrModel.of(freq, new double[]{1.5, 2.5, 0.6, 1.1, 3.3, 0.9});
        double[][] p = gtr.transitionProbabilities(0.4);
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                assertThat(freq[i] * p[i][j]).as("detailed balance (%d,%d)", i, j)
                        .isCloseTo(freq[j] * p[j][i], within(1e-9));
            }
        }
    }

    @Test
    void convergesToStationaryDistributionAtLargeT() {
        double[] freq = {0.4, 0.1, 0.1, 0.4};
        GtrModel gtr = GtrModel.of(freq, new double[]{1.0, 2.0, 3.0, 0.5, 1.5, 2.5});
        double[][] p = gtr.transitionProbabilities(1000.0);
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                assertThat(p[i][j]).as("P[%d][%d] at t->inf converges to pi[%d]", i, j, j).isCloseTo(freq[j], within(1e-6));
            }
        }
    }

    /**
     * With equal frequencies and equal rates, GTR degenerates to Jukes-Cantor,
     * which has a known closed form: P_ii(t) = 1/4 + 3/4*exp(-4t/3),
     * P_ij(t) = 1/4 - 1/4*exp(-4t/3) for i != j.
     */
    @Test
    void degeneratesToJukesCantorClosedFormWithEqualFreqsAndRates() {
        GtrModel jc = GtrModel.jukesCantor();
        double t = 0.2;
        double expectedSame = 0.25 + 0.75 * Math.exp(-4.0 * t / 3.0);
        double expectedDiff = 0.25 - 0.25 * Math.exp(-4.0 * t / 3.0);

        double[][] p = jc.transitionProbabilities(t);
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                assertThat(p[i][j]).isCloseTo(i == j ? expectedSame : expectedDiff, within(1e-9));
            }
        }
    }

    @Test
    void rejectsBaseFrequenciesNotSummingToOne() {
        org.assertj.core.api.Assertions.assertThatThrownBy(() ->
                GtrModel.of(new double[]{0.3, 0.3, 0.3, 0.3}, new double[]{1, 1, 1, 1, 1, 1}))
                .isInstanceOf(org.gvi.core.exception.GviInputException.class);
    }
}
