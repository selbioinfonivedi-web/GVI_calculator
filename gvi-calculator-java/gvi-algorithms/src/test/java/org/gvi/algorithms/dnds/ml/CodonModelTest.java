package org.gvi.algorithms.dnds.ml;

import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.within;

/**
 * Basic correctness properties of {@link CodonModel} (GY94), independent of
 * the ground-truth recovery check in {@link CodonMlFitterGroundTruthTest}:
 * P(t) must be a genuine stochastic matrix (rows sum to 1) at every t, must
 * approach the identity as t -&gt; 0, and must approach the stationary
 * distribution (every row converging to pi) as t -&gt; infinity -- standard
 * sanity checks for any continuous-time Markov substitution model.
 */
class CodonModelTest {

    private final double[] pi = uniformFrequencies();

    @Test
    void transitionProbabilityRowsSumToOneAtAnyBranchLength() {
        CodonModel model = CodonModel.of(pi, 2.5, 0.4);
        for (double t : new double[]{0.001, 0.05, 0.5, 2.0}) {
            double[][] p = model.transitionProbabilities(t);
            for (int i = 0; i < CodonAlphabet.SIZE; i++) {
                double rowSum = 0;
                for (int j = 0; j < CodonAlphabet.SIZE; j++) rowSum += p[i][j];
                assertThat(rowSum).as("row %d sum at t=%f", i, t).isCloseTo(1.0, within(1e-6));
            }
        }
    }

    @Test
    void approachesTheIdentityAsBranchLengthGoesToZero() {
        CodonModel model = CodonModel.of(pi, 2.5, 0.4);
        double[][] p = model.transitionProbabilities(1e-9);
        for (int i = 0; i < CodonAlphabet.SIZE; i++) {
            assertThat(p[i][i]).isCloseTo(1.0, within(1e-4));
        }
    }

    @Test
    void approachesTheStationaryDistributionAsBranchLengthGrows() {
        CodonModel model = CodonModel.of(pi, 2.5, 0.4);
        double[][] p = model.transitionProbabilities(50.0);
        for (int i = 0; i < CodonAlphabet.SIZE; i++) {
            for (int j = 0; j < CodonAlphabet.SIZE; j++) {
                assertThat(p[i][j]).as("P[%d][%d] at large t", i, j).isCloseTo(pi[j], within(1e-4));
            }
        }
    }

    @Test
    void higherOmegaMakesNonsynonymousSubstitutionsMoreLikelyOverAFixedTime() {
        CodonModel purifying = CodonModel.of(pi, 2.0, 0.1);
        CodonModel diversifying = CodonModel.of(pi, 2.0, 5.0);

        // TTT (Phe) -> TTA (Leu): a single-position, nonsynonymous change
        int from = CodonAlphabet.index("TTT");
        int to = CodonAlphabet.index("TTA");

        double t = 0.1;
        double pPurifying = purifying.transitionProbabilities(t)[from][to];
        double pDiversifying = diversifying.transitionProbabilities(t)[from][to];

        assertThat(pDiversifying).isGreaterThan(pPurifying);
    }

    private double[] uniformFrequencies() {
        // Reuse F3x4 with a fixed-seed synthetic sample rather than a literal uniform vector,
        // so the model is exercised under realistic (non-degenerate) frequencies.
        Random random = new Random(42L);
        List<String> sample = new ArrayList<>();
        char[] bases = {'A', 'C', 'G', 'T'};
        while (sample.size() < 3000) {
            StringBuilder codon = new StringBuilder(3);
            for (int pos = 0; pos < 3; pos++) codon.append(bases[random.nextInt(4)]);
            String c = codon.toString();
            if (CodonAlphabet.isSenseCodon(c)) sample.add(c);
        }
        return CodonFrequencies.f3x4(sample);
    }
}
