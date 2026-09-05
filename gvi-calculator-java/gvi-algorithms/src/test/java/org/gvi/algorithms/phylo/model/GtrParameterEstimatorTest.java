package org.gvi.algorithms.phylo.model;

import org.gvi.core.model.NucleotideSequence;
import org.gvi.core.model.SequenceAlignment;
import org.junit.jupiter.api.Test;

import java.util.List;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.within;

class GtrParameterEstimatorTest {

    /**
     * seq1="AAAA", seq2="ACGT": differences at pos1 (A,C), pos2 (A,G), pos3 (A,T); pos0 identical.
     * Base counts across both sequences (8 bases total): A=5, C=1, G=1, T=1.
     * Hand-computed expected frequencies and rates (with the +0.5 pseudocount) below.
     */
    @Test
    void matchesHandComputedFrequenciesAndRatesOnATinyAlignment() {
        SequenceAlignment aln = SequenceAlignment.of(List.of(
                new NucleotideSequence("seq1", "AAAA"),
                new NucleotideSequence("seq2", "ACGT")
        ));

        GtrParameterEstimator.Estimate est = GtrParameterEstimator.estimate(aln);

        // frequencies include a +0.5 pseudocount per base (denominator +4*0.5=+2) so a base that's
        // genuinely never observed still gets a small positive frequency, required by GtrModel's eigendecomposition
        double freqA = (5 + 0.5) / 10.0, freqC = (1 + 0.5) / 10.0, freqG = (1 + 0.5) / 10.0, freqT = (1 + 0.5) / 10.0;
        assertThat(est.baseFrequencies()[Nucleotide.A]).isCloseTo(freqA, within(1e-9));
        assertThat(est.baseFrequencies()[Nucleotide.C]).isCloseTo(freqC, within(1e-9));
        assertThat(est.baseFrequencies()[Nucleotide.G]).isCloseTo(freqG, within(1e-9));
        assertThat(est.baseFrequencies()[Nucleotide.T]).isCloseTo(freqT, within(1e-9));

        // observed unordered pair counts: AC=1, AG=1, AT=1, CG=0, CT=0, GT=0
        double expectedAC = (1 + 0.5) / (freqA * freqC);
        double expectedAG = (1 + 0.5) / (freqA * freqG);
        double expectedAT = (1 + 0.5) / (freqA * freqT);
        double expectedCG = (0 + 0.5) / (freqC * freqG);
        double expectedCT = (0 + 0.5) / (freqC * freqT);
        double expectedGT = (0 + 0.5) / (freqG * freqT);

        // exchangeabilityRates order: AC, AG, AT, CG, CT, GT
        assertThat(est.exchangeabilityRates()[0]).isCloseTo(expectedAC, within(1e-6));
        assertThat(est.exchangeabilityRates()[1]).isCloseTo(expectedAG, within(1e-6));
        assertThat(est.exchangeabilityRates()[2]).isCloseTo(expectedAT, within(1e-6));
        assertThat(est.exchangeabilityRates()[3]).isCloseTo(expectedCG, within(1e-6));
        assertThat(est.exchangeabilityRates()[4]).isCloseTo(expectedCT, within(1e-6));
        assertThat(est.exchangeabilityRates()[5]).isCloseTo(expectedGT, within(1e-6));
    }

    @Test
    void producesAValidGtrModel() {
        // sanity/integration check: the estimated parameters must be directly usable by GtrModel (positive rates, freqs sum to 1)
        SequenceAlignment aln = SequenceAlignment.of(List.of(
                new NucleotideSequence("seq1", "ACGTACGTACGT"),
                new NucleotideSequence("seq2", "ACGTACGAACGT"),
                new NucleotideSequence("seq3", "ACGTTCGTACGT")
        ));
        GtrParameterEstimator.Estimate est = GtrParameterEstimator.estimate(aln);
        GtrModel model = GtrModel.of(est.baseFrequencies(), est.exchangeabilityRates());
        double[][] p = model.transitionProbabilities(0.1);
        for (double[] row : p) {
            double sum = 0;
            for (double v : row) sum += v;
            assertThat(sum).isCloseTo(1.0, within(1e-6));
        }
    }
}
