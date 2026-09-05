package org.gvi.algorithms.phylo.model;

import org.gvi.algorithms.gd.GdMethod;
import org.gvi.algorithms.phylo.PhyloTreeFactory;
import org.gvi.core.model.NucleotideSequence;
import org.gvi.core.model.SequenceAlignment;
import org.junit.jupiter.api.Test;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static org.assertj.core.api.Assertions.assertThat;

class MlPhylogeneticOptimizerTest {

    private final MlPhylogeneticOptimizer optimizer = new MlPhylogeneticOptimizer();

    /**
     * All substitutions in this dataset are A<->G, which is a transition
     * (both purines) -- zero transversions anywhere. Under HKY85, the only
     * way to explain "many transitions, zero transversions" is a kappa
     * (transition/transversion rate ratio) well above 1. This directly
     * tests that kappa is genuinely being fit to the data's actual
     * ts/tv signal, not just sitting at its initial guess.
     */
    @Test
    void recoversAStrongTransitionBiasAsAHighKappa() {
        int length = 400;
        char[] base = new char[length];
        java.util.Arrays.fill(base, 'A');

        NucleotideSequence ref = new NucleotideSequence("ref", new String(base));
        NucleotideSequence q1 = mutate(base, 0, 15);   // 15 A->G transitions
        NucleotideSequence q2 = mutate(base, 15, 35);  // 20 A->G transitions

        SequenceAlignment aln = SequenceAlignment.of(List.of(ref, q1, q2), "ref");
        PhyloTreeFactory.BuiltTree built = PhyloTreeFactory.build(aln, GdMethod.JUKES_CANTOR);

        Map<String, String> sequences = new HashMap<>();
        for (NucleotideSequence s : List.of(ref, q1, q2)) sequences.put(s.getId(), s.getSequence());

        GtrParameterEstimator.Estimate estimate = GtrParameterEstimator.estimate(aln);
        RateParameterization hky85 = NamedSubstitutionModels.hky85();

        MlPhylogeneticOptimizer.Result result = optimizer.optimize(
                built.tree(), sequences, length, hky85, estimate.baseFrequencies(), 0, 1.0);

        assertThat(result.finalLogLikelihood()).isGreaterThanOrEqualTo(result.initialLogLikelihood());
        assertThat(result.optimizedTheta()).hasSize(1); // kappa
        assertThat(result.optimizedTheta()[0]).as("ML-fit kappa should reflect the strong transition bias").isGreaterThan(3.0);
        assertThat(result.optimizedAlpha()).isNull(); // gamma disabled (gammaCategoryCount=0)
    }

    @Test
    void jc69HasNoFreeRateParametersAndJustOptimizesBranchLengths() {
        int length = 300;
        char[] base = new char[length];
        java.util.Arrays.fill(base, 'A');
        NucleotideSequence ref = new NucleotideSequence("ref", new String(base));
        NucleotideSequence q1 = mutate(base, 0, 10);
        NucleotideSequence q2 = mutate(base, 10, 25);

        SequenceAlignment aln = SequenceAlignment.of(List.of(ref, q1, q2), "ref");
        PhyloTreeFactory.BuiltTree built = PhyloTreeFactory.build(aln, GdMethod.JUKES_CANTOR);
        Map<String, String> sequences = new HashMap<>();
        for (NucleotideSequence s : List.of(ref, q1, q2)) sequences.put(s.getId(), s.getSequence());

        MlPhylogeneticOptimizer.Result result = optimizer.optimize(
                built.tree(), sequences, length, NamedSubstitutionModels.jc69(), new double[]{0.25, 0.25, 0.25, 0.25}, 0, 1.0);

        assertThat(result.optimizedTheta()).isEmpty();
        assertThat(result.finalLogLikelihood()).isGreaterThanOrEqualTo(result.initialLogLikelihood());
    }

    @Test
    void gammaEnabledProducesAFiniteAlphaWithinBounds() {
        int length = 300;
        char[] base = new char[length];
        java.util.Arrays.fill(base, 'A');
        NucleotideSequence ref = new NucleotideSequence("ref", new String(base));
        NucleotideSequence q1 = mutate(base, 0, 10);
        NucleotideSequence q2 = mutate(base, 10, 25);

        SequenceAlignment aln = SequenceAlignment.of(List.of(ref, q1, q2), "ref");
        PhyloTreeFactory.BuiltTree built = PhyloTreeFactory.build(aln, GdMethod.JUKES_CANTOR);
        Map<String, String> sequences = new HashMap<>();
        for (NucleotideSequence s : List.of(ref, q1, q2)) sequences.put(s.getId(), s.getSequence());

        MlPhylogeneticOptimizer.Result result = optimizer.optimize(
                built.tree(), sequences, length, NamedSubstitutionModels.jc69(), new double[]{0.25, 0.25, 0.25, 0.25}, 4, 1.0);

        assertThat(result.optimizedAlpha()).isNotNull();
        assertThat(result.optimizedAlpha()).isBetween(0.02, 50.0);
        assertThat(result.finalGammaRates()).hasSize(4);
    }

    private static NucleotideSequence mutate(char[] base, int fromInclusive, int toExclusive) {
        char[] copy = base.clone();
        for (int i = fromInclusive; i < toExclusive; i++) copy[i] = 'G'; // A->G, a transition
        return new NucleotideSequence("q" + fromInclusive, new String(copy));
    }
}
