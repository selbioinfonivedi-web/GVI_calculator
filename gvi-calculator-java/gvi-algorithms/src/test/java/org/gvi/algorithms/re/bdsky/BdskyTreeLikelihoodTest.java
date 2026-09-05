package org.gvi.algorithms.re.bdsky;

import org.gvi.algorithms.gd.GdMethod;
import org.gvi.algorithms.phylo.PhyloNode;
import org.gvi.algorithms.phylo.PhyloTree;
import org.gvi.algorithms.phylo.PhyloTreeFactory;
import org.gvi.core.model.NucleotideSequence;
import org.gvi.core.model.SequenceAlignment;
import org.junit.jupiter.api.Test;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.within;

/**
 * Verifies the tree-recursion correctly COMPOSES the already-validated
 * propagator outputs ({@link BirthDeathSamplingPropagatorTest}) -- this
 * layer's job is the composition logic, so the expected value here is
 * built directly from independent calls to {@link BirthDeathSamplingPropagator},
 * not from {@link BdskyTreeLikelihood} itself.
 * <p>
 * Reflects the two bugs documented in {@link BdskyTreeLikelihood}'s class
 * javadoc, found via Monte Carlo cross-checks against a from-scratch
 * simulator (not by this test, which only checks the composition logic is
 * internally self-consistent): (1) p/logG must be evaluated at each node's
 * remaining time to the GLOBAL present, not per-edge-local elapsed time --
 * so dates are chosen here so P and Q are at DIFFERENT distances from the
 * present, deliberately exercising that; (2) the reference tip's own
 * contribution and its edge to M are dropped entirely (sampling-with-removal
 * is inconsistent with treating a sampled tip as also having descendants).
 * <p>
 * Uses a real 3-taxon tree built via the public alignment/NJ API (so the
 * topology -- reference rooted, with one child M that itself has two
 * children P and Q -- comes from real reconstruction, not a hand-built
 * mock), but assigns dates directly (independent of the substitution
 * branch lengths NJ produced) since the BDSKY likelihood only consumes
 * elapsed time between assigned dates, not substitution distance.
 */
class BdskyTreeLikelihoodTest {

    @Test
    void composesPropagatorOutputsCorrectlyOnGlobalRemainingTimeOnAMinimalThreeTaxonTree() {
        NucleotideSequence ref = new NucleotideSequence("ref", "A".repeat(1000));
        NucleotideSequence p = new NucleotideSequence("P", mutateAt(1000, 0));
        NucleotideSequence q = new NucleotideSequence("Q", mutateAt(1000, 1));
        SequenceAlignment alignment = SequenceAlignment.of(List.of(ref, p, q), "ref");
        PhyloTreeFactory.BuiltTree built = PhyloTreeFactory.build(alignment, GdMethod.JUKES_CANTOR);
        PhyloTree tree = built.tree();

        PhyloNode root = tree.root();
        assertThat(root.children()).hasSize(1);
        PhyloNode m = root.children().get(0);
        assertThat(m.children()).hasSize(2);
        PhyloNode p1 = m.children().get(0);
        PhyloNode p2 = m.children().get(1);

        Map<PhyloNode, Double> dates = new HashMap<>();
        dates.put(root, 0.0);
        dates.put(m, 2.0);
        dates.put(p1, 5.0); // R(p1) = presentTime(6) - 5 = 1
        dates.put(p2, 6.0); // R(p2) = presentTime(6) - 6 = 0 -- the most recent tip defines "the present"

        double lambda = 1.3, mu = 0.4, psi = 0.7;
        BdskyTreeLikelihood likelihood = new BdskyTreeLikelihood(lambda, mu, psi);
        double actual = likelihood.logLikelihood(root, dates);

        BirthDeathSamplingPropagator propagator = new BirthDeathSamplingPropagator(lambda, mu, psi);
        double logGAtM = propagator.integrate(4.0).logG();  // R(M) = 6 - 2 = 4
        double logGAtP1 = propagator.integrate(1.0).logG(); // R(p1) = 1
        // R(p2) = 0, and logG(0) = 0 always (see BirthDeathSamplingPropagatorTest), so that term drops out.
        // M contributes log(lambda) (a genuine branching event); reference contributes NOTHING (bug #2);
        // edge M->p1 contributes logG(R(M)) - logG(R(p1)); edge M->p2 contributes logG(R(M)) - logG(R(p2)=0);
        // p1 and p2 each contribute log(psi) as ordinary sampled tips.
        double expected = Math.log(lambda)
                + (logGAtM - logGAtP1) + Math.log(psi)
                + (logGAtM - 0.0) + Math.log(psi);

        assertThat(actual).isCloseTo(expected, within(1e-9));
    }

    private static String mutateAt(int length, int position) {
        char[] seq = new char[length];
        java.util.Arrays.fill(seq, 'A');
        seq[position] = 'G';
        return new String(seq);
    }
}
