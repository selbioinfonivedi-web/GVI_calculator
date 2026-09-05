package org.gvi.algorithms.phylo.model;

import org.gvi.algorithms.phylo.NeighborJoining;
import org.gvi.algorithms.phylo.PhyloEdge;
import org.gvi.algorithms.phylo.PhyloNode;
import org.gvi.algorithms.phylo.PhyloTree;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import static org.assertj.core.api.Assertions.assertThat;

class NniTopologySearchTest {

    private final NniTopologySearch search = new NniTopologySearch();
    private final TreeLikelihoodCalculator calculator = new TreeLikelihoodCalculator();

    /**
     * 5-taxon tree (A=reference/root, X=pruningRoot with children Y,Z; Y and Z
     * each a 2-taxon cherry), deliberately started in a low-likelihood
     * grouping: C and E are built near-identical (the two sequences with real
     * phylogenetic signal to cluster) but placed in *different* cherries
     * (Y={C,B}, Z={E,D}), while B and D are each distinctive outliers. Real
     * signal like this should make some neighboring topology strictly more
     * likely than the start, so the search should find and apply an
     * improving move (this deliberately does not hand-predict the exact
     * resulting bipartition -- with 2-taxon cherries, a single NNI move
     * always swaps a whole cherry, not individual members, so which specific
     * rearrangement wins is left to the likelihood surface; what's verified
     * is that real, structurally-valid progress happens).
     */
    @Test
    void findsAnImprovingMoveWhenStartingTopologyMisplacesSimilarSequences() {
        int length = 40;
        String a = "A".repeat(length);
        String cAndE = "A".repeat(length); // C and E: identical to the reference, to each other
        String b = "G".repeat(10) + "A".repeat(length - 10); // distinctive
        String d = "A".repeat(20) + "T".repeat(10) + "A".repeat(length - 30); // distinctive

        PhyloTree tree = buildFiveTaxonTree(0.05, 0.05, 0.05, 0.05, 0.05);
        Map<String, String> sequences = Map.of("A", a, "B", b, "C", cAndE, "D", d, "E", cAndE);
        GtrModel jc = GtrModel.jukesCantor();
        double[] gammaRates = {1.0};

        double initialLogL = calculator.logLikelihood(tree, sequences, length, jc, gammaRates);
        NniTopologySearch.Result result = search.search(tree, sequences, length, jc, gammaRates);

        assertThat(result.initialLogLikelihood()).isCloseTo(initialLogL, org.assertj.core.api.Assertions.within(1e-9));
        assertThat(result.movesAccepted()).isGreaterThan(0);
        assertThat(result.finalLogLikelihood()).isGreaterThan(result.initialLogLikelihood());

        double recomputed = calculator.logLikelihood(tree, sequences, length, jc, gammaRates);
        assertThat(recomputed).isCloseTo(result.finalLogLikelihood(), org.assertj.core.api.Assertions.within(1e-9));

        assertTreeStructurallyIntact(tree, Set.of("A", "B", "C", "D", "E"));
    }

    /**
     * When every sequence is identical, no topology is any more likely than another -- the search must be a
     * genuine no-op. Branch lengths are first optimized to convergence on the starting topology (via
     * {@link MlBranchLengthOptimizer}) so that every branch is already individually near its own optimum;
     * without that, a "no move" test would be confounded by NNI's per-candidate quick branch reoptimization
     * simply pulling under-optimized starting branches toward their optimum regardless of topology, which
     * is a real, separate effect and not evidence of a topology preference.
     */
    @Test
    void makesNoMoveWhenNoTopologyIsAnyBetterThanAnother() {
        int length = 20;
        String same = "A".repeat(length);
        PhyloTree tree = buildFiveTaxonTree(0.1, 0.1, 0.1, 0.1, 0.1);
        Map<String, String> sequences = Map.of("A", same, "B", same, "C", same, "D", same, "E", same);
        GtrModel jc = GtrModel.jukesCantor();
        double[] gammaRates = {1.0};

        new MlBranchLengthOptimizer().optimize(tree, sequences, length, jc, gammaRates);

        NniTopologySearch.Result result = search.search(tree, sequences, length, jc, gammaRates);

        assertThat(result.movesAccepted()).isZero();
        assertThat(result.finalLogLikelihood()).isCloseTo(result.initialLogLikelihood(), org.assertj.core.api.Assertions.within(1e-9));
    }

    /** A=0, B=1, C=2, D=3, E=4 (root); X=5 (pruningRoot), Y=6 (cherry B,C), Z=7 (cherry D,E). */
    private PhyloTree buildFiveTaxonTree(double tAX, double tXY, double tXZ, double tYB, double tYC) {
        List<PhyloEdge> edges = List.of(
                new PhyloEdge(0, 5, tAX),
                new PhyloEdge(5, 6, tXY),
                new PhyloEdge(5, 7, tXZ),
                new PhyloEdge(1, 6, tYB),
                new PhyloEdge(2, 6, tYC),
                new PhyloEdge(3, 7, 0.05),
                new PhyloEdge(4, 7, 0.05)
        );
        Map<Integer, String> labels = Map.of(0, "A", 1, "B", 2, "C", 3, "D", 4, "E");
        return PhyloTree.rootAt(new NeighborJoining.Result(edges, labels, List.of()), "A");
    }

    private void assertTreeStructurallyIntact(PhyloTree tree, Set<String> expectedTaxa) {
        Set<String> foundTaxa = new HashSet<>();
        List<PhyloNode> stack = new ArrayList<>();
        stack.add(tree.root());
        while (!stack.isEmpty()) {
            PhyloNode node = stack.remove(stack.size() - 1);
            if (node.isTaxon()) {
                assertThat(foundTaxa.add(node.label())).as("taxon %s should appear exactly once", node.label()).isTrue();
            } else if (!node.isRoot()) {
                assertThat(node.children()).as("every non-root internal node must remain strictly bifurcating").hasSize(2);
            }
            stack.addAll(node.children());
        }
        assertThat(foundTaxa).isEqualTo(expectedTaxa);
    }
}
