package org.gvi.algorithms.phylo.model;

import org.gvi.algorithms.phylo.NeighborJoining;
import org.gvi.algorithms.phylo.PhyloEdge;
import org.gvi.algorithms.phylo.PhyloTree;
import org.junit.jupiter.api.Test;

import java.util.List;
import java.util.Map;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.within;

class TreeLikelihoodCalculatorTest {

    private final TreeLikelihoodCalculator calc = new TreeLikelihoodCalculator();

    /**
     * Hand-verifiable star tree: A--(t1)--X, X--(t2)--B, X--(t3)--C, rooted
     * at A (the reference), single site: A='A', B='A', C='G'. Under JC69,
     * the exact (root-placement-invariant) likelihood of this tree is:
     * <pre>
     *   L = sum over ancestral states s at X of:
     *       0.25 * P(t1)[s,A] * P(t2)[s,A] * P(t3)[s,G]
     * </pre>
     * computed independently here from JC69's closed form and compared
     * against the pruning implementation -- this directly tests that the
     * root-taxon-is-observed-data restructuring described in
     * TreeLikelihoodCalculator's javadoc is handled correctly (a naive
     * pruning that treated the root as a free ancestral state would give a
     * different, wrong answer).
     */
    @Test
    void matchesHandComputedLikelihoodOnAThreeTaxonStarTree() {
        double t1 = 0.1, t2 = 0.2, t3 = 0.3;
        // Build the tree directly (bypassing NJ's distance-based reconstruction) so we control exact branch lengths:
        // taxa A=0, B=1, C=2; internal node X=3. Edges: A-X (t1), B-X (t2), C-X (t3).
        List<PhyloEdge> edges = List.of(
                new PhyloEdge(0, 3, t1),
                new PhyloEdge(1, 3, t2),
                new PhyloEdge(2, 3, t3)
        );
        Map<Integer, String> labels = Map.of(0, "A", 1, "B", 2, "C");
        NeighborJoining.Result njResult = new NeighborJoining.Result(edges, labels, List.of());
        PhyloTree tree = PhyloTree.rootAt(njResult, "A");

        Map<String, String> sequences = Map.of("A", "A", "B", "A", "C", "G");
        GtrModel jc = GtrModel.jukesCantor();

        double logL = calc.logLikelihood(tree, sequences, 1, jc, new double[]{1.0});

        double expectedL = handComputedStarTreeLikelihood(t1, t2, t3);
        assertThat(Math.exp(logL)).isCloseTo(expectedL, within(1e-9));
    }

    private double handComputedStarTreeLikelihood(double t1, double t2, double t3) {
        // states: 0=A,1=C,2=G,3=T. Observed: leaf-to-X edges see states {A(0), A(0), G(2)}.
        double total = 0;
        for (int s = 0; s < 4; s++) {
            double p1 = jcProb(t1, s, 0); // to A
            double p2 = jcProb(t2, s, 0); // to A
            double p3 = jcProb(t3, s, 2); // to G
            total += 0.25 * p1 * p2 * p3;
        }
        return total;
    }

    private double jcProb(double t, int from, int to) {
        double same = 0.25 + 0.75 * Math.exp(-4.0 * t / 3.0);
        double diff = 0.25 - 0.25 * Math.exp(-4.0 * t / 3.0);
        return from == to ? same : diff;
    }

    /**
     * A genuinely nested (non-star) 4-taxon tree: A--(1)--X--(3)--Y, X--(2)--B, Y--(1)--C, Y--(4)--D,
     * rooted at A. After rooting, X (pruningRoot) has children [B, Y], and Y has children [C, D] --
     * this exercises the recursive case (subtreeLikelihood calling itself), not just the shallow
     * star-tree case above. Verified against an independent nested-sum computation of the same
     * closed-form JC69 formula (sum over ancestral states at both X and Y).
     */
    @Test
    void matchesHandComputedLikelihoodOnANestedFourTaxonTree() {
        double tAX = 1, tXB = 2, tXY = 3, tYC = 1, tYD = 4;
        List<PhyloEdge> edges = List.of(
                new PhyloEdge(0, 4, tAX),  // A(0) -- X(4)
                new PhyloEdge(1, 4, tXB),  // B(1) -- X(4)
                new PhyloEdge(4, 5, tXY),  // X(4) -- Y(5)
                new PhyloEdge(2, 5, tYC),  // C(2) -- Y(5)
                new PhyloEdge(3, 5, tYD)   // D(3) -- Y(5)
        );
        Map<Integer, String> labels = Map.of(0, "A", 1, "B", 2, "C", 3, "D");
        PhyloTree tree = PhyloTree.rootAt(new NeighborJoining.Result(edges, labels, List.of()), "A");

        Map<String, String> sequences = Map.of("A", "A", "B", "A", "C", "G", "D", "T");
        GtrModel jc = GtrModel.jukesCantor();

        double logL = calc.logLikelihood(tree, sequences, 1, jc, new double[]{1.0});

        // independent nested-sum computation: ancestral states at X (sX) and Y (sY)
        double expected = 0;
        int obsA = 0, obsB = 0, obsC = 2, obsD = 3;
        for (int sX = 0; sX < 4; sX++) {
            double innerYSum = 0;
            for (int sY = 0; sY < 4; sY++) {
                innerYSum += jcProb(tXY, sX, sY) * jcProb(tYC, sY, obsC) * jcProb(tYD, sY, obsD);
            }
            expected += 0.25 * jcProb(tAX, sX, obsA) * jcProb(tXB, sX, obsB) * innerYSum;
        }

        assertThat(Math.exp(logL)).isCloseTo(expected, within(1e-9));
    }
}
