package org.gvi.algorithms.phylo;

import org.junit.jupiter.api.Test;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.within;

class NeighborJoiningTest {

    /**
     * Neighbor-Joining is provably exact (topology AND branch lengths) when
     * the input distance matrix is perfectly additive -- i.e. genuinely
     * derived from some tree (Studier & Keppler 1988). So instead of trying
     * to recall a textbook example from memory (risk of transcription
     * error), this test constructs the matrix itself from a known tree:
     * <p>
     * A --1-- X --3-- Y --1-- C
     *         |               |
     *         2               4
     *         |               |
     *         B               D
     * <p>
     * i.e. ((A:1,B:2):3,(C:1,D:4)) in Newick, with pairwise path-length
     * distances d(A,B)=3, d(A,C)=5, d(A,D)=8, d(B,C)=6, d(B,D)=9, d(C,D)=5
     * computed by hand from that tree. NJ must reconstruct it exactly.
     */
    @Test
    void reconstructsExactBranchLengthsFromAnAdditiveDistanceMatrix() {
        String[] labels = {"A", "B", "C", "D"};
        double[][] distances = {
                {0, 3, 5, 8},
                {3, 0, 6, 9},
                {5, 6, 0, 5},
                {8, 9, 5, 0}
        };

        NeighborJoining.Result result = NeighborJoining.build(labels, distances);
        PhyloTree tree = PhyloTree.rootAt(result, "A");

        // rooted at A (distance 0 to itself), root-to-tip to every other taxon
        // must equal the original matrix row for A -- this holds for ANY exactly
        // additive matrix and is a strong end-to-end correctness check.
        assertThat(tree.rootToTip("B")).isCloseTo(3.0, within(1e-9));
        assertThat(tree.rootToTip("C")).isCloseTo(5.0, within(1e-9));
        assertThat(tree.rootToTip("D")).isCloseTo(8.0, within(1e-9));

        // exactly 2 true branching points (internal nodes) for 4 taxa (n-2)
        assertThat(tree.branchPoints()).hasSize(2);
        assertThat(tree.taxonCount()).isEqualTo(4);
    }

    @Test
    void rootToTipIsInternallyConsistentRegardlessOfWhichTaxonIsChosenAsRoot() {
        // Same tree, but root at C instead of A this time -- root-to-tip to
        // every other taxon must equal the original matrix row for C.
        String[] labels = {"A", "B", "C", "D"};
        double[][] distances = {
                {0, 3, 5, 8},
                {3, 0, 6, 9},
                {5, 6, 0, 5},
                {8, 9, 5, 0}
        };
        NeighborJoining.Result result = NeighborJoining.build(labels, distances);
        PhyloTree tree = PhyloTree.rootAt(result, "C");

        assertThat(tree.rootToTip("A")).isCloseTo(5.0, within(1e-9));
        assertThat(tree.rootToTip("B")).isCloseTo(6.0, within(1e-9));
        assertThat(tree.rootToTip("D")).isCloseTo(5.0, within(1e-9));
    }

    @Test
    void worksWithFiveTaxaAdditiveMatrix() {
        // Star-like additive tree: center Z connects to internal X (parent of A,B) and directly to C, D, E.
        // A-X=1, B-X=2, X-Z=1, C-Z=2, D-Z=3, E-Z=1
        String[] labels = {"A", "B", "C", "D", "E"};
        double[][] distances = new double[5][5];
        // path distances computed by hand from the tree:
        // A-X=1, B-X=2, X-Z=1, C-Z=2, D-Z=3, E-Z=1
        double AX = 1, BX = 2, XZ = 1, CZ = 2, DZ = 3, EZ = 1;
        distances[0][1] = distances[1][0] = AX + BX;               // A-B
        distances[0][2] = distances[2][0] = AX + XZ + CZ;          // A-C
        distances[0][3] = distances[3][0] = AX + XZ + DZ;          // A-D
        distances[0][4] = distances[4][0] = AX + XZ + EZ;          // A-E
        distances[1][2] = distances[2][1] = BX + XZ + CZ;          // B-C
        distances[1][3] = distances[3][1] = BX + XZ + DZ;          // B-D
        distances[1][4] = distances[4][1] = BX + XZ + EZ;          // B-E
        distances[2][3] = distances[3][2] = CZ + DZ;               // C-D
        distances[2][4] = distances[4][2] = CZ + EZ;               // C-E
        distances[3][4] = distances[4][3] = DZ + EZ;               // D-E

        NeighborJoining.Result result = NeighborJoining.build(labels, distances);
        PhyloTree tree = PhyloTree.rootAt(result, "A");

        assertThat(tree.rootToTip("B")).isCloseTo(distances[0][1], within(1e-9));
        assertThat(tree.rootToTip("C")).isCloseTo(distances[0][2], within(1e-9));
        assertThat(tree.rootToTip("D")).isCloseTo(distances[0][3], within(1e-9));
        assertThat(tree.rootToTip("E")).isCloseTo(distances[0][4], within(1e-9));
        assertThat(tree.branchPoints()).hasSize(3); // n-2 = 3 for 5 taxa
    }
}
