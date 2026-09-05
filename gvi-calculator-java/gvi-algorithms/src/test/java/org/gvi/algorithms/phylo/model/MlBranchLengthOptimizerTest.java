package org.gvi.algorithms.phylo.model;

import org.gvi.algorithms.gd.GdMethod;
import org.gvi.algorithms.phylo.PhyloNode;
import org.gvi.algorithms.phylo.PhyloTree;
import org.gvi.algorithms.phylo.PhyloTreeFactory;
import org.gvi.core.model.NucleotideSequence;
import org.gvi.core.model.SequenceAlignment;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.within;

class MlBranchLengthOptimizerTest {

    private final MlBranchLengthOptimizer optimizer = new MlBranchLengthOptimizer();

    /**
     * A, B, C each carry a private, disjoint block of A->G substitutions
     * relative to an implicit all-A ancestor (star topology X), over a very
     * long alignment so p-distance ~= true branch length with negligible
     * multiple-hit noise -- same construction technique as
     * EvolutionaryRateCalculatorTest's tree-aware mu test. True generating
     * branch lengths: t(A-X)=0.001, t(B-X)=0.002, t(C-X)=0.003.
     * <p>
     * Branch lengths are then deliberately overwritten to a clearly wrong
     * constant (0.5) BEFORE optimizing, so this tests that the optimizer
     * itself -- not just NJ's already-good distance-based starting guess --
     * recovers values close to the true generating branch lengths.
     */
    @Test
    void recoversApproximatelyTrueBranchLengthsFromADeliberatelyWrongStartingPoint() {
        // Kept small (20,000bp, not 200,000) so the O(sites x branches x cycles x Brent-evals) optimization
        // loop finishes in a reasonable test runtime -- still long enough that p ~= true branch length closely.
        int length = 20_000;
        char[] base = new char[length];
        java.util.Arrays.fill(base, 'A');

        NucleotideSequence a = mutate("A", base, 100, 120);  // 20 positions -> t = 20/20000 = 0.001
        NucleotideSequence b = mutate("B", base, 0, 40);     // 40 positions -> t = 0.002
        NucleotideSequence c = mutate("C", base, 40, 100);   // 60 positions -> t = 0.003

        SequenceAlignment aln = SequenceAlignment.of(List.of(a, b, c), "A");
        PhyloTreeFactory.BuiltTree built = PhyloTreeFactory.build(aln, GdMethod.JUKES_CANTOR);
        PhyloTree tree = built.tree();

        // deliberately wreck the (already-good) NJ branch lengths before optimizing
        overwriteAllBranchLengths(tree, 0.5);

        Map<String, String> sequences = new HashMap<>();
        for (NucleotideSequence seq : List.of(a, b, c)) sequences.put(seq.getId(), seq.getSequence());

        GtrModel jc = GtrModel.jukesCantor();
        double[] singleCategory = {1.0};

        MlBranchLengthOptimizer.Result result = optimizer.optimize(tree, sequences, length, jc, singleCategory);

        // optimization must never make likelihood worse than where it started
        assertThat(result.finalLogLikelihood()).isGreaterThanOrEqualTo(result.initialLogLikelihood());
        // and starting from 0.5 everywhere, it must have improved substantially (0.5 is a terrible fit for ~0.001-0.003 true distances)
        assertThat(result.finalLogLikelihood()).isGreaterThan(result.initialLogLikelihood() + 100);

        Map<String, Double> branchLengths = new HashMap<>();
        collectTaxonBranchLengths(tree.root(), branchLengths);

        // A is the root taxon: its edge-to-the-rest-of-the-tree length is stored on its own
        // single child (the internal node), not on A itself -- see TreeLikelihoodCalculator's javadoc.
        double aBranchLength = tree.root().children().get(0).branchLength();

        assertThat(aBranchLength).isCloseTo(0.001, within(0.0005));
        assertThat(branchLengths.get("B")).isCloseTo(0.002, within(0.0005));
        assertThat(branchLengths.get("C")).isCloseTo(0.003, within(0.0005));
    }

    private void overwriteAllBranchLengths(PhyloTree tree, double value) {
        overwrite(tree.root(), value);
    }

    private void overwrite(PhyloNode node, double value) {
        for (PhyloNode child : node.children()) {
            child.setBranchLength(value);
            overwrite(child, value);
        }
    }

    private void collectTaxonBranchLengths(PhyloNode node, Map<String, Double> out) {
        if (node.isTaxon() && !node.isRoot()) {
            out.put(node.label(), node.branchLength());
        }
        for (PhyloNode child : node.children()) collectTaxonBranchLengths(child, out);
    }

    private static NucleotideSequence mutate(String id, char[] base, int fromInclusive, int toExclusive) {
        char[] copy = base.clone();
        for (int i = fromInclusive; i < toExclusive; i++) copy[i] = 'G';
        return new NucleotideSequence(id, new String(copy));
    }
}
