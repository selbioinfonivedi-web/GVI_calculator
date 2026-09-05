package org.gvi.algorithms.phylo.model;

import org.gvi.algorithms.phylo.PhyloNode;
import org.gvi.algorithms.phylo.PhyloTree;
import org.gvi.core.exception.GviComputationException;

import java.util.Arrays;
import java.util.Map;

/**
 * Felsenstein (1981) pruning algorithm: exact log-likelihood of a sequence
 * alignment given a fixed tree (topology + branch lengths) under a
 * GTR(+Gamma) substitution model. This is the objective function ML tree
 * search / branch-length optimization maximizes -- not an approximation of
 * anything, the standard exact algorithm real tools (RAxML, IQ-TREE,
 * PAML's baseml) use for this computation.
 * <p>
 * One important wrinkle specific to how {@link PhyloTree} is built here:
 * it's rooted AT an observed taxon (the reference sequence), not at a
 * hypothetical unobserved ancestor the way textbook pruning diagrams
 * usually show. That taxon's actual sequence data must be folded in as a
 * leaf constraint via its own edge -- treating it as a free summed-over
 * ancestral state (which a naive "just recurse from tree.root()" pruning
 * would do, since the root has exactly one child in this rooting scheme)
 * would silently throw away the reference's own data and give a wrong,
 * likelihood-inflating answer. Because GTR is time-reversible, the
 * likelihood is provably identical regardless of where the computation is
 * rooted, so this restructuring (compute from the root's child, folding
 * the root itself in as an extra pendant leaf on that same edge) is exact,
 * not an approximation.
 */
public final class TreeLikelihoodCalculator {

    public double logLikelihood(PhyloTree tree, Map<String, String> sequencesByLabel, int alignmentLength,
                                 GtrModel model, double[] gammaCategoryRates) {
        if (alignmentLength <= 0) {
            throw new GviComputationException("Cannot compute tree likelihood: alignment length must be positive");
        }
        double totalLogLikelihood = 0.0;
        for (int site = 0; site < alignmentLength; site++) {
            double siteLikelihood = 0.0;
            for (double rate : gammaCategoryRates) {
                siteLikelihood += siteLikelihoodForCategory(tree, sequencesByLabel, site, model, rate);
            }
            siteLikelihood /= gammaCategoryRates.length;
            // A site pattern that is essentially impossible under the model/tree (e.g. saturated data)
            // can numerically underflow to 0 -- clamp rather than let log() produce -Infinity and poison the sum.
            totalLogLikelihood += Math.log(Math.max(siteLikelihood, Double.MIN_NORMAL));
        }
        return totalLogLikelihood;
    }

    private double siteLikelihoodForCategory(PhyloTree tree, Map<String, String> sequences, int site, GtrModel model, double rate) {
        PhyloNode root = tree.root();
        if (root.children().size() != 1) {
            throw new GviComputationException("Expected a PhyloTree rooted at a taxon (exactly one child of root), found "
                    + root.children().size() + " -- tree was not built via PhyloTree.rootAt");
        }
        PhyloNode pruningRoot = root.children().get(0);
        double edgeToObservedRoot = pruningRoot.branchLength();

        double[] combined = filled(1.0);
        for (PhyloNode child : pruningRoot.children()) {
            multiplyIn(combined, edgeContribution(child, sequences, site, model, rate));
        }
        // fold in the reference taxon itself as an extra pendant leaf on the same edge (see class javadoc)
        double[] rootLeaf = leafPartialLikelihood(root, sequences, site);
        double[][] pRootEdge = model.transitionProbabilities(edgeToObservedRoot * rate);
        multiplyIn(combined, transform(pRootEdge, rootLeaf));

        double[] freq = model.baseFrequencies();
        double total = 0;
        for (int i = 0; i < 4; i++) total += freq[i] * combined[i];
        return total;
    }

    /** Likelihood contribution of the subtree rooted at {@code node}, propagated back through its own branch length. */
    private double[] edgeContribution(PhyloNode node, Map<String, String> sequences, int site, GtrModel model, double rate) {
        double[] nodeLikelihood = subtreeLikelihood(node, sequences, site, model, rate);
        double[][] p = model.transitionProbabilities(node.branchLength() * rate);
        return transform(p, nodeLikelihood);
    }

    private double[] subtreeLikelihood(PhyloNode node, Map<String, String> sequences, int site, GtrModel model, double rate) {
        if (node.children().isEmpty()) {
            return leafPartialLikelihood(node, sequences, site);
        }
        double[] combined = filled(1.0);
        for (PhyloNode child : node.children()) {
            multiplyIn(combined, edgeContribution(child, sequences, site, model, rate));
        }
        return combined;
    }

    private double[] leafPartialLikelihood(PhyloNode node, Map<String, String> sequences, int site) {
        String seq = sequences.get(node.label());
        if (seq == null) {
            throw new GviComputationException("No sequence data found for taxon '" + node.label() + "'");
        }
        return Nucleotide.partialLikelihood(seq.charAt(site));
    }

    private double[] transform(double[][] p, double[] childLikelihood) {
        double[] out = new double[4];
        for (int i = 0; i < 4; i++) {
            double sum = 0;
            for (int j = 0; j < 4; j++) sum += p[i][j] * childLikelihood[j];
            out[i] = sum;
        }
        return out;
    }

    private void multiplyIn(double[] a, double[] b) {
        for (int i = 0; i < a.length; i++) a[i] *= b[i];
    }

    private double[] filled(double value) {
        double[] a = new double[4];
        Arrays.fill(a, value);
        return a;
    }
}
