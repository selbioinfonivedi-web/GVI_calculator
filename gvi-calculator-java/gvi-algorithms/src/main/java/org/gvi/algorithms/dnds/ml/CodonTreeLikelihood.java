package org.gvi.algorithms.dnds.ml;

import org.gvi.algorithms.phylo.PhyloNode;
import org.gvi.algorithms.phylo.PhyloTree;
import org.gvi.core.exception.GviComputationException;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Felsenstein (1981) pruning at the codon level (61 states instead of the
 * usual 4 nucleotide states) under {@link CodonModel} -- the exact same
 * algorithm {@link org.gvi.algorithms.phylo.model.TreeLikelihoodCalculator}
 * uses for nucleotides, generalized in state count, and reusing the same
 * root-taxon-folding handling for this codebase's reference-rooted
 * {@link PhyloTree} convention (see that class's javadoc for why the root
 * taxon's own data has to be folded in explicitly rather than left as a
 * free summed ancestral state). No Gamma rate heterogeneity here: GY94's
 * M0 "one-ratio" model (what {@code codeml} calls it) applies one omega to
 * every site by construction, matching {@code codeml}'s own M0 baseline.
 * <p>
 * Each branch's P(t)=exp(Qt) matrix is an O(61^3) computation (the triple
 * sum in {@link CodonModel#transitionProbabilities}) but is the SAME for
 * every site given a fixed branch length -- so callers doing repeated
 * single-branch optimization (see {@link CodonMlFitter}) build the cache
 * once via {@link #buildTransitionMatrices} and only replace the ONE
 * changed branch's entry per candidate evaluation via
 * {@link #logLikelihoodWithMatrices}, rather than rebuilding every branch's
 * matrix on every evaluation -- a real, measured performance bug caught
 * during development: rebuilding per-site (not even per-eval) made a
 * 150-codon/5-taxon ground-truth test run 20+ CPU-minutes without finishing.
 */
final class CodonTreeLikelihood {

    double logLikelihood(PhyloTree tree, Map<String, List<String>> codonsByLabel, int siteCount, CodonModel model) {
        Map<PhyloNode, double[][]> matrices = buildTransitionMatrices(tree, model);
        return logLikelihoodWithMatrices(tree, codonsByLabel, siteCount, matrices, model.frequencies());
    }

    /** Builds every branch's P(t) matrix once for the given model -- call this whenever the MODEL (kappa/omega) changes. */
    Map<PhyloNode, double[][]> buildTransitionMatrices(PhyloTree tree, CodonModel model) {
        PhyloNode pruningRoot = pruningRoot(tree);
        Map<PhyloNode, double[][]> matrices = new HashMap<>();
        matrices.put(pruningRoot, model.transitionProbabilities(pruningRoot.branchLength()));
        cacheTransitionMatrices(pruningRoot, model, matrices);
        return matrices;
    }

    private void cacheTransitionMatrices(PhyloNode node, CodonModel model, Map<PhyloNode, double[][]> out) {
        for (PhyloNode child : node.children()) {
            out.put(child, model.transitionProbabilities(child.branchLength()));
            cacheTransitionMatrices(child, model, out);
        }
    }

    /** Log-likelihood using a caller-supplied, already-built matrix cache -- see class javadoc for why this matters for performance. */
    double logLikelihoodWithMatrices(PhyloTree tree, Map<String, List<String>> codonsByLabel, int siteCount,
                                      Map<PhyloNode, double[][]> transitionMatrices, double[] freq) {
        if (siteCount <= 0) {
            throw new GviComputationException("Cannot compute codon tree likelihood: site count must be positive");
        }
        PhyloNode root = tree.root();
        PhyloNode pruningRoot = pruningRoot(tree);

        double total = 0.0;
        for (int site = 0; site < siteCount; site++) {
            double siteLikelihood = siteLikelihood(root, pruningRoot, codonsByLabel, site, transitionMatrices, freq);
            total += Math.log(Math.max(siteLikelihood, Double.MIN_NORMAL));
        }
        return total;
    }

    private PhyloNode pruningRoot(PhyloTree tree) {
        PhyloNode root = tree.root();
        if (root.children().size() != 1) {
            throw new GviComputationException("Expected a PhyloTree rooted at a taxon (exactly one child of root), found "
                    + root.children().size() + " -- tree was not built via PhyloTree.rootAt");
        }
        return root.children().get(0);
    }

    private double siteLikelihood(PhyloNode root, PhyloNode pruningRoot, Map<String, List<String>> codons, int site,
                                   Map<PhyloNode, double[][]> transitionMatrices, double[] freq) {
        double[] combined = filled(1.0);
        for (PhyloNode child : pruningRoot.children()) {
            multiplyIn(combined, edgeContribution(child, codons, site, transitionMatrices));
        }
        double[] rootLeaf = leafPartialLikelihood(root, codons, site);
        multiplyIn(combined, transform(transitionMatrices.get(pruningRoot), rootLeaf));

        double total = 0;
        for (int i = 0; i < CodonAlphabet.SIZE; i++) total += freq[i] * combined[i];
        return total;
    }

    private double[] edgeContribution(PhyloNode node, Map<String, List<String>> codons, int site, Map<PhyloNode, double[][]> transitionMatrices) {
        double[] nodeLikelihood = subtreeLikelihood(node, codons, site, transitionMatrices);
        return transform(transitionMatrices.get(node), nodeLikelihood);
    }

    private double[] subtreeLikelihood(PhyloNode node, Map<String, List<String>> codons, int site, Map<PhyloNode, double[][]> transitionMatrices) {
        if (node.children().isEmpty()) {
            return leafPartialLikelihood(node, codons, site);
        }
        double[] combined = filled(1.0);
        for (PhyloNode child : node.children()) {
            multiplyIn(combined, edgeContribution(child, codons, site, transitionMatrices));
        }
        return combined;
    }

    private double[] leafPartialLikelihood(PhyloNode node, Map<String, List<String>> codons, int site) {
        List<String> seq = codons.get(node.label());
        if (seq == null) {
            throw new GviComputationException("No codon data found for taxon '" + node.label() + "'");
        }
        String codon = seq.get(site);
        double[] out = new double[CodonAlphabet.SIZE];
        if (CodonAlphabet.isSenseCodon(codon)) {
            out[CodonAlphabet.index(codon)] = 1.0;
        } else {
            // stop codon (alignment/indel noise) or gap/ambiguity code: uninformative at this site for this taxon
            Arrays.fill(out, 1.0);
        }
        return out;
    }

    private double[] transform(double[][] p, double[] childLikelihood) {
        int n = CodonAlphabet.SIZE;
        double[] out = new double[n];
        for (int i = 0; i < n; i++) {
            double[] pRow = p[i];
            double sum = 0;
            for (int j = 0; j < n; j++) sum += pRow[j] * childLikelihood[j];
            out[i] = sum;
        }
        return out;
    }

    private void multiplyIn(double[] a, double[] b) {
        for (int i = 0; i < a.length; i++) a[i] *= b[i];
    }

    private double[] filled(double value) {
        double[] a = new double[CodonAlphabet.SIZE];
        Arrays.fill(a, value);
        return a;
    }
}
