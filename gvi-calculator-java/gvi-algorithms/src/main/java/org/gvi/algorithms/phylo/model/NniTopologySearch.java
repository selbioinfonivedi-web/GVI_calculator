package org.gvi.algorithms.phylo.model;

import org.gvi.algorithms.phylo.PhyloNode;
import org.gvi.algorithms.phylo.PhyloTree;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * Nearest-Neighbor-Interchange (NNI) topology search: hill-climbs away from
 * Neighbor-Joining's distance-based starting topology toward a
 * higher-likelihood one, under a fixed GTR(+Gamma) model.
 * <p>
 * For every internal edge (u, v) -- u's other child w, v's two children
 * v1/v2 -- there are exactly two alternative topologies reachable by a
 * single NNI move: swap w with v1, or swap w with v2 (swapping w with u's
 * own parent-side subtree is not a separate move; it produces the same two
 * bipartitions, so scanning every qualifying edge this way already covers
 * the full NNI neighborhood without needing to special-case the path back
 * toward the root). {@link PhyloTree} is rooted at an observed taxon (see
 * {@link TreeLikelihoodCalculator}'s class javadoc); edges considered here
 * are always strictly between two internal (non-taxon) nodes, so this never
 * touches the reference taxon's own attachment point.
 * <p>
 * Each round scores every candidate move by reoptimizing only the two
 * branch lengths the swap directly touches (a cheap, local Brent pass --
 * not a full-tree cycle) before comparing likelihoods; this keeps candidate
 * scoring fast enough to scan the whole neighborhood, at the cost of a
 * slightly noisy signal for moves whose benefit depends on jointly
 * adjusting branches further away. The single best-improving move per round
 * is applied for real (reoptimizing its branches again, since the scoring
 * pass was reverted), and rounds repeat until no move improves the
 * likelihood or {@link #MAX_ROUNDS} is reached. This is a real, standard
 * (if simplified -- no fast-NNI heuristics, no simultaneous multi-move
 * rounds) topology search, not a placeholder: it can and does change the
 * tree shape away from NJ's when a better-likelihood rearrangement exists.
 * <p>
 * Restricted to small taxon counts ({@link #MAX_TAXA_FOR_NNI}) because cost
 * scales with (internal edges) x (2 candidate moves) x (rounds), each
 * requiring its own tree traversal -- unlike branch-length or
 * substitution-parameter optimization, which scale with taxa/sites alone.
 */
public final class NniTopologySearch {

    public static final int MAX_TAXA_FOR_NNI = 15;
    private static final int MAX_ROUNDS = 8;
    // A candidate must beat the current best by more than numerical noise from the per-candidate Brent
    // re-optimization (not just any epsilon > 0) or ties/near-ties would be "accepted" as real improvements.
    private static final double MIN_IMPROVEMENT = 1e-4;

    public record Result(double finalLogLikelihood, double initialLogLikelihood, int roundsUsed, int movesAccepted) {
    }

    private record Candidate(PhyloNode u, PhyloNode v, PhyloNode w, PhyloNode vChild) {
    }

    public Result search(PhyloTree tree, Map<String, String> sequencesByLabel, int alignmentLength,
                          GtrModel model, double[] gammaCategoryRates) {
        TreeLikelihoodCalculator calculator = new TreeLikelihoodCalculator();
        double initialLogL = calculator.logLikelihood(tree, sequencesByLabel, alignmentLength, model, gammaCategoryRates);
        double currentLogL = initialLogL;
        int round = 0;
        int movesAccepted = 0;

        for (; round < MAX_ROUNDS; round++) {
            Candidate bestCandidate = null;
            double bestLogL = currentLogL;

            for (Candidate candidate : collectCandidates(tree)) {
                double wBranch = candidate.w.branchLength();
                double vChildBranch = candidate.vChild.branchLength();

                applySwap(candidate);
                reoptimizeTouchedBranches(tree, candidate, sequencesByLabel, alignmentLength, model, gammaCategoryRates, calculator);
                double candidateLogL = calculator.logLikelihood(tree, sequencesByLabel, alignmentLength, model, gammaCategoryRates);

                revertSwap(candidate);
                candidate.w.setBranchLength(wBranch);
                candidate.vChild.setBranchLength(vChildBranch);

                if (candidateLogL > bestLogL + MIN_IMPROVEMENT) {
                    bestLogL = candidateLogL;
                    bestCandidate = candidate;
                }
            }

            if (bestCandidate == null) {
                break;
            }
            applySwap(bestCandidate);
            reoptimizeTouchedBranches(tree, bestCandidate, sequencesByLabel, alignmentLength, model, gammaCategoryRates, calculator);
            currentLogL = calculator.logLikelihood(tree, sequencesByLabel, alignmentLength, model, gammaCategoryRates);
            movesAccepted++;
        }

        return new Result(currentLogL, initialLogL, round, movesAccepted);
    }

    private List<Candidate> collectCandidates(PhyloTree tree) {
        List<Candidate> candidates = new ArrayList<>();
        for (PhyloNode u : tree.branchPoints()) {
            if (u.children().size() != 2) {
                continue; // defensive: every non-root internal node should be strictly bifurcating
            }
            for (PhyloNode v : new ArrayList<>(u.children())) {
                if (v.isTaxon() || v.children().size() != 2) {
                    continue; // v must itself be an internal node for this edge to have an NNI neighborhood
                }
                PhyloNode w = otherChild(u, v);
                for (PhyloNode vChild : new ArrayList<>(v.children())) {
                    candidates.add(new Candidate(u, v, w, vChild));
                }
            }
        }
        return candidates;
    }

    private PhyloNode otherChild(PhyloNode parent, PhyloNode exclude) {
        for (PhyloNode child : parent.children()) {
            if (child != exclude) {
                return child;
            }
        }
        throw new IllegalStateException("Internal node has fewer than 2 children -- not a valid NNI edge");
    }

    private void applySwap(Candidate c) {
        c.w.reparent(c.v);
        c.vChild.reparent(c.u);
    }

    private void revertSwap(Candidate c) {
        c.w.reparent(c.u);
        c.vChild.reparent(c.v);
    }

    private void reoptimizeTouchedBranches(PhyloTree tree, Candidate c, Map<String, String> sequences, int alignmentLength,
                                            GtrModel model, double[] gammaRates, TreeLikelihoodCalculator calculator) {
        // Two quick passes: each branch's optimum depends on the other's current value, so one pass each is
        // not quite a joint optimum, but is enough to give a fair likelihood signal without full-tree cost.
        for (int pass = 0; pass < 2; pass++) {
            MlBranchLengthOptimizer.optimizeSingleBranch(tree, c.w, sequences, alignmentLength, model, gammaRates, calculator);
            MlBranchLengthOptimizer.optimizeSingleBranch(tree, c.vChild, sequences, alignmentLength, model, gammaRates, calculator);
        }
    }
}
