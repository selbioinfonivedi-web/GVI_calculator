package org.gvi.algorithms.phylo.model;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.univariate.BrentOptimizer;
import org.apache.commons.math3.optim.univariate.SearchInterval;
import org.apache.commons.math3.optim.univariate.UnivariateObjectiveFunction;
import org.apache.commons.math3.optim.univariate.UnivariatePointValuePair;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.gvi.algorithms.phylo.PhyloNode;
import org.gvi.algorithms.phylo.PhyloTree;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * Maximum-likelihood branch-length optimization on a FIXED topology
 * (Felsenstein pruning likelihood as the objective, GTR+Gamma as the
 * model): cycles through every branch, 1D-optimizing its length via
 * Brent's method while holding all others fixed, repeating until the
 * total log-likelihood stops improving meaningfully. This is the standard
 * "branch length optimization" step real ML tools interleave with
 * topology search -- here used on its own, on the Neighbor-Joining
 * starting topology, as a real accuracy upgrade over NJ's plain
 * distance-based branch lengths (which don't jointly account for the rest
 * of the tree the way a likelihood-based estimate does).
 */
public final class MlBranchLengthOptimizer {

    private static final double MIN_BRANCH_LENGTH = 1e-8;
    private static final double MAX_BRANCH_LENGTH = 10.0;
    private static final double CONVERGENCE_THRESHOLD = 1e-4;
    private static final int MAX_CYCLES = 20;
    private static final int MAX_EVALUATIONS_PER_BRANCH = 100;

    public record Result(double finalLogLikelihood, double initialLogLikelihood, int cyclesUsed) {
    }

    public Result optimize(PhyloTree tree, Map<String, String> sequencesByLabel, int alignmentLength,
                            GtrModel model, double[] gammaCategoryRates) {
        TreeLikelihoodCalculator calculator = new TreeLikelihoodCalculator();
        List<PhyloNode> branches = collectBranches(tree);

        double initialLogL = calculator.logLikelihood(tree, sequencesByLabel, alignmentLength, model, gammaCategoryRates);
        double previousLogL = initialLogL;
        int cycle = 0;

        for (; cycle < MAX_CYCLES; cycle++) {
            oneCycleOverAllBranches(tree, branches, sequencesByLabel, alignmentLength, model, gammaCategoryRates, calculator);
            double newLogL = calculator.logLikelihood(tree, sequencesByLabel, alignmentLength, model, gammaCategoryRates);
            boolean converged = (newLogL - previousLogL) < CONVERGENCE_THRESHOLD;
            previousLogL = newLogL;
            if (converged) {
                cycle++;
                break;
            }
        }

        return new Result(previousLogL, initialLogL, cycle);
    }

    /** One Brent-optimization pass over every branch (package-visible so {@link MlPhylogeneticOptimizer} can interleave it with substitution-parameter optimization). */
    static void oneCycleOverAllBranches(PhyloTree tree, List<PhyloNode> branches, Map<String, String> sequences, int alignmentLength,
                                         GtrModel model, double[] gammaRates, TreeLikelihoodCalculator calculator) {
        for (PhyloNode branch : branches) {
            optimizeSingleBranch(tree, branch, sequences, alignmentLength, model, gammaRates, calculator);
        }
    }

    /** Package-visible so {@link NniTopologySearch} can re-optimize just the handful of branches a candidate move touched, without paying for a full-tree cycle. */
    static void optimizeSingleBranch(PhyloTree tree, PhyloNode branch, Map<String, String> sequences, int alignmentLength,
                                      GtrModel model, double[] gammaRates, TreeLikelihoodCalculator calculator) {
        UnivariateFunction objective = t -> {
            branch.setBranchLength(Math.max(MIN_BRANCH_LENGTH, Math.min(MAX_BRANCH_LENGTH, t)));
            return calculator.logLikelihood(tree, sequences, alignmentLength, model, gammaRates);
        };

        BrentOptimizer optimizer = new BrentOptimizer(1e-5, 1e-10);
        double startingPoint = Math.max(MIN_BRANCH_LENGTH, Math.min(MAX_BRANCH_LENGTH, branch.branchLength()));
        UnivariatePointValuePair result = optimizer.optimize(
                new MaxEval(MAX_EVALUATIONS_PER_BRANCH),
                new UnivariateObjectiveFunction(objective),
                GoalType.MAXIMIZE,
                new SearchInterval(MIN_BRANCH_LENGTH, MAX_BRANCH_LENGTH, startingPoint));
        branch.setBranchLength(result.getPoint());
    }

    /** Every node except the root itself (the root taxon's own "branch" is represented by its child's branchLength -- see TreeLikelihoodCalculator). */
    static List<PhyloNode> collectBranches(PhyloTree tree) {
        List<PhyloNode> nodes = new ArrayList<>();
        collect(tree.root(), nodes);
        return nodes;
    }

    private static void collect(PhyloNode node, List<PhyloNode> out) {
        for (PhyloNode child : node.children()) {
            out.add(child);
            collect(child, out);
        }
    }
}
