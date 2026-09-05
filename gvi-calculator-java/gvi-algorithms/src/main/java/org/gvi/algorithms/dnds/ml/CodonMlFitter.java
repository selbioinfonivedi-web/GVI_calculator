package org.gvi.algorithms.dnds.ml;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.univariate.BrentOptimizer;
import org.apache.commons.math3.optim.univariate.SearchInterval;
import org.apache.commons.math3.optim.univariate.UnivariateObjectiveFunction;
import org.apache.commons.math3.optim.univariate.UnivariatePointValuePair;
import org.gvi.algorithms.phylo.PhyloNode;
import org.gvi.algorithms.phylo.PhyloTree;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * Joint maximum-likelihood fit of branch lengths, kappa, and omega under
 * {@link CodonModel} -- the same coordinate-ascent pattern
 * {@link org.gvi.algorithms.phylo.model.MlPhylogeneticOptimizer} uses for
 * nucleotide GTR(+Gamma): cycle through (every branch length, kappa,
 * omega), 1D Brent-optimizing each while holding the rest fixed, repeat
 * until the log-likelihood stops improving. This is what makes omega
 * (dN/dS) a genuine maximum-likelihood estimate under an explicit codon
 * substitution process -- {@code codeml}'s M0 "one-ratio" model -- instead
 * of the counting-based Nei-Gojobori statistic used elsewhere in this
 * project. Coordinate ascent provably never decreases the likelihood per
 * step but is not a joint multivariate optimizer, so (as with the
 * nucleotide optimizer) it can still settle into a local optimum.
 */
final class CodonMlFitter {

    private static final double MIN_BRANCH_LENGTH = 1e-8;
    private static final double MAX_BRANCH_LENGTH = 10.0;
    private static final double MIN_KAPPA = 0.05;
    private static final double MAX_KAPPA = 20.0;
    private static final double MIN_OMEGA = 1e-3;
    private static final double MAX_OMEGA = 20.0;
    private static final double CONVERGENCE_THRESHOLD = 1e-3;
    private static final int MAX_CYCLES = 10;
    private static final int MAX_EVALUATIONS_PER_PARAMETER = 40;

    record Result(double omega, double kappa, double logLikelihood, double initialLogLikelihood, int cyclesUsed) {
    }

    Result optimize(PhyloTree tree, Map<String, List<String>> codonsByLabel, int siteCount, double[] pi,
                     double initialKappa, double initialOmega) {
        CodonTreeLikelihood calculator = new CodonTreeLikelihood();
        List<PhyloNode> branches = collectBranches(tree);

        double[] kappaHolder = {initialKappa};
        double[] omegaHolder = {initialOmega};
        CodonModel[] modelHolder = {CodonModel.of(pi, kappaHolder[0], omegaHolder[0])};

        double initialLogL = calculator.logLikelihood(tree, codonsByLabel, siteCount, modelHolder[0]);
        double previousLogL = initialLogL;
        int cycle = 0;

        for (; cycle < MAX_CYCLES; cycle++) {
            // Built once per cycle (the model is fixed for the whole branch-length round); each branch's own
            // Brent search below only replaces THAT branch's entry per candidate length, not the whole cache --
            // see CodonTreeLikelihood's javadoc for why that distinction matters for performance at 61 states.
            final CodonModel cycleModel = modelHolder[0];
            Map<PhyloNode, double[][]> matrices = calculator.buildTransitionMatrices(tree, cycleModel);
            double[] freq = cycleModel.frequencies();
            for (PhyloNode branch : branches) {
                UnivariateFunction objective = t -> {
                    double clamped = Math.max(MIN_BRANCH_LENGTH, Math.min(MAX_BRANCH_LENGTH, t));
                    branch.setBranchLength(clamped);
                    matrices.put(branch, cycleModel.transitionProbabilities(clamped));
                    return calculator.logLikelihoodWithMatrices(tree, codonsByLabel, siteCount, matrices, freq);
                };
                double best = brentMaximize(objective, MIN_BRANCH_LENGTH, MAX_BRANCH_LENGTH, branch.branchLength());
                branch.setBranchLength(best);
                matrices.put(branch, cycleModel.transitionProbabilities(best));
            }

            UnivariateFunction kappaObjective = k -> {
                CodonModel candidate = CodonModel.of(pi, k, omegaHolder[0]);
                return calculator.logLikelihood(tree, codonsByLabel, siteCount, candidate);
            };
            kappaHolder[0] = brentMaximize(kappaObjective, MIN_KAPPA, MAX_KAPPA, kappaHolder[0]);

            UnivariateFunction omegaObjective = w -> {
                CodonModel candidate = CodonModel.of(pi, kappaHolder[0], w);
                return calculator.logLikelihood(tree, codonsByLabel, siteCount, candidate);
            };
            omegaHolder[0] = brentMaximize(omegaObjective, MIN_OMEGA, MAX_OMEGA, omegaHolder[0]);

            modelHolder[0] = CodonModel.of(pi, kappaHolder[0], omegaHolder[0]);
            double newLogL = calculator.logLikelihood(tree, codonsByLabel, siteCount, modelHolder[0]);
            boolean converged = (newLogL - previousLogL) < CONVERGENCE_THRESHOLD;
            previousLogL = newLogL;
            if (converged) {
                cycle++;
                break;
            }
        }

        return new Result(omegaHolder[0], kappaHolder[0], previousLogL, initialLogL, cycle);
    }

    private double brentMaximize(UnivariateFunction objective, double lower, double upper, double startingPoint) {
        BrentOptimizer optimizer = new BrentOptimizer(1e-4, 1e-9);
        double clampedStart = Math.max(lower, Math.min(upper, startingPoint));
        UnivariatePointValuePair result = optimizer.optimize(
                new MaxEval(MAX_EVALUATIONS_PER_PARAMETER),
                new UnivariateObjectiveFunction(objective),
                GoalType.MAXIMIZE,
                new SearchInterval(lower, upper, clampedStart));
        return result.getPoint();
    }

    private List<PhyloNode> collectBranches(PhyloTree tree) {
        List<PhyloNode> nodes = new ArrayList<>();
        collect(tree.root(), nodes);
        return nodes;
    }

    private void collect(PhyloNode node, List<PhyloNode> out) {
        for (PhyloNode child : node.children()) {
            out.add(child);
            collect(child, out);
        }
    }
}
