package org.gvi.algorithms.phylo.model;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.univariate.BrentOptimizer;
import org.apache.commons.math3.optim.univariate.SearchInterval;
import org.apache.commons.math3.optim.univariate.UnivariateObjectiveFunction;
import org.apache.commons.math3.optim.univariate.UnivariatePointValuePair;
import org.gvi.algorithms.phylo.PhyloNode;
import org.gvi.algorithms.phylo.PhyloTree;

import java.util.List;
import java.util.Map;
import java.util.Random;

/**
 * Joint maximum-likelihood optimization of branch lengths AND substitution
 * model parameters (the rate parameters of a {@link RateParameterization}
 * -- e.g. kappa for HKY85, the 5 free GTR rates -- plus, if Gamma rate
 * heterogeneity is enabled, the alpha shape parameter) via coordinate
 * ascent: cycle through (all branch lengths, each rate parameter, alpha),
 * 1D Brent-optimizing each while holding everything else fixed, repeating
 * until the total log-likelihood stops improving.
 * <p>
 * This is what makes the rate parameters and alpha genuinely ML-estimated
 * rather than empirical/user-supplied (the gap flagged in
 * {@link MlBranchLengthOptimizer} and {@link EvolutionaryRateCalculator}'s
 * earlier version) -- coordinate ascent is not a full joint multivariate
 * optimization (it can't move diagonally in parameter space in one step),
 * but it provably never decreases the likelihood at each step and is what
 * many real ML phylogenetics tools actually do in practice (alternating
 * branch-length and model-parameter optimization phases).
 */
public final class MlPhylogeneticOptimizer {

    private static final double CONVERGENCE_THRESHOLD = 1e-4;
    private static final int MAX_CYCLES = 15;
    private static final int MAX_EVALUATIONS_PER_PARAMETER = 80;
    private static final double MIN_ALPHA = 0.02;
    private static final double MAX_ALPHA = 50.0;

    public record Result(double finalLogLikelihood, double initialLogLikelihood, int cyclesUsed,
                          double[] optimizedTheta, Double optimizedAlpha, GtrModel finalModel, double[] finalGammaRates) {
    }

    /**
     * @param gammaCategoryCount 0 or 1 to disable Gamma rate heterogeneity (single rate=1.0, alpha not optimized);
     *                           >=2 to enable it with that many discretized categories, jointly ML-optimizing alpha
     * @param initialAlpha       starting value for alpha when Gamma is enabled (ignored otherwise)
     */
    public Result optimize(PhyloTree tree, Map<String, String> sequencesByLabel, int alignmentLength,
                            RateParameterization parameterization, double[] empiricalFrequencies,
                            int gammaCategoryCount, double initialAlpha) {
        TreeLikelihoodCalculator calculator = new TreeLikelihoodCalculator();
        List<PhyloNode> branches = MlBranchLengthOptimizer.collectBranches(tree);
        boolean gammaEnabled = gammaCategoryCount >= 2;

        double[] theta = parameterization.initialTheta().clone();
        double[] alphaHolder = {initialAlpha};

        GtrModel model = parameterization.buildModel(theta, empiricalFrequencies);
        double[] gammaRates = gammaEnabled ? DiscreteGammaRates.categories(alphaHolder[0], gammaCategoryCount) : new double[]{1.0};

        double initialLogL = calculator.logLikelihood(tree, sequencesByLabel, alignmentLength, model, gammaRates);
        double previousLogL = initialLogL;
        int cycle = 0;

        for (; cycle < MAX_CYCLES; cycle++) {
            MlBranchLengthOptimizer.oneCycleOverAllBranches(tree, branches, sequencesByLabel, alignmentLength, model, gammaRates, calculator);

            for (int i = 0; i < theta.length; i++) {
                final int idx = i;
                final double[] finalGammaRates = gammaRates;
                UnivariateFunction objective = value -> {
                    theta[idx] = value;
                    GtrModel candidate = parameterization.buildModel(theta, empiricalFrequencies);
                    return calculator.logLikelihood(tree, sequencesByLabel, alignmentLength, candidate, finalGammaRates);
                };
                theta[idx] = brentMaximize(objective, parameterization.lowerBounds()[i], parameterization.upperBounds()[i], theta[idx]);
            }
            model = parameterization.buildModel(theta, empiricalFrequencies);

            if (gammaEnabled) {
                final GtrModel currentModel = model;
                UnivariateFunction objective = a -> {
                    double[] rates = DiscreteGammaRates.categories(a, gammaCategoryCount);
                    return calculator.logLikelihood(tree, sequencesByLabel, alignmentLength, currentModel, rates);
                };
                alphaHolder[0] = brentMaximize(objective, MIN_ALPHA, MAX_ALPHA, alphaHolder[0]);
                gammaRates = DiscreteGammaRates.categories(alphaHolder[0], gammaCategoryCount);
            }

            double newLogL = calculator.logLikelihood(tree, sequencesByLabel, alignmentLength, model, gammaRates);
            boolean converged = (newLogL - previousLogL) < CONVERGENCE_THRESHOLD;
            previousLogL = newLogL;
            if (converged) {
                cycle++;
                break;
            }
        }

        return new Result(previousLogL, initialLogL, cycle, theta, gammaEnabled ? alphaHolder[0] : null, model, gammaRates);
    }

    /**
     * Same joint optimization as {@link #optimize}, but coordinate ascent only
     * guarantees convergence to a <em>local</em> optimum -- from a bad starting
     * point it can settle short of the true maximum-likelihood parameters. This
     * runs the default start plus {@code additionalRandomRestarts} more from
     * randomized starting points (rate parameters and, if enabled, alpha drawn
     * uniformly from their bounds; branch lengths reset to the tree's original
     * lengths each time so restarts don't inherit a previous attempt's
     * branch-length state), and keeps whichever run reached the highest final
     * log-likelihood. This does not guarantee the global optimum either -- no
     * finite number of restarts can -- but it materially reduces the chance
     * that a single unlucky starting point is mistaken for "the" ML estimate.
     */
    public Result optimizeWithRestarts(PhyloTree tree, Map<String, String> sequencesByLabel, int alignmentLength,
                                        RateParameterization parameterization, double[] empiricalFrequencies,
                                        int gammaCategoryCount, double initialAlpha,
                                        int additionalRandomRestarts, long seed) {
        List<PhyloNode> branches = MlBranchLengthOptimizer.collectBranches(tree);
        double[] originalBranchLengths = snapshotBranchLengths(branches);

        Result best = optimize(tree, sequencesByLabel, alignmentLength, parameterization, empiricalFrequencies, gammaCategoryCount, initialAlpha);
        double[] bestBranchLengths = snapshotBranchLengths(branches);

        Random random = new Random(seed);
        boolean gammaEnabled = gammaCategoryCount >= 2;
        for (int r = 0; r < additionalRandomRestarts; r++) {
            restoreBranchLengths(branches, originalBranchLengths);
            double[] randomTheta = new double[parameterization.initialTheta().length];
            for (int i = 0; i < randomTheta.length; i++) {
                double lo = parameterization.lowerBounds()[i];
                double hi = parameterization.upperBounds()[i];
                randomTheta[i] = lo + random.nextDouble() * (hi - lo);
            }
            RateParameterization restartParameterization = new RateParameterization(parameterization.modelName(),
                    parameterization.useEmpiricalFrequencies(), randomTheta, parameterization.lowerBounds(),
                    parameterization.upperBounds(), parameterization.thetaToRates());
            double randomAlpha = gammaEnabled ? MIN_ALPHA + random.nextDouble() * (MAX_ALPHA - MIN_ALPHA) : initialAlpha;

            Result candidate = optimize(tree, sequencesByLabel, alignmentLength, restartParameterization,
                    empiricalFrequencies, gammaCategoryCount, randomAlpha);
            if (candidate.finalLogLikelihood() > best.finalLogLikelihood()) {
                best = candidate;
                bestBranchLengths = snapshotBranchLengths(branches);
            }
        }

        restoreBranchLengths(branches, bestBranchLengths);
        return best;
    }

    private double[] snapshotBranchLengths(List<PhyloNode> branches) {
        double[] snapshot = new double[branches.size()];
        for (int i = 0; i < branches.size(); i++) snapshot[i] = branches.get(i).branchLength();
        return snapshot;
    }

    private void restoreBranchLengths(List<PhyloNode> branches, double[] snapshot) {
        for (int i = 0; i < branches.size(); i++) branches.get(i).setBranchLength(snapshot[i]);
    }

    private double brentMaximize(UnivariateFunction objective, double lower, double upper, double startingPoint) {
        BrentOptimizer optimizer = new BrentOptimizer(1e-5, 1e-10);
        double clampedStart = Math.max(lower, Math.min(upper, startingPoint));
        UnivariatePointValuePair result = optimizer.optimize(
                new MaxEval(MAX_EVALUATIONS_PER_PARAMETER),
                new UnivariateObjectiveFunction(objective),
                GoalType.MAXIMIZE,
                new SearchInterval(lower, upper, clampedStart));
        return result.getPoint();
    }
}
