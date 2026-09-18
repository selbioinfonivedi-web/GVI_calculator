package org.gvi.algorithms.mu;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.univariate.BrentOptimizer;
import org.apache.commons.math3.optim.univariate.SearchInterval;
import org.apache.commons.math3.optim.univariate.UnivariateObjectiveFunction;
import org.apache.commons.math3.optim.univariate.UnivariatePointValuePair;
import org.apache.commons.math3.stat.regression.SimpleRegression;
import org.gvi.algorithms.phylo.PhyloNode;
import org.gvi.algorithms.phylo.PhyloTree;
import org.gvi.core.exception.GviComputationException;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Uncorrelated lognormal relaxed molecular clock -- the same generative model BEAST2's UCLD clock
 * uses (Drummond, Ho, Phillips &amp; Rambaut 2006, "Relaxed phylogenetics and dating with
 * confidence"), fit here by maximum likelihood via coordinate ascent rather than MCMC. Every other
 * {@link MuEstimationMethod} assumes one shared rate across the whole tree (a strict clock); this
 * one instead lets each branch draw its own rate from a lognormal distribution, which matters
 * whenever different lineages plausibly evolve at different speeds (e.g. a fast-evolving outbreak
 * clade alongside slowly-diverging background diversity) -- a strict-clock fit forced onto such
 * data reports one compromise rate that describes no actual lineage well.
 * <p>
 * <b>Model.</b> Given a tree with (initially unknown) internal node dates and known tip dates, each
 * edge e has elapsed time {@code dt_e = date(child) - date(parent)} and an observed branch length
 * {@code b_e} (substitutions/site, from the same Neighbor-Joining topology {@link
 * LeastSquaresDatingEstimator} uses). Under a strict clock, every edge's implied rate
 * {@code b_e / dt_e} would be the same number; under this relaxed clock, {@code log(b_e / dt_e)} is
 * instead treated as an i.i.d. draw from {@code Normal(log(rate0), sigma^2)} -- i.e. the per-branch
 * rate is lognormally distributed with log-scale mean {@code log(rate0)} and log-scale standard
 * deviation {@code sigma}. {@code rate0} is therefore the fitted clock's <em>median</em> branch
 * rate (not the arithmetic mean -- the two coincide only when sigma=0), reported as {@link
 * MuResult#muSubPerSiteYear()} the same as every other method's single rate; {@code sigma}
 * determines the rate coefficient of variation reported alongside it, {@code sqrt(exp(sigma^2)-1)},
 * the standard summary BEAST itself reports for a UCLD clock (values noticeably above 0 indicate
 * real among-lineage rate variation; a fit that converges to sigma near 0 is behaviourally a strict
 * clock and {@link EvolutionaryRateCalculator#computeLeastSquaresDating} is the simpler, cheaper
 * estimator for that case).
 * <p>
 * <b>Fit.</b> Coordinate ascent, alternating:
 * <ol>
 *   <li>Closed-form update of {@code (rate0, sigma)}: with every node date fixed, this is exactly
 *       maximum-likelihood fitting a Normal distribution to the observed {@code log(b_e/dt_e)}
 *       values -- the sample mean and (population) standard deviation.</li>
 *   <li>Per-node date update: with {@code (rate0, sigma)} and every other node's date fixed, each
 *       internal node's date is found by 1-D maximization (Brent's method, the same pattern already
 *       used for every other ML fit in this codebase, e.g. {@code MlPhylogeneticOptimizer},
 *       {@code BdskyMlFitter}) of the log-likelihood contribution from just that node's parent edge
 *       and child edges, over the temporal-precedence-constrained interval {@code [date(parent),
 *       min(date(children))]}.</li>
 * </ol>
 * This is a real accuracy upgrade over forcing one shared rate onto genuinely heterogeneous
 * lineages, but it is a point estimate via optimization, not a full Bayesian posterior with credible
 * intervals the way BEAST's own MCMC gives -- and it still assumes the fixed NJ topology is correct,
 * exactly as {@link LeastSquaresDatingEstimator} does.
 */
public final class RelaxedClockMlEstimator {

    private static final int MAX_ITERATIONS = 200;
    private static final double CONVERGENCE_TOLERANCE = 1e-10;
    private static final double MIN_SIGMA = 1e-6;
    private static final int MAX_EVALUATIONS_PER_NODE = 200;
    /**
     * Finite stand-in for "no lower bound" on the root's own date, since Brent's method needs a
     * finite search interval. Set generously relative to the data's own timespan so it is never the
     * binding constraint in practice -- a root fit that actually wants to be this old has far bigger
     * problems than an arbitrary search bound.
     */
    private static final double ROOT_LOWER_BOUND_MARGIN_FACTOR = 100.0;
    private static final double ROOT_LOWER_BOUND_MARGIN_MIN_YEARS = 100.0;

    public record Estimate(double rate0SubPerSiteYear, double sigma, Map<PhyloNode, Double> nodeDates,
                           double logLikelihood, int iterationsUsed) {

        /** BEAST's own summary statistic for a relaxed clock's among-branch rate variation. */
        public double coefficientOfVariation() {
            return Math.sqrt(Math.exp(sigma * sigma) - 1.0);
        }
    }

    public Estimate estimate(PhyloTree tree, Map<String, Double> tipDecimalYears) {
        List<PhyloNode> allNodes = new ArrayList<>();
        collect(tree.root(), allNodes);

        List<PhyloNode> internalNodes = new ArrayList<>();
        Map<PhyloNode, Double> dates = new HashMap<>();
        for (PhyloNode n : allNodes) {
            if (n.isTaxon()) {
                Double d = tipDecimalYears.get(n.label());
                if (d == null) {
                    throw new GviComputationException("Relaxed-clock dating needs a date for every taxon; missing for '" + n.label() + "'");
                }
                dates.put(n, d);
            } else {
                internalNodes.add(n);
            }
        }
        if (internalNodes.isEmpty()) {
            throw new GviComputationException("Relaxed-clock dating needs at least one internal branch point");
        }

        double rootLowerBoundMargin = Math.max(ROOT_LOWER_BOUND_MARGIN_MIN_YEARS,
                ROOT_LOWER_BOUND_MARGIN_FACTOR * timeSpan(tipDecimalYears));

        double seedRate = seedDates(allNodes, dates, internalNodes);
        if (seedRate <= 0) {
            throw new GviComputationException("Relaxed-clock dating: initial root-to-tip rate estimate is non-positive; temporal/genetic signal too weak to seed the optimizer");
        }

        double rate0 = seedRate;
        double sigma = 0.5; // a generic, mild starting spread; the first closed-form update replaces it immediately
        double previousLogL = logLikelihood(allNodes, dates, rate0, sigma);

        int iteration = 0;
        for (; iteration < MAX_ITERATIONS; iteration++) {
            double[] fit = fitRate0AndSigma(allNodes, dates);
            rate0 = fit[0];
            sigma = fit[1];

            for (PhyloNode node : internalNodes) {
                updateNodeDate(node, dates, rate0, sigma, rootLowerBoundMargin);
            }

            double logL = logLikelihood(allNodes, dates, rate0, sigma);
            boolean converged = Math.abs(logL - previousLogL) < CONVERGENCE_TOLERANCE * Math.max(1.0, Math.abs(previousLogL));
            previousLogL = logL;
            if (converged) {
                iteration++;
                break;
            }
        }

        return new Estimate(rate0, sigma, dates, previousLogL, iteration);
    }

    private void collect(PhyloNode node, List<PhyloNode> out) {
        out.add(node);
        for (PhyloNode child : node.children()) collect(child, out);
    }

    private double timeSpan(Map<String, Double> tipDecimalYears) {
        double min = tipDecimalYears.values().stream().mapToDouble(Double::doubleValue).min().orElse(0.0);
        double max = tipDecimalYears.values().stream().mapToDouble(Double::doubleValue).max().orElse(0.0);
        return Math.max(0.0, max - min);
    }

    /** Same seeding trick as {@link LeastSquaresDatingEstimator}: root-to-tip regression, then interpolate internal dates along it. */
    private double seedDates(List<PhyloNode> allNodes, Map<PhyloNode, Double> dates, List<PhyloNode> internalNodes) {
        SimpleRegression clock = new SimpleRegression();
        for (PhyloNode n : allNodes) {
            if (n.isTaxon()) clock.addData(dates.get(n), n.rootToTipDistance());
        }
        double rate = clock.getSlope();
        double intercept = clock.getIntercept();
        double rootYear = rate > 0 ? -intercept / rate : dates.values().stream().mapToDouble(Double::doubleValue).min().orElse(0.0);

        for (PhyloNode n : internalNodes) {
            dates.put(n, rootYear + n.rootToTipDistance() / (rate > 0 ? rate : 1.0));
        }
        return rate;
    }

    /** Closed-form MLE of a Normal fit to every edge's observed log-rate: sample mean and population standard deviation. */
    private double[] fitRate0AndSigma(List<PhyloNode> allNodes, Map<PhyloNode, Double> dates) {
        List<Double> logRates = new ArrayList<>();
        for (PhyloNode n : allNodes) {
            if (n.parent() == null) continue;
            double dt = dates.get(n) - dates.get(n.parent());
            if (dt <= 0) continue; // uninformative (zero/negative-duration) edge, same guard as LeastSquaresDatingEstimator
            double observedRate = n.branchLength() / dt;
            if (observedRate <= 0) continue; // a zero branch length carries no rate information either
            logRates.add(Math.log(observedRate));
        }
        if (logRates.isEmpty()) {
            throw new GviComputationException("Relaxed-clock dating: no informative (positive-duration, positive-length) edges remain to fit a rate");
        }
        double mean = logRates.stream().mapToDouble(Double::doubleValue).average().orElse(0.0);
        double variance = logRates.stream().mapToDouble(v -> (v - mean) * (v - mean)).sum() / logRates.size();
        double sigma = Math.max(MIN_SIGMA, Math.sqrt(variance));
        return new double[]{Math.exp(mean), sigma};
    }

    /** 1-D maximization of the log-likelihood contribution touching one node, over the temporal-precedence-constrained interval. */
    private void updateNodeDate(PhyloNode node, Map<PhyloNode, Double> dates, double rate0, double sigma, double rootLowerBoundMargin) {
        PhyloNode parent = node.parent();
        List<PhyloNode> children = node.children();

        double upperBound = children.stream().mapToDouble(dates::get).min().orElse(Double.POSITIVE_INFINITY);
        double lowerBound = parent != null ? dates.get(parent) : upperBound - rootLowerBoundMargin;
        if (!(lowerBound < upperBound)) return; // numerically pinned this round; leave the date as-is rather than search an empty/inverted interval

        // Shrink the search interval a hair inside its true bounds: at either endpoint some touching
        // edge's elapsed time hits exactly zero, whose log-rate is -Infinity -- a valid but useless
        // evaluation point for Brent to waste evaluations probing.
        double epsilon = (upperBound - lowerBound) * 1e-9;
        double searchLower = lowerBound + epsilon;
        double searchUpper = upperBound - epsilon;
        if (!(searchLower < searchUpper)) return;

        UnivariateFunction objective = t -> localLogLikelihood(node, parent, children, dates, t, rate0, sigma);
        double currentDate = Math.max(searchLower, Math.min(searchUpper, dates.get(node)));

        BrentOptimizer optimizer = new BrentOptimizer(1e-8, 1e-12);
        UnivariatePointValuePair result = optimizer.optimize(
                new MaxEval(MAX_EVALUATIONS_PER_NODE),
                new UnivariateObjectiveFunction(objective),
                GoalType.MAXIMIZE,
                new SearchInterval(searchLower, searchUpper, currentDate));
        dates.put(node, result.getPoint());
    }

    private double localLogLikelihood(PhyloNode node, PhyloNode parent, List<PhyloNode> children,
                                       Map<PhyloNode, Double> dates, double candidateDate, double rate0, double sigma) {
        double logL = 0.0;
        if (parent != null) {
            double dt = candidateDate - dates.get(parent);
            logL += edgeLogLikelihood(node.branchLength(), dt, rate0, sigma);
        }
        for (PhyloNode child : children) {
            double dt = dates.get(child) - candidateDate;
            logL += edgeLogLikelihood(child.branchLength(), dt, rate0, sigma);
        }
        return logL;
    }

    private double logLikelihood(List<PhyloNode> allNodes, Map<PhyloNode, Double> dates, double rate0, double sigma) {
        double logL = 0.0;
        for (PhyloNode n : allNodes) {
            if (n.parent() == null) continue;
            double dt = dates.get(n) - dates.get(n.parent());
            logL += edgeLogLikelihood(n.branchLength(), dt, rate0, sigma);
        }
        return logL;
    }

    /**
     * log-density of one edge's observed rate under Normal(log(rate0), sigma^2) in log-rate space.
     * <p>
     * A degenerate edge (zero/negative elapsed time, or a zero branch length -- NJ clamps negative
     * branch lengths to exactly zero, a documented, real occurrence, not a rare edge case) carries no
     * rate information rather than being "impossible": it contributes 0 (i.e. is skipped), the same
     * treatment {@link #fitRate0AndSigma} already gives such edges. Returning {@code
     * Double.NEGATIVE_INFINITY} instead would be wrong twice over -- it would permanently poison
     * every downstream sum this feeds (the coordinate-ascent convergence check compares two
     * log-likelihoods by subtraction, and {@code -Infinity - (-Infinity) = NaN}, so {@code
     * Math.abs(...) < tolerance} is never true and the fit silently runs to {@link #MAX_ITERATIONS}
     * every time instead of converging), and it would drown out a touching node's OTHER, genuinely
     * informative edge in the same local sum ({@code -Infinity + anything = -Infinity}).
     */
    private double edgeLogLikelihood(double branchLength, double dt, double rate0, double sigma) {
        if (dt <= 0 || branchLength <= 0) return 0.0;
        double logRate = Math.log(branchLength / dt);
        double z = (logRate - Math.log(rate0)) / sigma;
        return -Math.log(sigma) - 0.5 * Math.log(2 * Math.PI) - 0.5 * z * z;
    }
}
