package org.gvi.algorithms.mu;

import org.apache.commons.math3.stat.regression.SimpleRegression;
import org.gvi.algorithms.phylo.PhyloNode;
import org.gvi.algorithms.phylo.PhyloTree;
import org.gvi.core.exception.GviComputationException;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Native least-squares divergence-time dating -- the same estimation
 * problem To, Jung, Ly-Trong, Minh &amp; von Haeseler (2016, "Fast dating
 * using least-squares") solve as the LSD2 algorithm behind IQ-TREE's
 * {@code --date} option: given a fixed tree topology with branch lengths
 * in substitution units and known tip collection dates, jointly estimate a
 * single strict-clock rate mu AND every internal node's divergence date,
 * minimizing the total squared error between each edge's observed branch
 * length and {@code mu * elapsedYears} on that edge, subject to every
 * node's date being no later than any of its children's (temporal
 * precedence along the tree).
 * <p>
 * This is a real accuracy upgrade over {@link EvolutionaryRateCalculator}'s
 * simpler root-to-tip regression: root-to-tip regression fits ONE straight
 * line through (date, cumulative-root-to-tip-distance) pairs, implicitly
 * treating every tip as an independent, equally-weighted observation of
 * the same clock -- it does not use each individual EDGE's own length,
 * and does not explicitly enforce that internal node dates respect the
 * tree's own branching order (though in practice a fitted straight line
 * rarely violates it badly). This estimator instead uses every edge
 * directly and enforces temporal precedence explicitly, the same problem
 * LSD2 solves -- though via a different numerical method: LSD2 uses a
 * specialized O(n) algorithm; this uses constrained coordinate ascent
 * (closed-form per step, see below), the same general pattern already used
 * elsewhere in this codebase (e.g. {@code MlPhylogeneticOptimizer}) rather
 * than LSD2's own internal algorithm. Both solve the identical objective,
 * so both should converge to the same (or a very close) answer on
 * well-behaved data.
 * <p>
 * Each coordinate-ascent step has a closed-form solution (derived directly
 * by setting the relevant partial derivative of the sum-of-squares
 * objective to zero -- no black-box general optimizer needed):
 * <ul>
 *   <li>Given fixed node dates, the rate mu that minimizes
 *       {@code sum((branchLength - mu*elapsed)^2)} over all edges is the
 *       weighted-least-squares slope-through-the-origin,
 *       {@code mu = sum(b_e * dt_e) / sum(dt_e^2)}.</li>
 *   <li>Given a fixed rate and every other node's date, a single internal
 *       node x's optimal date (minimizing just the sum-of-squares terms
 *       touching x: its one parent edge and each of its child edges) is
 *       {@code t_x = (b_px + mu*t_p - sum(b_xc) + mu*sum(t_c)) / (mu*(1+childCount))}
 *       for a non-root node, or the corresponding parent-free form for the
 *       root -- both derived by differentiating the local sum of squares
 *       and solving the resulting linear equation in t_x.</li>
 * </ul>
 * After each closed-form update, the node's date is clamped to
 * {@code [parent's current date, min(children's current dates)]} to
 * enforce temporal precedence (a standard projected-coordinate-descent
 * technique for a constrained convex problem).
 */
public final class LeastSquaresDatingEstimator {

    private static final int MAX_ITERATIONS = 500;
    private static final double CONVERGENCE_TOLERANCE = 1e-12;

    public record Estimate(double muSubPerSiteYear, Map<PhyloNode, Double> nodeDates, double weightedSse, int iterationsUsed) {
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
                    throw new GviComputationException("Least-squares dating needs a date for every taxon; missing for '" + n.label() + "'");
                }
                dates.put(n, d);
            } else {
                internalNodes.add(n);
            }
        }
        if (internalNodes.isEmpty()) {
            throw new GviComputationException("Least-squares dating needs at least one internal branch point");
        }

        double mu = initialize(allNodes, dates, internalNodes);
        if (mu <= 0) {
            throw new GviComputationException("Least-squares dating: initial root-to-tip rate estimate is non-positive; temporal/genetic signal too weak to seed the optimizer");
        }

        double previousSse = weightedSse(allNodes, dates, mu);
        int iteration = 0;
        for (; iteration < MAX_ITERATIONS; iteration++) {
            mu = updateMu(allNodes, dates);
            for (PhyloNode node : internalNodes) {
                updateNodeDate(node, dates, mu);
            }
            double sse = weightedSse(allNodes, dates, mu);
            boolean converged = Math.abs(previousSse - sse) < CONVERGENCE_TOLERANCE * Math.max(1.0, previousSse);
            previousSse = sse;
            if (converged) {
                iteration++;
                break;
            }
        }

        return new Estimate(mu, dates, previousSse, iteration);
    }

    private void collect(PhyloNode node, List<PhyloNode> out) {
        out.add(node);
        for (PhyloNode child : node.children()) collect(child, out);
    }

    /** Seeds mu via ordinary root-to-tip regression, then interpolates every internal date along its (monotonic, since NJ clamps negative branch lengths to 0) root-to-tip distance -- a feasible starting point for the constrained optimization that follows. */
    private double initialize(List<PhyloNode> allNodes, Map<PhyloNode, Double> dates, List<PhyloNode> internalNodes) {
        SimpleRegression clock = new SimpleRegression();
        for (PhyloNode n : allNodes) {
            if (n.isTaxon()) clock.addData(dates.get(n), n.rootToTipDistance());
        }
        double mu = clock.getSlope();
        double intercept = clock.getIntercept();
        double rootYear = mu > 0 ? -intercept / mu : dates.values().stream().mapToDouble(Double::doubleValue).min().orElse(0.0);

        for (PhyloNode n : internalNodes) {
            dates.put(n, rootYear + n.rootToTipDistance() / (mu > 0 ? mu : 1.0));
        }
        return mu;
    }

    /** Closed-form weighted-least-squares rate: minimizes sum((b_e - mu*dt_e)^2) over all edges. */
    private double updateMu(List<PhyloNode> allNodes, Map<PhyloNode, Double> dates) {
        double numerator = 0.0, denominator = 0.0;
        for (PhyloNode n : allNodes) {
            if (n.parent() == null) continue;
            double dt = dates.get(n) - dates.get(n.parent());
            if (dt <= 0) continue; // a zero-length interval carries no rate information; skip rather than divide by zero
            numerator += n.branchLength() * dt;
            denominator += dt * dt;
        }
        if (denominator <= 0) {
            throw new GviComputationException("Least-squares dating: no informative (positive-duration) edges remain to fit a rate");
        }
        return numerator / denominator;
    }

    /** Closed-form date update for one internal node, clamped to respect temporal precedence against its (currently fixed) parent and children. */
    private void updateNodeDate(PhyloNode node, Map<PhyloNode, Double> dates, double mu) {
        PhyloNode parent = node.parent();
        List<PhyloNode> children = node.children();

        double weightSum = children.size() + (parent != null ? 1 : 0);
        if (weightSum == 0 || mu <= 0) return;

        double rhs = 0.0;
        if (parent != null) {
            rhs += node.branchLength() + mu * dates.get(parent);
        }
        for (PhyloNode child : children) {
            rhs += mu * dates.get(child) - child.branchLength();
        }
        double proposed = rhs / (mu * weightSum);

        double lowerBound = parent != null ? dates.get(parent) : Double.NEGATIVE_INFINITY;
        double upperBound = children.stream().mapToDouble(dates::get).min().orElse(Double.POSITIVE_INFINITY);
        double clamped = Math.max(lowerBound, Math.min(upperBound, proposed));
        dates.put(node, clamped);
    }

    private double weightedSse(List<PhyloNode> allNodes, Map<PhyloNode, Double> dates, double mu) {
        double sse = 0.0;
        for (PhyloNode n : allNodes) {
            if (n.parent() == null) continue;
            double dt = dates.get(n) - dates.get(n.parent());
            double residual = n.branchLength() - mu * dt;
            sse += residual * residual;
        }
        return sse;
    }
}
