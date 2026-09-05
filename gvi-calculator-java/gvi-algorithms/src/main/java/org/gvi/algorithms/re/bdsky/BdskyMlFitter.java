package org.gvi.algorithms.re.bdsky;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.univariate.BrentOptimizer;
import org.apache.commons.math3.optim.univariate.SearchInterval;
import org.apache.commons.math3.optim.univariate.UnivariateObjectiveFunction;
import org.apache.commons.math3.optim.univariate.UnivariatePointValuePair;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.gvi.algorithms.phylo.PhyloNode;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * Maximum-likelihood fit of the constant-rate birth-death-sampling process
 * ({@link BdskyTreeLikelihood}) via coordinate ascent -- the same general
 * pattern already used elsewhere in this codebase for ML phylogenetics
 * ({@code MlPhylogeneticOptimizer}: cycle through parameters, 1D-optimize
 * each via Brent's method holding the others fixed, repeat to
 * convergence). This is a real accuracy upgrade over the phylodynamic
 * fallback's Euler-Lotka growth-rate conversion -- it uses the actual
 * birth-death-sampling tree likelihood (the same generative model BEAST2's
 * BDSKY package maximizes via MCMC) rather than a summary growth rate --
 * but it is a point estimate via optimization, not a full Bayesian
 * posterior with credible intervals the way BDSKY's own MCMC gives.
 * <p>
 * Parameterized as (Re, delta, s) -- reproductive number, become-
 * uninfectious rate, sampling proportion -- matching BEAST2's own BDSKY
 * parameterization (not the raw (lambda, mu, psi) rates) since it makes
 * every parameter's natural bound explicit and directly gives the
 * quantity actually wanted: {@code lambda = Re*delta, psi = s*delta,
 * mu = (1-s)*delta}.
 */
final class BdskyMlFitter {

    private static final int MAX_CYCLES = 50;
    private static final double CONVERGENCE_TOLERANCE = 1e-8;
    private static final int MAX_EVALUATIONS_PER_PARAMETER = 200;

    static final double RE_MIN = 1e-3, RE_MAX = 20.0;
    /**
     * How close to a search bound counts as "pinned to it", expressed as a fraction of the
     * search interval's width.
     * <p>
     * A maximum-likelihood fit that comes to rest against the edge of its own search interval
     * has not found an interior optimum -- the likelihood was still increasing when the
     * optimizer ran out of room, so the value reported is a property of the bracket, not of
     * the data. Callers must treat such a fit as unidentifiable rather than as an estimate.
     * <p>
     * The tolerance is interval-relative rather than machine-epsilon because Brent's method
     * converges to <em>near</em> a bound, not exactly onto it: a real fit on this project's
     * corpus stopped at 19.99997674 against a bound of 20.0, which any tolerance tighter than
     * ~2.3e-5 would have waved through as a genuine estimate of Re = 20.
     */
    private static final double BOUND_TOLERANCE_FRACTION = 1e-4;

    /** True when {@code re} came to rest against either edge of the search interval. */
    public static boolean reIsAtSearchBound(double re) {
        double margin = (RE_MAX - RE_MIN) * BOUND_TOLERANCE_FRACTION;
        return (re - RE_MIN) <= margin || (RE_MAX - re) <= margin;
    }
    private static final double DELTA_MIN = 1e-3, DELTA_MAX = 2000.0;
    private static final double S_MIN = 1e-6, S_MAX = 1.0 - 1e-6;

    record Result(double re, double becomeUninfectiousRate, double samplingProportion, double logLikelihood, int cyclesUsed) {
    }

    Result fit(PhyloNode root, Map<PhyloNode, Double> dates, double initialRe, double initialDelta, double initialS) {
        double re = clamp(initialRe, RE_MIN, RE_MAX);
        double delta = clamp(initialDelta, DELTA_MIN, DELTA_MAX);
        double s = clamp(initialS, S_MIN, S_MAX);
        double previousLogL = logLikelihoodAt(root, dates, re, delta, s);

        int cycle = 0;
        for (; cycle < MAX_CYCLES; cycle++) {
            double fixedDelta = delta, fixedS = s;
            re = brentMaximize(r -> logLikelihoodAt(root, dates, r, fixedDelta, fixedS), RE_MIN, RE_MAX, re);
            double fixedRe = re;
            delta = brentMaximize(d -> logLikelihoodAt(root, dates, fixedRe, d, fixedS), DELTA_MIN, DELTA_MAX, delta);
            double fixedDelta2 = delta;
            s = brentMaximize(sv -> logLikelihoodAt(root, dates, fixedRe, fixedDelta2, sv), S_MIN, S_MAX, s);

            double logL = logLikelihoodAt(root, dates, re, delta, s);
            boolean converged = Math.abs(logL - previousLogL) < CONVERGENCE_TOLERANCE * Math.max(1.0, Math.abs(previousLogL));
            previousLogL = logL;
            if (converged) {
                cycle++;
                break;
            }
        }
        return new Result(re, delta, s, previousLogL, cycle);
    }




    /** One point on the Re profile: the best log-likelihood attainable at that Re. */
    record ProfilePoint(double re, double samplingProportion, double logLikelihood) {
    }

    /**
     * Profile likelihood over Re with the become-uninfectious rate held at the value implied by the
     * caller's generation time, maximizing over the sampling proportion at each Re.
     * <p>
     * This replaces a free three-parameter coordinate ascent over (Re, delta, s), which did not work:
     * the likelihood has a ridge along which delta and s trade off against Re at constant
     * log-likelihood, and axis-parallel coordinate ascent slides along it to wherever it happens to
     * start. Measured on a simulated tree with a true Re of 1.5, four different starting points all
     * converged to a log-likelihood of -31.72258 while reporting Re of 1.06, 1.15, 1.50 and 1.75 --
     * the number returned was a property of the initial guess, not of the data.
     * <p>
     * Pinning delta is not a modelling convenience; it is the one piece of the triple actually known
     * from outside the tree. delta is the rate of becoming uninfectious, i.e. the reciprocal of the
     * infectious period, which is what the bundled generation-time table supplies. Leaving it free
     * let a fit initialized at 73/yr wander to 10/yr on data generated at 0.8/yr.
     * <p>
     * The profile is reported with a likelihood-ratio interval rather than as a bare point estimate,
     * because on realistic tree sizes it is wide -- honestly so.
     */
    List<ProfilePoint> profileOverRe(PhyloNode root, Map<PhyloNode, Double> dates, double fixedDelta,
                                     double reMin, double reMax, int reSteps,
                                     double sMin, double sMax) {
        double delta = clamp(fixedDelta, DELTA_MIN, DELTA_MAX);
        List<ProfilePoint> profile = new ArrayList<>(reSteps);
        for (int i = 0; i < reSteps; i++) {
            double re = reSteps == 1 ? reMin : reMin + (reMax - reMin) * i / (double) (reSteps - 1);
            profile.add(profileAt(root, dates, delta, re, sMin, sMax));
        }
        return profile;
    }

    /**
     * Best log-likelihood attainable at one Re, maximizing over the sampling proportion.
     * <p>
     * s is found with Brent rather than a dense scan: a 100-point scan at every Re point made a
     * single estimate cost roughly 24,000 likelihood evaluations and about 80 seconds, which is not
     * a usable runtime for an index the pipeline computes on every dataset.
     */
    private ProfilePoint profileAt(PhyloNode root, Map<PhyloNode, Double> dates, double delta,
                                   double re, double sMin, double sMax) {
        double lo = clamp(sMin, S_MIN, S_MAX), hi = clamp(sMax, S_MIN, S_MAX);
        double bestS = Double.NaN, bestLogL = Double.NEGATIVE_INFINITY;
        // A handful of seeds first: the surface in s is not always unimodal, and Brent from a
        // single start can settle into the wrong basin.
        for (double seed : new double[]{lo, lo + (hi - lo) * 0.25, lo + (hi - lo) * 0.5,
                                        lo + (hi - lo) * 0.75, hi}) {
            double ll = logLikelihoodAt(root, dates, re, delta, seed);
            if (Double.isFinite(ll) && ll > bestLogL) { bestLogL = ll; bestS = seed; }
        }
        if (!Double.isFinite(bestLogL)) return new ProfilePoint(re, Double.NaN, Double.NEGATIVE_INFINITY);
        try {
            double refined = brentMaximize(sv -> logLikelihoodAt(root, dates, re, delta, sv), lo, hi, bestS);
            double ll = logLikelihoodAt(root, dates, re, delta, refined);
            if (Double.isFinite(ll) && ll > bestLogL) { bestLogL = ll; bestS = refined; }
        } catch (RuntimeException e) {
            // Brent could not bracket a maximum in s at this Re; the seeded best stands.
        }
        return new ProfilePoint(re, bestS, bestLogL);
    }

    /**
     * Two-stage profile: a coarse sweep to find the peak, then a fine sweep across the region the
     * likelihood-ratio interval can occupy. Same resolution as a single dense grid where it matters,
     * at a fraction of the evaluations.
     */
    List<ProfilePoint> refinedProfileOverRe(PhyloNode root, Map<PhyloNode, Double> dates, double fixedDelta,
                                            double reMin, double reMax, int coarseSteps, int fineSteps,
                                            double sMin, double sMax, double logLikelihoodDrop) {
        List<ProfilePoint> coarse = profileOverRe(root, dates, fixedDelta, reMin, reMax, coarseSteps, sMin, sMax);
        ProfilePoint peak = null;
        for (ProfilePoint p : coarse) if (peak == null || p.logLikelihood() > peak.logLikelihood()) peak = p;
        if (peak == null || !Double.isFinite(peak.logLikelihood())) return coarse;

        // Widen to whatever the coarse grid says is inside the interval, then one step either side
        // so the fine sweep can find the crossing rather than starting on top of it.
        double coarseStep = (reMax - reMin) / (coarseSteps - 1);
        double lo = reMax, hi = reMin;
        for (ProfilePoint p : coarse) {
            if (p.logLikelihood() >= peak.logLikelihood() - logLikelihoodDrop) {
                lo = Math.min(lo, p.re());
                hi = Math.max(hi, p.re());
            }
        }
        lo = Math.max(reMin, lo - coarseStep);
        hi = Math.min(reMax, hi + coarseStep);
        if (hi <= lo) return coarse;

        List<ProfilePoint> fine = profileOverRe(root, dates, fixedDelta, lo, hi, fineSteps, sMin, sMax);
        // Keep the coarse points outside the refined window so a caller can still see that the
        // likelihood falls away on both sides -- that is what makes an edge-running interval visible.
        List<ProfilePoint> merged = new ArrayList<>(coarse.size() + fine.size());
        for (ProfilePoint p : coarse) if (p.re() < lo || p.re() > hi) merged.add(p);
        merged.addAll(fine);
        merged.sort(java.util.Comparator.comparingDouble(ProfilePoint::re));
        return merged;
    }

    private double logLikelihoodAt(PhyloNode root, Map<PhyloNode, Double> dates, double re, double delta, double s) {
        double lambda = re * delta;
        double psi = s * delta;
        double mu = (1 - s) * delta;
        return new BdskyTreeLikelihood(lambda, mu, psi).logLikelihoodConditionedOnSurvival(root, dates);
    }

    private double brentMaximize(UnivariateFunction objective, double lower, double upper, double startingPoint) {
        BrentOptimizer optimizer = new BrentOptimizer(1e-6, 1e-10);
        double clampedStart = clamp(startingPoint, lower, upper);
        UnivariatePointValuePair result = optimizer.optimize(
                new MaxEval(MAX_EVALUATIONS_PER_PARAMETER),
                new UnivariateObjectiveFunction(objective),
                GoalType.MAXIMIZE,
                new SearchInterval(lower, upper, clampedStart));
        return result.getPoint();
    }

    private double clamp(double value, double lower, double upper) {
        return Math.max(lower, Math.min(upper, value));
    }
}
