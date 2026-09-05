package org.gvi.algorithms.re.bdsky;

import org.gvi.algorithms.phylo.PhyloNode;

import java.util.Map;

/**
 * The birth-death-sampling tree likelihood (Stadler 2010) -- the same
 * generative model BEAST2's BDSKY package maximizes via MCMC -- evaluated
 * via a post-order recursion over {@link BirthDeathSamplingPropagator}'s
 * p(tau)/logG(tau) propagators, for a FIXED, already-time-calibrated tree
 * (dates from {@code LeastSquaresDatingEstimator}), fitting the
 * constant-rate case via maximum likelihood rather than full Bayesian MCMC.
 * <p>
 * <b>Two real bugs were found and fixed here via independent Monte Carlo
 * cross-checks against a from-scratch Gillespie simulator</b> (see the
 * project README's Validation section and {@code GPropagatorMonteCarloTest}
 * / {@code Q1VerificationTest} for the derivations) -- documented in detail
 * because both are easy to reintroduce:
 * <p>
 * <b>1. p(tau)/logG(tau) must be evaluated on a GLOBAL time axis, not
 * per-edge-local elapsed time.</b> p(tau) means "probability of no sampled
 * descendants BY THE PRESENT, for a lineage tau time before it" -- tau is
 * remaining time <i>to the present</i>, a single shared reference point for
 * the whole tree, not "time since this specific edge started". An earlier
 * version called {@code propagator.integrate(edgeLength)} fresh for every
 * edge, implicitly treating each edge's own start as if it were "the
 * present" (p(0)=1) -- correct only for edges ending at the most recent
 * tip, wrong for every other edge. The fix: compute
 * {@code R(node) = presentTime - date(node)} (present = the latest date in
 * the tree) for every node, and an edge's logG contribution is
 * {@code logG(R(parent)) - logG(R(child))} -- both evaluated via fresh
 * {@code propagator.integrate(R)} calls from the SAME global tau=0 (since
 * p/logG are autonomous functions of elapsed duration, G(R_parent)/G(R_child)
 * -- the correct edge factor -- equals exactly this log-difference).
 * <p>
 * <b>2. The reference/root cannot be modeled as "sampled, with its lineage
 * continuing down to the first real branch point".</b> This project's
 * simulator (and the standard BDSKY convention) uses sampling WITH REMOVAL:
 * once a lineage is sampled, it cannot have further descendants. But this
 * codebase's {@code PhyloTree} is rooted AT the reference taxon (an
 * "outgroup-style" rooting -- see {@code PhyloTree.rootAt}), giving it
 * exactly one child (M) in the data structure -- treating the
 * reference-to-M edge as "the reference's own hidden continuation" (as an
 * earlier version did) contradicts sampling-with-removal: a lineage can't
 * be both sampled AND have descendants. The fix used here: treat M (the
 * root's one child -- always a genuine two-child branch point, or the
 * whole tree is just 2 taxa) as the effective origin for the likelihood,
 * and do not evaluate the reference's own ψ/edge contribution at all. This
 * is a deliberate, documented simplification (it doesn't use 100% of the
 * data -- the reference's own sampling event is dropped) rather than
 * introducing a new estimated "virtual origin date" parameter, which would
 * also require changes to the dating step and was judged higher-risk to
 * get right than the accuracy given up.
 */
final class BdskyTreeLikelihood {

    private final double lambda;
    private final double mu;
    private final double psi;
    private final BirthDeathSamplingPropagator propagator;

    BdskyTreeLikelihood(double lambda, double mu, double psi) {
        this.lambda = lambda;
        this.mu = mu;
        this.psi = psi;
        this.propagator = new BirthDeathSamplingPropagator(lambda, mu, psi);
    }

    /**
     * Divides the raw tree likelihood by {@code (1 - p(R(origin)))} -- the
     * probability the process, starting at the effective origin, produces
     * at least one sampled descendant by the present (BEAST2 BDSKY's
     * {@code conditionOnSurvival}). Verified empirically to make only a
     * small correction in typical regimes (see README) -- included for
     * correctness, not because it was the dominant source of the bugs
     * described in this class's javadoc.
     */
    double logLikelihoodConditionedOnSurvival(PhyloNode root, Map<PhyloNode, Double> dates) {
        PhyloNode origin = effectiveOrigin(root);
        double raw = subtreeLogLikelihood(origin, dates);
        double presentTime = dates.values().stream().mapToDouble(Double::doubleValue).max().orElseThrow();
        double originR = presentTime - dates.get(origin);
        double pNoSample = propagator.integrate(originR).p();
        double survivalProbability = 1.0 - pNoSample;
        return raw - Math.log(survivalProbability);
    }

    double logLikelihood(PhyloNode root, Map<PhyloNode, Double> dates) {
        return subtreeLogLikelihood(effectiveOrigin(root), dates);
    }

    /** See bug #2 in the class javadoc: start from the root's one child (a genuine 2-child branch point) rather than the reference tip itself. */
    private PhyloNode effectiveOrigin(PhyloNode root) {
        if (root.isTaxon() && root.children().size() == 1) {
            return root.children().get(0);
        }
        return root;
    }

    private double subtreeLogLikelihood(PhyloNode node, Map<PhyloNode, Double> dates) {
        double presentTime = dates.values().stream().mapToDouble(Double::doubleValue).max().orElseThrow();
        return subtreeLogLikelihood(node, dates, presentTime);
    }

    private double subtreeLogLikelihood(PhyloNode node, Map<PhyloNode, Double> dates, double presentTime) {
        double total = node.isTaxon() ? Math.log(psi) : Math.log(lambda);
        double rParent = presentTime - dates.get(node);
        for (PhyloNode child : node.children()) {
            double rChild = presentTime - dates.get(child);
            double logGParent = propagator.integrate(rParent).logG();
            double logGChild = propagator.integrate(rChild).logG();
            total += (logGParent - logGChild) + subtreeLogLikelihood(child, dates, presentTime);
        }
        return total;
    }

    double lambda() {
        return lambda;
    }

    double mu() {
        return mu;
    }

    double psi() {
        return psi;
    }
}
