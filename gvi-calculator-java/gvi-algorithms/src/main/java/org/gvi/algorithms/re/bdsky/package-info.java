/**
 * Native birth-death-sampling process machinery (the same generative model
 * BEAST2's BDSKY package uses for Re), fit via maximum likelihood rather
 * than full Bayesian MCMC -- STATUS: validated, wired into the pipeline
 * as the opt-in {@code --bdsky-re} tier.
 * <p>
 * Independently verified layer by layer, each against a from-scratch
 * source of truth (never against this project's own other outputs):
 * <ul>
 *   <li>{@link org.gvi.algorithms.re.bdsky.BirthDeathSamplingPropagator} --
 *       p(tau)/logG(tau), checked against two closed-form special cases
 *       AND (for generic, non-degenerate rates) direct Monte Carlo
 *       simulation ({@code MonteCarloVerificationTest}).</li>
 *   <li>The relationship between logG and the probability of an exact
 *       sample count -- {@code Q1VerificationTest} -- checked via an
 *       independently-derived inhomogeneous ODE against Monte Carlo.</li>
 *   <li>{@link org.gvi.algorithms.re.bdsky.BdskyTreeLikelihood} -- the
 *       tree-recursion composition, checked against a hand-computed value
 *       on a minimal tree.</li>
 *   <li>{@link org.gvi.algorithms.re.bdsky.BdskyMlFitter} -- end-to-end,
 *       via a full Gillespie birth-death-sampling simulation -> real JC69
 *       sequence evolution -> real NJ reconstruction -> real LSD dating ->
 *       ML fit, checked against the known true Re
 *       ({@code BdskyMlFitterGroundTruthTest}).</li>
 * </ul>
 * <p>
 * <b>Two real bugs were found and fixed during this validation</b> (both
 * documented in detail in {@link org.gvi.algorithms.re.bdsky.BdskyTreeLikelihood}'s
 * class javadoc, since they're easy to reintroduce):
 * <ol>
 *   <li>p(tau)/logG(tau) must be evaluated on a GLOBAL remaining-time-to-
 *       present axis, not per-edge-local elapsed time -- caught when a
 *       direct Monte Carlo check of psi*g(tau) (the density actually used
 *       at every tree event) came out ~2x off from simulation, even though
 *       p(tau) alone was independently confirmed correct.</li>
 *   <li>The reference/root (this codebase roots trees at a sampled taxon,
 *       not a true bifurcating MRCA) cannot be modeled as "sampled AND
 *       having continuing descendants" -- inconsistent with
 *       sampling-with-removal. Fixed by starting the likelihood at the
 *       root's one child (a genuine branch point) and dropping the
 *       reference's own contribution -- a documented, deliberate small
 *       loss of information rather than introducing a new estimated
 *       "virtual origin date" parameter.</li>
 * </ol>
 * Both were found via INDEPENDENT Monte Carlo cross-checks against a
 * from-scratch Gillespie simulator (not by loosening a tolerance until an
 * end-to-end test passed) -- see the README's Validation section for the
 * full diagnostic trail. After both fixes, the end-to-end ground-truth
 * test shows a genuine, well-behaved likelihood (a single interior maximum,
 * consistent maximum reached from multiple starting points) rather than
 * the runaway-to-implausible-values behavior that flagged the bugs in the
 * first place.
 * <p>
 * Known remaining limitation, consistent with the wider BDSKY literature:
 * become-uninfectious-rate (delta) and sampling-proportion (s) are not
 * individually well-identified from tree shape alone without informative
 * priors (real BEAST2 BDSKY uses a concentrated Beta prior on s for
 * exactly this reason) -- this ML fitter has no priors, so delta/s can
 * vary across different starting points while landing on essentially the
 * SAME maximum likelihood and a similar Re. Re itself -- the actual
 * quantity this tier is for -- is comparatively well recovered.
 */
package org.gvi.algorithms.re.bdsky;
