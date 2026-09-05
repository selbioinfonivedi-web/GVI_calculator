/**
 * A native maximum-likelihood codon-substitution dN/dS estimator --
 * Goldman &amp; Yang (1994), GY94 -- the same core continuous-time Markov
 * model PAML's {@code codeml} builds its M0 "one-ratio" analysis on: a
 * single instantaneous-rate matrix over the 61 sense codons, parameterized
 * by kappa (transition/transversion ratio), omega (dN/dS itself), and F3x4
 * codon frequencies, with branch lengths + kappa + omega jointly
 * maximum-likelihood fit via coordinate ascent (Felsenstein pruning as the
 * likelihood, Brent's method per parameter -- the same established pattern
 * {@link org.gvi.algorithms.phylo.model.MlPhylogeneticOptimizer} uses for
 * nucleotide GTR(+Gamma)).
 * <p>
 * <b>What this gives you that {@link org.gvi.algorithms.dnds.DnDsCalculator}'s
 * default Nei-Gojobori counting method doesn't:</b> proper correction for
 * multiple substitutions at a site via an explicit CTMC (rather than a
 * post-hoc Jukes-Cantor distance correction), and a single ML omega jointly
 * estimated across the whole tree rather than pairwise/pooled counting.
 * <p>
 * <b>What this deliberately does NOT attempt</b> -- the part of
 * {@code codeml} that is genuinely a much larger undertaking and remains
 * out of scope: site-class mixture models (M1a/M2a nearly-neutral vs
 * positive-selection, M7/M8 Beta-distributed omega classes), branch-site
 * tests, likelihood-ratio tests between nested models, and Bayes Empirical
 * Bayes per-site posterior classification. Those are what let real
 * {@code codeml} runs detect selection acting on a handful of sites against
 * a genome-wide background of purifying selection; this model gives one
 * omega for the whole gene (codeml's M0 baseline), not per-site resolution
 * -- {@link org.gvi.algorithms.dnds.DnDsCalculator}'s sliding-window scan
 * remains this project's (non-ML) answer to site-level resolution.
 * <p>
 * <b>Validation:</b> {@code CodonModelTest} checks basic CTMC correctness
 * properties (P(t) row-stochasticity, identity at t-&gt;0, convergence to
 * the stationary distribution at t-&gt;infinity). {@code CodonMlFitterGroundTruthTest}
 * simulates codon evolution under KNOWN kappa/omega via an independent
 * Gillespie (Doob-Gillespie SSA) simulator that reimplements the GY94 rate
 * formula directly (never calling {@link org.gvi.algorithms.dnds.ml.CodonModel}
 * or the fitter itself) and checks the ML fit recovers values close to the
 * true generating parameters, for both a purifying-selection and a
 * positive-selection scenario.
 * <p>
 * <b>Cost / caps:</b> {@link org.gvi.algorithms.dnds.ml.MlCodonDnDsEstimator}
 * is capped at 15 taxa and 300 codons -- much tighter than the nucleotide
 * ML paths, because the per-site pruning cost scales with the SQUARE of the
 * state count (61 vs 4, ~232x more expensive per site), and every Brent
 * evaluation during branch-length optimization still recomputes the whole
 * tree's likelihood (the same computationally simple, established pattern
 * {@link org.gvi.algorithms.phylo.model.MlBranchLengthOptimizer} uses for
 * nucleotides -- not reimplemented as an incremental/partial update here).
 * Genes/datasets exceeding these caps fall back to Nei-Gojobori automatically
 * (see {@code GviPipeline}'s wiring of {@code --ml-dnds}).
 */
package org.gvi.algorithms.dnds.ml;
