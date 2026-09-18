package org.gvi.algorithms.mu;

import org.apache.commons.math3.stat.regression.SimpleRegression;
import org.gvi.algorithms.gd.GdMethod;
import org.gvi.algorithms.gd.GdResult;
import org.gvi.algorithms.gd.GeneticDistanceCalculator;
import org.gvi.algorithms.phylo.BootstrapSupportCalculator;
import org.gvi.algorithms.phylo.PhyloNode;
import org.gvi.algorithms.phylo.PhyloTree;
import org.gvi.algorithms.phylo.PhyloTreeFactory;
import org.gvi.algorithms.phylo.model.GtrParameterEstimator;
import org.gvi.algorithms.phylo.model.MlPhylogeneticOptimizer;
import org.gvi.algorithms.phylo.model.NamedSubstitutionModels;
import org.gvi.algorithms.phylo.model.NniTopologySearch;
import org.gvi.algorithms.phylo.model.RateParameterization;
import org.gvi.algorithms.phylo.model.SubstitutionModelSelector;
import org.gvi.core.exception.GviComputationException;
import org.gvi.core.model.NucleotideSequence;
import org.gvi.core.model.SequenceAlignment;
import org.gvi.core.util.TemporalUtil;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;

/**
 * Index 1 - Evolutionary Rate (mu), Section 5.1. Two estimators:
 * <p>
 * {@link #compute} -- baseline: distance from each dated sequence straight
 * to the reference, regressed against sampling date. Simple and fast, but
 * treats every sequence as independently diverged from the reference, which
 * ignores that the sequences share ancestry with each other too.
 * <p>
 * {@link #computeTreeAware} -- builds an actual Neighbor-Joining tree
 * (Saitou &amp; Nei 1987, via {@link PhyloTreeFactory}) rooted at the
 * reference, and regresses each dated taxon's cumulative root-to-tip branch
 * length (i.e. distance along the real inferred lineage, not a straight
 * line to the reference) against its sampling date. This is the standard
 * TempEst/root-to-tip approach as actually practiced (it operates on a
 * tree, not raw pairwise distances) and is more accurate whenever there's
 * meaningful shared branching structure among the samples. O(n^3) from the
 * underlying NJ construction, so it's only attempted up to
 * {@link PhyloTreeFactory#MAX_TAXA_FOR_TREE} taxa -- callers should catch
 * the exception and fall back to {@link #compute} for larger datasets.
 */
public final class EvolutionaryRateCalculator {

    private final GeneticDistanceCalculator distanceCalculator = new GeneticDistanceCalculator();

    /**
     * Baseline (non-tree) estimator's equivalent of {@link PhyloTreeFactory#temporalAnchorId}:
     * distance-to-reference has the exact same rooting-sensitivity problem as root-to-tip
     * distance on a tree, since "reference" plays the same role of a fixed comparison point.
     * Returns the earliest-dated sequence in the alignment, or the declared pipeline reference
     * when no sequence is dated or it already is the earliest.
     */
    private NucleotideSequence temporalAnchor(SequenceAlignment alignment, List<String> diagnostics) {
        NucleotideSequence declaredReference = alignment.getReference();
        NucleotideSequence earliest = null;
        double earliestYear = Double.POSITIVE_INFINITY;
        for (NucleotideSequence seq : alignment.getSequences()) {
            if (seq.getCollectionDate().isEmpty()) continue;
            double year = TemporalUtil.toDecimalYear(seq.getCollectionDate().get());
            if (year < earliestYear) {
                earliestYear = year;
                earliest = seq;
            }
        }
        if (earliest == null || earliest.getId().equals(declaredReference.getId())) {
            return declaredReference;
        }
        diagnostics.add(String.format(
                "Using '%s' (earliest collection date, %.2f) as the temporal anchor for this mu regression instead of "
                        + "the pipeline reference '%s' -- distance-to-reference from an arbitrary or recent --reference-id "
                        + "can otherwise invert or discard a real molecular-clock signal.",
                earliest.getId(), earliestYear, declaredReference.getId()));
        return earliest;
    }

    public MuResult compute(SequenceAlignment alignment) {
        return compute(alignment, GenomeType.UNSPECIFIED);
    }

    public MuResult compute(SequenceAlignment alignment, GenomeType genomeType) {
        List<String> diagnostics = new ArrayList<>();
        NucleotideSequence reference = temporalAnchor(alignment, diagnostics);
        List<double[]> points = new ArrayList<>();
        int saturatedExcluded = 0;
        int noDateExcluded = 0;

        for (NucleotideSequence seq : alignment.getSequences()) {
            if (seq.getCollectionDate().isEmpty()) {
                noDateExcluded++;
                continue;
            }
            double year = TemporalUtil.toDecimalYear(seq.getCollectionDate().get());
            GdResult gd = distanceCalculator.compute(reference, seq, GdMethod.JUKES_CANTOR);
            if (gd.saturated() || gd.distance() == null) {
                saturatedExcluded++;
                continue;
            }
            points.add(new double[]{year, gd.distance()});
        }

        if (noDateExcluded > 0) {
            diagnostics.add(noDateExcluded + " sequence(s) excluded: no collection_date in metadata");
        }
        if (saturatedExcluded > 0) {
            diagnostics.add(saturatedExcluded + " sequence(s) excluded: Jukes-Cantor distance to reference saturated");
        }

        return regressAndBuildResult(points, MuEstimationMethod.PAIRWISE_TO_REFERENCE, diagnostics, genomeType);
    }

    public MuResult computeTreeAware(SequenceAlignment alignment) {
        return computeTreeAware(alignment, false, BootstrapSupportCalculator.DEFAULT_REPLICATES, GenomeType.UNSPECIFIED);
    }

    /**
     * Same tree-aware estimator, optionally also running a Felsenstein
     * (1985) nonparametric bootstrap ({@link BootstrapSupportCalculator})
     * on the underlying NJ topology and appending a support summary to the
     * diagnostics -- the same "bootstrap support %" figure RAxML/IQ-TREE
     * report, absent from this project's trees until now. Opt-in (not the
     * default {@link #computeTreeAware(SequenceAlignment)} path) because it
     * multiplies tree-building cost by {@code bootstrapReplicates}; skipped
     * with an explanatory diagnostic above
     * {@link BootstrapSupportCalculator#MAX_TAXA_FOR_BOOTSTRAP} taxa.
     */
    public MuResult computeTreeAware(SequenceAlignment alignment, boolean withBootstrapSupport, int bootstrapReplicates) {
        return computeTreeAware(alignment, withBootstrapSupport, bootstrapReplicates, GenomeType.UNSPECIFIED);
    }

    public MuResult computeTreeAware(SequenceAlignment alignment, boolean withBootstrapSupport, int bootstrapReplicates, GenomeType genomeType) {
        List<String> rootingDiagnostics = new ArrayList<>();
        String rootId = PhyloTreeFactory.temporalAnchorId(alignment, rootingDiagnostics);
        PhyloTreeFactory.BuiltTree built = PhyloTreeFactory.build(alignment, GdMethod.JUKES_CANTOR, rootId);
        List<String> diagnostics = new ArrayList<>(rootingDiagnostics);
        diagnostics.addAll(built.diagnostics());
        if (withBootstrapSupport) {
            appendBootstrapSupportDiagnostic(alignment, built.tree(), bootstrapReplicates, diagnostics);
        }
        return computeTreeAware(alignment, built.tree(), diagnostics, MuEstimationMethod.TREE_ROOT_TO_TIP, genomeType);
    }

    /**
     * Native least-squares divergence-time dating -- the same estimation
     * problem LSD2 (To, Jung, Ly-Trong, Minh &amp; von Haeseler 2016) solves,
     * behind IQ-TREE's {@code --date} option: jointly fit a single
     * strict-clock rate AND every internal node's date directly against
     * every edge's own branch length (not one aggregate root-to-tip
     * regression line), enforcing that no node is dated later than its own
     * children. See {@link LeastSquaresDatingEstimator} for the derivation
     * and closed-form per-step update rules. Needs every sequence in the
     * alignment to have a collection date (unlike {@link #computeTreeAware},
     * which tolerates some undated sequences by simply excluding them from
     * the regression).
     */
    public MuResult computeLeastSquaresDating(SequenceAlignment alignment) {
        return computeLeastSquaresDating(alignment, GenomeType.UNSPECIFIED);
    }

    public MuResult computeLeastSquaresDating(SequenceAlignment alignment, GenomeType genomeType) {
        List<String> diagnostics = new ArrayList<>();
        String rootId = PhyloTreeFactory.temporalAnchorId(alignment, diagnostics);
        PhyloTreeFactory.BuiltTree built = PhyloTreeFactory.build(alignment, GdMethod.JUKES_CANTOR, rootId);
        diagnostics.addAll(built.diagnostics());

        Map<String, Double> tipDates = new HashMap<>();
        int noDateExcluded = 0;
        for (NucleotideSequence seq : alignment.getSequences()) {
            if (seq.getCollectionDate().isEmpty()) {
                noDateExcluded++;
                continue;
            }
            tipDates.put(seq.getId(), TemporalUtil.toDecimalYear(seq.getCollectionDate().get()));
        }
        if (noDateExcluded > 0) {
            throw new GviComputationException("Least-squares dating needs every sequence to have a collection date ("
                    + noDateExcluded + " missing); exclude undated sequences first or use the tree-aware estimator instead");
        }

        LeastSquaresDatingEstimator.Estimate estimate = new LeastSquaresDatingEstimator().estimate(built.tree(), tipDates);

        double meanBranchLength = estimate.nodeDates().keySet().stream()
                .filter(n -> n.parent() != null).mapToDouble(PhyloNode::branchLength).average().orElse(0.0);
        double ssTot = estimate.nodeDates().keySet().stream().filter(n -> n.parent() != null)
                .mapToDouble(n -> Math.pow(n.branchLength() - meanBranchLength, 2)).sum();
        double rSquared = ssTot > 0 ? Math.max(0.0, 1.0 - estimate.weightedSse() / ssTot) : 0.0;

        diagnostics.add(String.format("Least-squares dating (To et al. 2016 LSD2-equivalent objective): jointly fit mu and every "
                        + "internal node's date via constrained coordinate ascent (%d iteration(s), weighted SSE=%.6g) against every "
                        + "edge directly, instead of a single root-to-tip regression line.",
                estimate.iterationsUsed(), estimate.weightedSse()));
        if (rSquared < 0.3) {
            diagnostics.add(String.format("Weak fit (R^2=%.3f < 0.3): this rate estimate is unreliable", rSquared));
        }

        double mu = estimate.muSubPerSiteYear();
        String category;
        if (mu <= 0) {
            category = "N/A (non-positive rate estimate -- check temporal signal / possible sampling artifact)";
        } else {
            var band = MuReferenceTable.classify(mu, genomeType);
            category = MuReferenceTable.describe(band, genomeType);
        }

        double timeSpanYears = tipDates.values().stream().mapToDouble(Double::doubleValue).max().orElse(0.0)
                - tipDates.values().stream().mapToDouble(Double::doubleValue).min().orElse(0.0);
        return new MuResult(mu, rSquared, tipDates.size(), timeSpanYears, MuEstimationMethod.LEAST_SQUARES_DATING, category, diagnostics);
    }

    /**
     * Uncorrelated lognormal relaxed clock -- see {@link RelaxedClockMlEstimator} for the model and
     * fitting method. Like {@link #computeLeastSquaresDating}, needs every sequence dated (a relaxed
     * clock's per-branch rates are even less identifiable than a strict clock's single rate when
     * dates are missing, not more).
     */
    public MuResult computeRelaxedClock(SequenceAlignment alignment) {
        return computeRelaxedClock(alignment, GenomeType.UNSPECIFIED);
    }

    public MuResult computeRelaxedClock(SequenceAlignment alignment, GenomeType genomeType) {
        List<String> diagnostics = new ArrayList<>();
        String rootId = PhyloTreeFactory.temporalAnchorId(alignment, diagnostics);
        PhyloTreeFactory.BuiltTree built = PhyloTreeFactory.build(alignment, GdMethod.JUKES_CANTOR, rootId);
        diagnostics.addAll(built.diagnostics());

        Map<String, Double> tipDates = new HashMap<>();
        int noDateExcluded = 0;
        for (NucleotideSequence seq : alignment.getSequences()) {
            if (seq.getCollectionDate().isEmpty()) {
                noDateExcluded++;
                continue;
            }
            tipDates.put(seq.getId(), TemporalUtil.toDecimalYear(seq.getCollectionDate().get()));
        }
        if (noDateExcluded > 0) {
            throw new GviComputationException("Relaxed-clock dating needs every sequence to have a collection date ("
                    + noDateExcluded + " missing); exclude undated sequences first or use the tree-aware estimator instead");
        }

        RelaxedClockMlEstimator.Estimate estimate = new RelaxedClockMlEstimator().estimate(built.tree(), tipDates);

        double meanBranchLength = estimate.nodeDates().keySet().stream()
                .filter(n -> n.parent() != null).mapToDouble(PhyloNode::branchLength).average().orElse(0.0);
        double residualSse = 0.0, totalSse = 0.0;
        for (PhyloNode n : estimate.nodeDates().keySet()) {
            if (n.parent() == null) continue;
            double dt = estimate.nodeDates().get(n) - estimate.nodeDates().get(n.parent());
            double predicted = estimate.rate0SubPerSiteYear() * dt;
            residualSse += Math.pow(n.branchLength() - predicted, 2);
            totalSse += Math.pow(n.branchLength() - meanBranchLength, 2);
        }
        double rSquared = totalSse > 0 ? Math.max(0.0, 1.0 - residualSse / totalSse) : 0.0;

        diagnostics.add(String.format(java.util.Locale.ROOT,
                "Uncorrelated lognormal relaxed clock (Drummond et al. 2006 UCLD-equivalent model, fit by maximum "
                        + "likelihood): %d iteration(s), log-likelihood=%.6g, rate coefficient of variation=%.3f "
                        + "(0 = behaviourally a strict clock; the reported rate is the fitted clock's median branch rate).",
                estimate.iterationsUsed(), estimate.logLikelihood(), estimate.coefficientOfVariation()));

        double mu = estimate.rate0SubPerSiteYear();
        List<DateRandomizationTest.TemporalPoint> drtPoints = new ArrayList<>();
        for (NucleotideSequence seq : alignment.getSequences()) {
            Double year = tipDates.get(seq.getId());
            if (year != null) drtPoints.add(new DateRandomizationTest.TemporalPoint(year, built.tree().rootToTip(seq.getId())));
        }
        DateRandomizationTest.Result drt = new DateRandomizationTest().run(drtPoints, mu);
        diagnostics.add(drt.explanation());

        String category;
        if (mu <= 0) {
            category = "N/A (non-positive rate estimate -- check temporal signal / possible sampling artifact)";
        } else {
            var band = MuReferenceTable.classify(mu, genomeType);
            category = MuReferenceTable.describe(band, genomeType);
        }

        double timeSpanYears = tipDates.values().stream().mapToDouble(Double::doubleValue).max().orElse(0.0)
                - tipDates.values().stream().mapToDouble(Double::doubleValue).min().orElse(0.0);
        return new MuResult(mu, rSquared, tipDates.size(), timeSpanYears, MuEstimationMethod.RELAXED_CLOCK_ML, category,
                diagnostics, drt.passed(), estimate.coefficientOfVariation());
    }

    private void appendBootstrapSupportDiagnostic(SequenceAlignment alignment, PhyloTree tree, int replicates, List<String> diagnostics) {
        if (alignment.size() > BootstrapSupportCalculator.MAX_TAXA_FOR_BOOTSTRAP) {
            diagnostics.add("Bootstrap support skipped: " + alignment.size() + " taxa exceeds the "
                    + BootstrapSupportCalculator.MAX_TAXA_FOR_BOOTSTRAP
                    + "-taxon cap (cost multiplies tree-building by the replicate count)");
            return;
        }
        BootstrapSupportCalculator.Support support = new BootstrapSupportCalculator(replicates, 42L)
                .compute(alignment, GdMethod.JUKES_CANTOR, tree);
        diagnostics.addAll(support.diagnostics());
        if (support.supportByClade().isEmpty()) {
            diagnostics.add("Bootstrap support: no informative internal splits to test (tree too small)");
            return;
        }
        var values = support.supportByClade().values();
        double min = values.stream().mapToDouble(Double::doubleValue).min().orElse(0);
        double max = values.stream().mapToDouble(Double::doubleValue).max().orElse(0);
        double mean = values.stream().mapToDouble(Double::doubleValue).average().orElse(0);
        long weak = values.stream().filter(v -> v < 70.0).count();
        diagnostics.add(String.format(
                "Bootstrap support (Felsenstein 1985, %d replicate(s), %d internal split(s)): min=%.0f%%, mean=%.0f%%, max=%.0f%%%s",
                support.successfulReplicates(), values.size(), min, mean, max,
                weak > 0 ? String.format("; %d split(s) below the 70%% Hillis & Bull (1993) well-supported threshold -- treat that part of the topology cautiously", weak) : ""));
    }

    /** Reuses an already-built tree (e.g. shared with the LTT Re estimator) instead of constructing a new one. */
    public MuResult computeTreeAware(SequenceAlignment alignment, PhyloTree tree, List<String> diagnostics) {
        return computeTreeAware(alignment, tree, diagnostics, MuEstimationMethod.TREE_ROOT_TO_TIP, GenomeType.UNSPECIFIED);
    }

    private MuResult computeTreeAware(SequenceAlignment alignment, PhyloTree tree, List<String> diagnostics, MuEstimationMethod method, GenomeType genomeType) {
        List<double[]> points = new ArrayList<>();
        int noDateExcluded = 0;
        for (NucleotideSequence seq : alignment.getSequences()) {
            if (seq.getCollectionDate().isEmpty()) {
                noDateExcluded++;
                continue;
            }
            double year = TemporalUtil.toDecimalYear(seq.getCollectionDate().get());
            points.add(new double[]{year, tree.rootToTip(seq.getId())});
        }
        if (noDateExcluded > 0) {
            diagnostics.add(noDateExcluded + " sequence(s) excluded: no collection_date in metadata");
        }
        return regressAndBuildResult(points, method, diagnostics, genomeType);
    }

    /** Above this taxon count or alignment length, the ML branch-length optimization is not attempted automatically -- it is O(cycles x branches x Brent-evals x sites) and gets slow fast; this is an opt-in, not a default. */
    public static final int MAX_TAXA_FOR_ML = 40;
    public static final int MAX_SITES_FOR_ML = 50_000;
    private static final int DEFAULT_GAMMA_CATEGORIES = 4;
    // Coordinate ascent only guarantees a *local* optimum; a few randomized restarts materially reduce (not
    // eliminate) the chance a single unlucky starting point is mistaken for the true ML fit. Kept modest since
    // each restart re-runs the full coordinate-ascent cycle from scratch.
    private static final int ML_RESTART_COUNT = 2;
    private static final long ML_RESTART_SEED = 42L;

    /**
     * Highest-accuracy (and by far the slowest) mu estimator, using GTR
     * with empirical rate parameters and a fixed/no Gamma shape. Kept for
     * backward compatibility; prefer {@link #computeGtrGammaAware(SequenceAlignment, String, Double)}
     * which genuinely ML-fits the substitution rate parameters (and alpha,
     * if Gamma is enabled) instead of using empirical/fixed values, and
     * can automatically select the best-fitting model instead of assuming
     * GTR is the right one for this data.
     */
    public MuResult computeGtrGammaAware(SequenceAlignment alignment, Double gammaAlpha) {
        return computeGtrGammaAware(alignment, "gtr", gammaAlpha);
    }

    /**
     * Highest-accuracy (and by far the slowest) mu estimator: builds the NJ
     * topology, then jointly maximum-likelihood-optimizes every branch
     * length AND the substitution model's own rate parameters (and, if
     * Gamma is enabled, its alpha shape parameter) via coordinate ascent
     * ({@link MlPhylogeneticOptimizer}) -- a real accuracy upgrade over both
     * NJ's plain distance-based branch lengths and using empirical/assumed
     * substitution parameters. This is opt-in, not run automatically by
     * {@link #computeTreeAware}, because it is dramatically more expensive
     * (many Felsenstein-pruning likelihood evaluations per parameter per
     * refinement cycle -- and {@code modelName="auto"} multiplies that by
     * up to 6 candidate models). Explicitly capped at {@link #MAX_TAXA_FOR_ML}
     * taxa and {@link #MAX_SITES_FOR_ML} alignment sites; callers should
     * catch the exception and fall back to {@link #computeTreeAware} beyond that.
     *
     * @param modelName  one of jc69, f81, k80, hky85, tn93, gtr, or "auto" to fit all
     *                   6 and pick the best by AIC (see {@link SubstitutionModelSelector}) --
     *                   directly answers "what if this data doesn't actually fit a GTR model".
     * @param gammaAlpha if supplied, enables Gamma rate heterogeneity, ML-optimized starting
     *                   from this value; {@code null} skips Gamma entirely (single rate category).
     */
    public MuResult computeGtrGammaAware(SequenceAlignment alignment, String modelName, Double gammaAlpha) {
        return computeGtrGammaAware(alignment, modelName, gammaAlpha, GenomeType.UNSPECIFIED);
    }

    /** Same as {@link #computeGtrGammaAware(SequenceAlignment, String, Double)}, classifying the resulting mu against the given genome type instead of always assuming RNA. */
    public MuResult computeGtrGammaAware(SequenceAlignment alignment, String modelName, Double gammaAlpha, GenomeType genomeType) {
        if (alignment.size() > MAX_TAXA_FOR_ML) {
            throw new GviComputationException("High-accuracy ML mu skipped: " + alignment.size()
                    + " taxa exceeds the " + MAX_TAXA_FOR_ML + "-taxon cap for ML optimization; use computeTreeAware instead");
        }
        if (alignment.length() > MAX_SITES_FOR_ML) {
            throw new GviComputationException("High-accuracy ML mu skipped: alignment length " + alignment.length()
                    + " exceeds the " + MAX_SITES_FOR_ML + "-site cap for ML optimization; use computeTreeAware instead");
        }

        List<String> diagnostics = new ArrayList<>();
        int gammaCategories = gammaAlpha != null ? DEFAULT_GAMMA_CATEGORIES : 0;
        double initialAlpha = gammaAlpha != null ? gammaAlpha : 1.0;

        PhyloTree tree;
        Map<String, String> sequencesByLabel = new HashMap<>();
        for (NucleotideSequence seq : alignment.getSequences()) sequencesByLabel.put(seq.getId(), seq.getSequence());
        boolean nniEligible = alignment.size() <= NniTopologySearch.MAX_TAXA_FOR_NNI;

        if ("auto".equalsIgnoreCase(modelName)) {
            SubstitutionModelSelector.SelectionResult selection = new SubstitutionModelSelector()
                    .select(alignment, NamedSubstitutionModels.all(), gammaCategories,
                            PhyloTreeFactory.temporalAnchorId(alignment, diagnostics));
            tree = selection.treeUnderBestModel();
            StringBuilder comparison = new StringBuilder("Model selection (AIC, lower=better): ");
            for (var fit : selection.allFits()) {
                comparison.append(String.format("%s(AIC=%.1f) ", fit.modelName(), fit.aic()));
            }
            diagnostics.add(comparison.toString().trim());
            diagnostics.add(String.format("Selected %s: logL=%.2f, %d free parameter(s)%s",
                    selection.best().modelName(), selection.best().logLikelihood(), selection.best().freeParameters(),
                    selection.best().optimizedAlpha() != null ? String.format(", gamma alpha=%.3f", selection.best().optimizedAlpha()) : ""));

            if (nniEligible) {
                NniTopologySearch.Result nni = new NniTopologySearch().search(
                        tree, sequencesByLabel, alignment.length(), selection.bestModel(), selection.bestGammaRates());
                if (nni.movesAccepted() > 0) {
                    diagnostics.add(String.format("NNI topology search (auto): %d move(s) accepted over %d round(s), "
                                    + "log-likelihood %.2f -> %.2f -- refitting model parameters on the new topology",
                            nni.movesAccepted(), nni.roundsUsed(), nni.initialLogLikelihood(), nni.finalLogLikelihood()));
                    RateParameterization winningParameterization = withInitialTheta(
                            NamedSubstitutionModels.byName(selection.best().modelName()), selection.best().optimizedTheta());
                    GtrParameterEstimator.Estimate estimate = GtrParameterEstimator.estimate(alignment);
                    double refitAlpha = selection.best().optimizedAlpha() != null ? selection.best().optimizedAlpha() : initialAlpha;
                    MlPhylogeneticOptimizer.Result refit = new MlPhylogeneticOptimizer().optimize(tree, sequencesByLabel,
                            alignment.length(), winningParameterization, estimate.baseFrequencies(), gammaCategories, refitAlpha);
                    diagnostics.add(String.format("Post-NNI refit: log-likelihood -> %.2f over %d cycle(s)",
                            refit.finalLogLikelihood(), refit.cyclesUsed()));
                } else {
                    diagnostics.add("NNI topology search (auto): no rearrangement improved on the Neighbor-Joining topology");
                }
            } else {
                diagnostics.add("NNI topology search skipped: " + alignment.size() + " taxa exceeds the "
                        + NniTopologySearch.MAX_TAXA_FOR_NNI + "-taxon cap");
            }
        } else {
            RateParameterization parameterization = NamedSubstitutionModels.byName(modelName);
            GtrParameterEstimator.Estimate estimate = GtrParameterEstimator.estimate(alignment);
            diagnostics.addAll(estimate.diagnostics());

            PhyloTreeFactory.BuiltTree built = PhyloTreeFactory.build(alignment, GdMethod.JUKES_CANTOR,
                    PhyloTreeFactory.temporalAnchorId(alignment, diagnostics));
            diagnostics.addAll(built.diagnostics());
            tree = built.tree();

            MlPhylogeneticOptimizer.Result result = new MlPhylogeneticOptimizer().optimizeWithRestarts(tree, sequencesByLabel,
                    alignment.length(), parameterization, estimate.baseFrequencies(), gammaCategories, initialAlpha,
                    ML_RESTART_COUNT, ML_RESTART_SEED);
            diagnostics.add(String.format("%s ML fit (%d restart(s)): log-likelihood %.2f -> %.2f over %d cycle(s)%s%s",
                    parameterization.modelName(), ML_RESTART_COUNT, result.initialLogLikelihood(), result.finalLogLikelihood(), result.cyclesUsed(),
                    result.optimizedTheta().length > 0 ? "; rate params=" + java.util.Arrays.toString(result.optimizedTheta()) : "",
                    result.optimizedAlpha() != null ? String.format("; gamma alpha=%.3f", result.optimizedAlpha()) : ""));

            if (nniEligible) {
                NniTopologySearch.Result nni = new NniTopologySearch().search(
                        tree, sequencesByLabel, alignment.length(), result.finalModel(), result.finalGammaRates());
                if (nni.movesAccepted() > 0) {
                    diagnostics.add(String.format("NNI topology search: %d move(s) accepted over %d round(s), "
                                    + "log-likelihood %.2f -> %.2f -- refitting model parameters on the new topology",
                            nni.movesAccepted(), nni.roundsUsed(), nni.initialLogLikelihood(), nni.finalLogLikelihood()));
                    RateParameterization refitParameterization = withInitialTheta(parameterization, result.optimizedTheta());
                    double refitAlpha = result.optimizedAlpha() != null ? result.optimizedAlpha() : initialAlpha;
                    MlPhylogeneticOptimizer.Result refit = new MlPhylogeneticOptimizer().optimize(tree, sequencesByLabel,
                            alignment.length(), refitParameterization, estimate.baseFrequencies(), gammaCategories, refitAlpha);
                    diagnostics.add(String.format("Post-NNI refit: log-likelihood -> %.2f over %d cycle(s)",
                            refit.finalLogLikelihood(), refit.cyclesUsed()));
                } else {
                    diagnostics.add("NNI topology search: no rearrangement improved on the Neighbor-Joining topology");
                }
            } else {
                diagnostics.add("NNI topology search skipped: " + alignment.size() + " taxa exceeds the "
                        + NniTopologySearch.MAX_TAXA_FOR_NNI + "-taxon cap");
            }
        }

        return computeTreeAware(alignment, tree, diagnostics, MuEstimationMethod.GTR_GAMMA_ML_BRANCH_LENGTHS, genomeType);
    }

    private RateParameterization withInitialTheta(RateParameterization base, double[] initialTheta) {
        return new RateParameterization(base.modelName(), base.useEmpiricalFrequencies(), initialTheta,
                base.lowerBounds(), base.upperBounds(), base.thetaToRates());
    }

    private MuResult regressAndBuildResult(List<double[]> points, MuEstimationMethod method, List<String> diagnostics, GenomeType genomeType) {
        TreeSet<Double> distinctYears = new TreeSet<>();
        for (double[] p : points) distinctYears.add(p[0]);
        if (distinctYears.size() < 2) {
            throw new GviComputationException(
                    "Cannot estimate mu: need at least 2 distinct collection dates with usable distances, found "
                            + distinctYears.size() + ". Provide a temporally spread dataset.");
        }

        SimpleRegression regression = new SimpleRegression();
        for (double[] p : points) regression.addData(p[0], p[1]);

        double mu = regression.getSlope();
        double rSquared = Double.isNaN(regression.getRSquare()) ? 0.0 : regression.getRSquare();
        double timeSpan = distinctYears.last() - distinctYears.first();

        if (rSquared < 0.3) {
            diagnostics.add(String.format(
                    "Temporal signal is weak (R^2=%.3f < 0.3) -- treat this mu estimate as unreliable, not a precise rate", rSquared));
        } else {
            diagnostics.add(String.format("Temporal signal R^2=%.3f over %.2f years (n=%d)", rSquared, timeSpan, points.size()));
        }

        // Date randomization (Ramsden et al. 2009; Duchene et al. 2015). R^2 alone cannot tell a
        // real molecular clock from several independently-introduced lineages pooled into one
        // alignment: the latter produces a tight regression whose slope has nothing to do with
        // elapsed time. Shuffling the dates tests that directly.
        List<DateRandomizationTest.TemporalPoint> drtPoints = new ArrayList<>();
        for (double[] p : points) drtPoints.add(new DateRandomizationTest.TemporalPoint(p[0], p[1]));
        DateRandomizationTest.Result drt = new DateRandomizationTest().run(drtPoints, mu);
        diagnostics.add(drt.explanation());

        String category;
        if (mu <= 0) {
            category = "N/A (non-positive rate estimate -- check temporal signal / possible sampling artifact)";
        } else {
            var band = MuReferenceTable.classify(mu, genomeType);
            category = MuReferenceTable.describe(band, genomeType);
        }

        return new MuResult(mu, rSquared, points.size(), timeSpan, method, category, diagnostics, drt.passed());
    }
}
