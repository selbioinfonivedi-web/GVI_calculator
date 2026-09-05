package org.gvi.algorithms.re.bdsky;

import org.gvi.algorithms.gd.GdMethod;
import org.gvi.algorithms.mu.LeastSquaresDatingEstimator;
import org.gvi.algorithms.phylo.PhyloTreeFactory;
import org.gvi.algorithms.re.ReMethod;
import org.gvi.algorithms.re.ReReferenceTable;
import org.gvi.algorithms.re.ReResult;
import org.gvi.core.exception.GviComputationException;
import org.gvi.core.model.ConfidenceInterval;
import org.gvi.core.model.NucleotideSequence;
import org.gvi.core.model.SequenceAlignment;
import org.gvi.core.util.TemporalUtil;

import java.time.LocalDate;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;

/**
 * Public entry point for the native birth-death-sampling (BDSKY-equivalent)
 * Re estimator -- see this package's {@code package-info.java} for the full
 * validation history, including the two real bugs found and fixed via
 * independent Monte Carlo cross-checks before this was wired in. Builds an
 * NJ tree, dates every node via {@link LeastSquaresDatingEstimator} (needs
 * every sequence dated), then fits (Re, becomeUninfectiousRate,
 * samplingProportion) via {@link BdskyMlFitter}.
 * <p>
 * A real accuracy upgrade over the phylodynamic LTT/Euler-Lotka fallback
 * for datasets with no case-incidence data -- it uses the actual
 * birth-death-sampling tree likelihood (the same generative model BEAST2's
 * BDSKY maximizes) rather than a growth-rate summary statistic. Still not
 * BDSKY itself: a maximum-likelihood point estimate via coordinate ascent,
 * not full Bayesian MCMC with credible intervals, and
 * becomeUninfectiousRate/samplingProportion are not individually
 * well-identified without informative priors (a documented BDSKY
 * literature limitation, not specific to this implementation) -- Re itself
 * is the well-recovered quantity.
 */
public final class BdskyReEstimator {

    /** Likelihood evaluation is O(taxa) ODE integrations, evaluated many times per coordinate-ascent cycle -- capped tighter than the general 300-taxon NJ limit to keep a fit tractable. */
    public static final int MAX_TAXA_FOR_BDSKY = 60;

    /** Range and resolution of the Re profile grid. */
    private static final double RE_PROFILE_MIN = 0.10, RE_PROFILE_MAX = 12.0;
    private static final int RE_PROFILE_COARSE_STEPS = 60;
    private static final int RE_PROFILE_FINE_STEPS = 80;
    /**
     * Sampling proportion is profiled over a band rather than fitted freely. A value near 1 means
     * essentially every infection in the outbreak was sequenced, which no surveillance system
     * achieves; leaving the upper bound at 1 let the optimizer buy likelihood by driving s to 0.96
     * and moving Re to match.
     */
    private static final double SAMPLING_PROPORTION_MIN = 0.002, SAMPLING_PROPORTION_MAX = 0.80;
    /** chi-square(1) 95% quantile halved: the standard profile-likelihood interval cut-off. */
    private static final double PROFILE_CI_LOG_LIKELIHOOD_DROP = 1.92;

    public ReResult compute(SequenceAlignment alignment, double generationTimeDays) {
        if (alignment.size() > MAX_TAXA_FOR_BDSKY) {
            throw new GviComputationException("BDSKY-ML Re skipped: " + alignment.size() + " taxa exceeds the "
                    + MAX_TAXA_FOR_BDSKY + "-taxon cap; use the phylodynamic fallback instead");
        }

        long datedCount = alignment.getSequences().stream().filter(s -> s.getCollectionDate().isPresent()).count();
        if (datedCount < alignment.size()) {
            throw new GviComputationException("BDSKY-ML Re needs every sequence dated ("
                    + (alignment.size() - datedCount) + " missing); use the phylodynamic fallback instead");
        }

        List<String> diagnostics = new ArrayList<>();
        String rootId = PhyloTreeFactory.temporalAnchorId(alignment, diagnostics);
        PhyloTreeFactory.BuiltTree built = PhyloTreeFactory.build(alignment, GdMethod.JUKES_CANTOR, rootId);
        diagnostics.addAll(built.diagnostics());

        Map<String, Double> tipDates = new HashMap<>();
        for (NucleotideSequence seq : alignment.getSequences()) {
            tipDates.put(seq.getId(), TemporalUtil.toDecimalYear(seq.getCollectionDate().get()));
        }
        LeastSquaresDatingEstimator.Estimate dating = new LeastSquaresDatingEstimator().estimate(built.tree(), tipDates);

        // delta is the reciprocal of the infectious period, which the generation-time table supplies.
        // It is held here, not merely used as a starting guess: see BdskyMlFitter.profileOverRe.
        double delta = 365.25 / generationTimeDays;
        List<BdskyMlFitter.ProfilePoint> profile = new BdskyMlFitter().refinedProfileOverRe(
                built.tree().root(), dating.nodeDates(), delta,
                RE_PROFILE_MIN, RE_PROFILE_MAX, RE_PROFILE_COARSE_STEPS, RE_PROFILE_FINE_STEPS,
                SAMPLING_PROPORTION_MIN, SAMPLING_PROPORTION_MAX, PROFILE_CI_LOG_LIKELIHOOD_DROP);

        BdskyMlFitter.ProfilePoint peak = profile.stream()
                .max(Comparator.comparingDouble(BdskyMlFitter.ProfilePoint::logLikelihood))
                .orElseThrow(() -> new GviComputationException("birth-death profile produced no evaluable points"));
        if (!Double.isFinite(peak.logLikelihood())) {
            throw new GviComputationException("the birth-death likelihood is not finite anywhere on the Re profile for "
                    + "this tree -- the dated topology is degenerate (typically identical sequences or zero temporal spread)");
        }

        // Likelihood-ratio interval: the set of Re whose profile log-likelihood is within
        // 1.92 of the peak, i.e. chi-square(1) at 95% halved. Reported rather than a bare
        // point estimate because on realistic tree sizes it is wide, and the width IS the finding.
        double threshold = peak.logLikelihood() - PROFILE_CI_LOG_LIKELIHOOD_DROP;
        double lower = Double.NaN, upper = Double.NaN;
        for (BdskyMlFitter.ProfilePoint pt : profile) {
            if (pt.logLikelihood() >= threshold) {
                if (Double.isNaN(lower)) lower = pt.re();
                upper = pt.re();
            }
        }

        // An interval running to the edge of the profile grid means the likelihood had not turned
        // over yet: the data do not bound Re on that side, and a point estimate drawn from inside
        // such an interval is the grid's property, not the tree's.
        double gridStep = (RE_PROFILE_MAX - RE_PROFILE_MIN) / (RE_PROFILE_COARSE_STEPS - 1);
        if (lower <= RE_PROFILE_MIN + gridStep || upper >= RE_PROFILE_MAX - gridStep) {
            throw new GviComputationException(String.format(
                    "the birth-death profile likelihood does not bound Re from this tree: the 95%% likelihood-ratio "
                            + "interval [%.2f, %.2f] runs to the edge of the searched range [%.1f, %.1f], so the "
                            + "likelihood was still rising at the boundary -- typically too few taxa, too little "
                            + "temporal spread, or an alignment pooling unrelated lineages",
                    lower, upper, RE_PROFILE_MIN, RE_PROFILE_MAX));
        }

        diagnostics.add(String.format("Native birth-death-sampling profile-likelihood fit (BDSKY-equivalent generative model, "
                        + "maximum likelihood rather than Bayesian MCMC): Re=%.4f, 95%% likelihood-ratio interval [%.2f, %.2f]. "
                        + "becomeUninfectiousRate held at %.4f/yr, the reciprocal of the supplied %.1f-day generation time; "
                        + "samplingProportion profiled over [%.3f, %.2f] and best at %.3f; log-likelihood=%.4f over %d Re grid points. "
                        + "delta is held rather than fitted because (Re, delta, s) has a likelihood ridge along which a free "
                        + "three-parameter fit returns whatever its starting point happened to be near.",
                peak.re(), lower, upper, delta, generationTimeDays,
                SAMPLING_PROPORTION_MIN, SAMPLING_PROPORTION_MAX, peak.samplingProportion(),
                peak.logLikelihood(), profile.size()));
        diagnostics.add(String.format("Re is sensitive to the generation time: because delta = 365.25/generationTime scales "
                        + "the whole birth-death process, a generation time wrong by a factor of k moves Re by roughly the "
                        + "same factor. The %.1f-day value used here came from --generation-time-days or the bundled table; "
                        + "check it is right for this pathogen before quoting the number.", generationTimeDays));
        if ((upper - lower) > peak.re()) {
            diagnostics.add(String.format("Wide interval: the 95%% range [%.2f, %.2f] is broader than the point estimate "
                    + "itself, so this tree constrains Re only loosely. Treat the interval, not the point, as the result.",
                    lower, upper));
        }

        var band = ReReferenceTable.classify(peak.re());
        String category = band.category() + " (" + band.transmission() + "); " + band.controlMeasure();

        LocalDate lastDate = alignment.getSequences().stream()
                .map(NucleotideSequence::getCollectionDate).filter(Optional::isPresent).map(Optional::get)
                .max(LocalDate::compareTo).orElse(LocalDate.now());

        Double doublingTimeDays = null;
        double netGrowthRate = (peak.re() - 1.0) * delta; // Re=lambda/delta, so lambda-delta = (Re-1)*delta
        if (netGrowthRate > 0) {
            doublingTimeDays = (Math.log(2) / netGrowthRate) * 365.25;
        }

        return new ReResult(ReMethod.BDSKY_ML, lastDate, peak.re(),
                new ConfidenceInterval(lower, upper, 0.95), List.of(), doublingTimeDays, category, diagnostics);
    }
}
