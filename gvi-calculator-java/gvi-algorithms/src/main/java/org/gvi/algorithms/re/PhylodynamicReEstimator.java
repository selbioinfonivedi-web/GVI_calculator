package org.gvi.algorithms.re;

import org.apache.commons.math3.stat.regression.SimpleRegression;
import org.gvi.algorithms.gd.GdMethod;
import org.gvi.algorithms.phylo.PhyloNode;
import org.gvi.algorithms.phylo.PhyloTree;
import org.gvi.algorithms.phylo.PhyloTreeFactory;
import org.gvi.core.exception.GviComputationException;
import org.gvi.core.model.NucleotideSequence;
import org.gvi.core.model.SequenceAlignment;
import org.gvi.core.util.TemporalUtil;

import java.time.LocalDate;
import java.util.ArrayList;
import java.util.List;
import java.util.TreeSet;

/**
 * Index 2 - Effective Reproduction Number, fallback estimator (Section
 * 5.2b) used when no case-incidence series is available, only dated
 * sequences. Implements the classic lineages-through-time (LTT) growth-rate
 * method (Pybus, Rambaut &amp; Harvey, 2000): build a genealogy, count how
 * many lineages are present at each point in time as the tree is traversed
 * from root to present, and regress ln(lineage count) against time -- the
 * slope is the exponential growth rate r.
 * <p>
 * Concretely:
 * <ol>
 *   <li>Build a Neighbor-Joining tree ({@link PhyloTreeFactory}) rooted at
 *       the reference.</li>
 *   <li>Fit a root-to-tip regression (distance vs. sampling date) on that
 *       same tree to get a local substitution rate mu and the implied
 *       calendar time of the root (a strict-clock time-scaling of the tree,
 *       the same core idea TempEst/root-to-tip dating uses to place dates on
 *       internal nodes, just without its more elaborate outlier handling).</li>
 *   <li>Use mu to convert every internal branch point's root-to-tip
 *       substitution distance into an estimated calendar time.</li>
 *   <li>Count lineages through time across those branch points and regress
 *       ln(count) against estimated time -- Pybus &amp; Rambaut's method.</li>
 *   <li>Convert the growth rate r to Re via the Euler-Lotka renewal
 *       equation (Wallinga &amp; Lipsitch 2007): {@code 1 = Re * sum_u w(u) *
 *       e^(-r*u)} over the discretized serial-interval distribution w(u) --
 *       the same {@link SerialInterval} machinery the Cori et al. incidence
 *       estimator uses, so both Re paths are grounded in the same
 *       assumption about how long an infectious interval actually is,
 *       instead of the fallback using a cruder ad-hoc conversion. This
 *       exactly generalizes the older {@code Re = 1 + r*generationTime}
 *       linear approximation (a first-order Taylor expansion of this same
 *       equation, verified to agree with it in the small-r limit and to
 *       diverge -- more accurately -- for faster growth/decline).</li>
 * </ol>
 * This is a real, citable, tractable method -- but simpler than full
 * coalescent-likelihood or birth-death skyline inference (e.g. BEAST's
 * BDSKY, Stadler et al. 2013), which jointly estimate the tree, its
 * time-calibration, and a full birth-death-sampling process likelihood by
 * MCMC (with credible intervals) instead of point-estimating each step in
 * sequence via a simpler, decades-older growth-rate method. Flagged in
 * diagnostics as an approximation, and only used when no case-incidence
 * data is available.
 */
public final class PhylodynamicReEstimator {

    /**
     * The true serial-interval standard deviation is pathogen-specific and
     * not knowable from a single generation-time scalar (the CLI/GUI only
     * expose a mean via {@code --generation-time-days}) -- this ratio
     * matches {@link SerialInterval#defaultProfile()}'s own mean=5/sd=2
     * shape (sd/mean = 0.4) so a custom generation time gets a
     * proportionally-shaped Gamma distribution rather than an arbitrary one.
     */
    private static final double SERIAL_INTERVAL_SD_TO_MEAN_RATIO = 0.4;

    /**
     * Lower bound on the discretization horizon, retained so a short (viral) generation time
     * produces exactly the same 20-day grid this estimator has always used.
     */
    private static final int MIN_SERIAL_INTERVAL_T_MAX_DAYS = 20;

    /**
     * How many standard deviations past the mean the serial-interval grid must reach.
     * <p>
     * This horizon used to be a hardcoded 20 days for every pathogen, which silently destroyed
     * the Euler-Lotka conversion for anything slower than a respiratory virus. {@link SerialInterval}
     * renormalizes its weights over [1, tMax], so a Gamma with a 457-day mean evaluated on a 20-day
     * grid does not become a truncated 457-day distribution -- it becomes a ~20-day one. The
     * conversion Re = 1 / sum_u w(u) e^(-r u) then reduces to roughly 1 + r*(20/365), which is
     * within a few percent of 1.0 for any realistic growth rate. Measured on simulated
     * birth-death-sampling trees with known Re, the estimator returned 1.03, 1.07 and 1.03 for
     * true Re of 1.5, 3.0 and 1.125 -- it could not separate a tripling epidemic from a static one.
     * <p>
     * Six standard deviations past the mean captures the Gamma's mass to better than 1e-4 at the
     * sd/mean ratio used here, so the discretized weights represent the distribution actually asked for.
     */
    private static final double SERIAL_INTERVAL_SD_HORIZON = 6.0;

    /** Discretization horizon for a serial interval of the given shape, in days. */
    static int serialIntervalTMax(double meanDays, double sdDays) {
        int needed = (int) Math.ceil(meanDays + SERIAL_INTERVAL_SD_HORIZON * sdDays);
        return Math.max(MIN_SERIAL_INTERVAL_T_MAX_DAYS, needed);
    }

    public ReResult compute(SequenceAlignment alignment, double generationTimeDays) {
        long datedCount = alignment.getSequences().stream().filter(s -> s.getCollectionDate().isPresent()).count();
        if (datedCount < 5) {
            throw new GviComputationException(
                    "Phylodynamic Re fallback needs at least 5 dated sequences, got " + datedCount);
        }

        List<String> diagnostics = new ArrayList<>();
        String rootId = PhyloTreeFactory.temporalAnchorId(alignment, diagnostics);
        PhyloTreeFactory.BuiltTree built = PhyloTreeFactory.build(alignment, GdMethod.JUKES_CANTOR, rootId);
        PhyloTree tree = built.tree();
        diagnostics.addAll(built.diagnostics());

        // Step 2: root-to-tip regression on this tree to time-scale it (same idea as tree-aware mu).
        SimpleRegression clockRegression = new SimpleRegression();
        TreeSet<Double> distinctYears = new TreeSet<>();
        for (NucleotideSequence seq : alignment.getSequences()) {
            if (seq.getCollectionDate().isEmpty()) continue;
            double year = TemporalUtil.toDecimalYear(seq.getCollectionDate().get());
            clockRegression.addData(year, tree.rootToTip(seq.getId()));
            distinctYears.add(year);
        }
        if (distinctYears.size() < 2) {
            throw new GviComputationException("Phylodynamic Re fallback needs at least 2 distinct collection dates, got " + distinctYears.size());
        }
        double mu = clockRegression.getSlope();
        if (mu <= 0) {
            throw new GviComputationException(
                    "Phylodynamic Re fallback: root-to-tip regression gives a non-positive substitution rate (mu=" + mu
                            + "), so branch lengths cannot be converted to calendar time. Temporal/genetic signal may be too weak.");
        }
        double intercept = clockRegression.getIntercept();
        double rootYear = -intercept / mu; // year at which root-to-tip distance extrapolates to 0

        // Step 3+4: convert every branch point's substitution distance to calendar time, build LTT, regress.
        List<PhyloNode> branchPoints = tree.branchPoints();
        if (branchPoints.size() < 3) {
            throw new GviComputationException(
                    "Phylodynamic Re fallback needs at least 3 tree branch points for a meaningful LTT regression, got " + branchPoints.size());
        }
        List<Double> branchTimes = new ArrayList<>();
        for (PhyloNode bp : branchPoints) {
            branchTimes.add(rootYear + bp.rootToTipDistance() / mu);
        }
        branchTimes.sort(Double::compareTo);

        SimpleRegression ltt = new SimpleRegression();
        ltt.addData(rootYear, Math.log(1)); // one lineage at the root itself
        int lineageCount = 1;
        for (double t : branchTimes) {
            lineageCount++;
            ltt.addData(t, Math.log(lineageCount));
        }

        double rPerYear = ltt.getSlope();
        double rSquared = Double.isNaN(ltt.getRSquare()) ? 0.0 : ltt.getRSquare();

        double serialIntervalSd = generationTimeDays * SERIAL_INTERVAL_SD_TO_MEAN_RATIO;
        SerialInterval si = new SerialInterval(generationTimeDays, serialIntervalSd,
                serialIntervalTMax(generationTimeDays, serialIntervalSd));
        double re = reFromGrowthRate(rPerYear, si);

        Double doublingTimeDays = null;
        if (rPerYear > 0) {
            doublingTimeDays = (Math.log(2) / rPerYear) * 365.25;
        }

        diagnostics.add("Phylodynamic fallback (Pybus & Rambaut lineages-through-time method): no incidence series "
                + "supplied. Tree time-scaled with mu=" + String.format("%.6g", mu) + " sub/site/yr (root ~"
                + String.format("%.2f", rootYear) + "); growth rate fit over " + branchPoints.size()
                + " branch points. Assumes a strict molecular clock and roughly constant sampling proportion over "
                + "time -- treat as an approximation, not a substitute for case-incidence-based Re or full "
                + "coalescent/birth-death phylodynamic inference (e.g. BEAST).");
        diagnostics.add(String.format("Growth rate r=%.4f/yr converted to Re via the Euler-Lotka renewal equation "
                        + "(Wallinga & Lipsitch 2007) against a Gamma(mean=%.1f, sd=%.1f) serial interval -- generalizes "
                        + "the older linear Re~=1+r*generationTime approximation (exact for r=0, increasingly more "
                        + "accurate than it for faster growth/decline).",
                rPerYear, generationTimeDays, generationTimeDays * SERIAL_INTERVAL_SD_TO_MEAN_RATIO));
        if (rSquared < 0.3) {
            diagnostics.add(String.format("Weak lineages-through-time signal (R^2=%.3f < 0.3): this Re estimate is unreliable", rSquared));
        }

        var band = ReReferenceTable.classify(re);
        String category = band.category() + " (" + band.transmission() + "); " + band.controlMeasure();

        LocalDate lastDate = alignment.getSequences().stream()
                .map(NucleotideSequence::getCollectionDate).filter(java.util.Optional::isPresent).map(java.util.Optional::get)
                .max(LocalDate::compareTo).orElse(LocalDate.now());

        return new ReResult(ReMethod.PHYLODYNAMIC_FALLBACK, lastDate, re, null, List.of(), doublingTimeDays, category, diagnostics);
    }

    /**
     * Euler-Lotka renewal equation (Wallinga &amp; Lipsitch 2007): for
     * exponential growth I(t) = I0*e^(r*t) under the renewal equation
     * I(t) = Re * integral(I(t-a) * w(a) da), substituting and cancelling
     * I0*e^(r*t) gives 1 = Re * integral(e^(-r*a) * w(a) da) --
     * discretized over the serial interval's daily weights, this is
     * Re = 1 / sum_u( w(u) * e^(-r*u) ). r is converted from per-year (the
     * LTT regression's units) to per-day (the serial interval's units)
     * first. Verified: r=0 gives Re=1 exactly (Sum w(u)=1), and a
     * first-order Taylor expansion for small r reduces to the classic
     * Re ~= 1 + r*meanSerialInterval approximation.
     */
    static double reFromGrowthRate(double rPerYear, SerialInterval si) {
        double rPerDay = rPerYear / 365.25;
        double denominator = 0.0;
        for (int u = 1; u <= si.tMax(); u++) {
            denominator += si.weight(u) * Math.exp(-rPerDay * u);
        }
        return 1.0 / denominator;
    }
}
