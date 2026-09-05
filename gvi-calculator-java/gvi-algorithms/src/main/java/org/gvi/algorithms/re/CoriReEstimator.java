package org.gvi.algorithms.re;

import org.apache.commons.math3.distribution.GammaDistribution;
import org.gvi.core.exception.GviComputationException;
import org.gvi.core.model.ConfidenceInterval;
import org.gvi.core.model.IncidencePoint;

import java.time.LocalDate;
import java.time.temporal.ChronoUnit;
import java.util.ArrayList;
import java.util.List;

/**
 * Index 2 - Effective Reproduction Number, incidence-based estimator
 * (Section 5.2a): Cori, Ferguson, Fraser &amp; Cauchemez (2013), the same
 * method underlying the widely-used EpiEstim R package. Closed-form
 * Gamma-Poisson conjugate posterior, no MCMC required.
 */
public final class CoriReEstimator {

    private final double priorShape; // EpiEstim default: mean=5, sd=5 -> shape=1
    private final double priorRate;  // -> rate=0.2
    private final int windowDays;

    public CoriReEstimator() {
        this(1.0, 0.2, 7);
    }

    public CoriReEstimator(double priorShape, double priorRate, int windowDays) {
        this.priorShape = priorShape;
        this.priorRate = priorRate;
        this.windowDays = windowDays;
    }

    public ReResult compute(List<IncidencePoint> incidence, SerialInterval si) {
        if (incidence.size() < windowDays + si.tMax()) {
            throw new GviComputationException("Incidence series too short for Re estimation: need at least "
                    + (windowDays + si.tMax()) + " days of data (window=" + windowDays + " + serial-interval tMax=" + si.tMax()
                    + "), got " + incidence.size());
        }
        LocalDate startDate = incidence.get(0).date();
        long spanDays = ChronoUnit.DAYS.between(startDate, incidence.get(incidence.size() - 1).date()) + 1;
        double[] dailyIncidence = new double[(int) spanDays];
        for (IncidencePoint p : incidence) {
            int idx = (int) ChronoUnit.DAYS.between(startDate, p.date());
            dailyIncidence[idx] += p.newCases();
        }

        List<ReDailyEstimate> series = new ArrayList<>();
        int firstEstimableDay = si.tMax() + windowDays - 1;
        for (int t = firstEstimableDay; t < dailyIncidence.length; t++) {
            double sumI = 0, sumLambda = 0;
            for (int s = t - windowDays + 1; s <= t; s++) {
                sumI += dailyIncidence[s];
                double lambdaS = 0;
                for (int u = 1; u <= si.tMax() && (s - u) >= 0; u++) {
                    lambdaS += dailyIncidence[s - u] * si.weight(u);
                }
                sumLambda += lambdaS;
            }
            if (sumLambda <= 0) continue; // no infectious pressure in window; Re undefined here, skip rather than divide by zero

            double posteriorShape = priorShape + sumI;
            double posteriorRate = priorRate + sumLambda;
            double mean = posteriorShape / posteriorRate;
            GammaDistribution posterior = new GammaDistribution(posteriorShape, 1.0 / posteriorRate);
            double lower = posterior.inverseCumulativeProbability(0.025);
            double upper = posterior.inverseCumulativeProbability(0.975);

            series.add(new ReDailyEstimate(startDate.plusDays(t), mean, new ConfidenceInterval(lower, upper, 0.95)));
        }

        if (series.isEmpty()) {
            throw new GviComputationException("No day in the incidence series had sufficient prior infectious pressure to estimate Re");
        }

        ReDailyEstimate latest = series.get(series.size() - 1);
        var band = ReReferenceTable.classify(latest.mean());
        String category = band.category() + " (" + band.transmission() + "); " + band.controlMeasure();

        List<String> diagnostics = new ArrayList<>();
        diagnostics.add("Cori et al. (2013) sliding-window estimator, window=" + windowDays + " days, "
                + series.size() + " daily estimate(s) produced from " + dailyIncidence.length + " days of incidence data");

        return new ReResult(ReMethod.CORI_INCIDENCE, latest.date(), latest.mean(), latest.credibleInterval(),
                series, null, category, diagnostics);
    }
}
