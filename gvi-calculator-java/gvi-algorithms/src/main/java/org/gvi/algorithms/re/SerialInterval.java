package org.gvi.algorithms.re;

import org.apache.commons.math3.distribution.GammaDistribution;
import org.gvi.core.exception.GviInputException;

/**
 * Discretized serial-interval distribution w(u), u=1..tMax days, used as
 * the infectiousness profile in the Cori et al. (2013) Re estimator
 * (Section 5.2a). Modeled as a Gamma(mean, sd) distribution sampled at
 * integer days and renormalized to sum to 1 -- a standard simplified
 * discretization (the exact Cori/Nishiura integral form is more intricate;
 * this density-sampling approximation is accurate to within a fraction of
 * a percent for mean/sd in the ranges typical of respiratory pathogens and
 * is far less error-prone to implement correctly).
 */
public final class SerialInterval {

    private final double[] weights; // index 0 == day 1

    public SerialInterval(double meanDays, double sdDays, int tMax) {
        if (meanDays <= 0 || sdDays <= 0) {
            throw new GviInputException("Serial interval mean and sd must both be positive (got mean=" + meanDays + ", sd=" + sdDays + ")");
        }
        double shape = (meanDays * meanDays) / (sdDays * sdDays);
        double scale = (sdDays * sdDays) / meanDays;
        GammaDistribution gamma = new GammaDistribution(shape, scale);

        double[] raw = new double[tMax];
        double sum = 0;
        for (int u = 1; u <= tMax; u++) {
            raw[u - 1] = gamma.density(u);
            sum += raw[u - 1];
        }
        if (sum <= 0) {
            throw new GviInputException("Serial interval distribution has no mass in [1," + tMax + "] days -- widen tMax or check mean/sd");
        }
        this.weights = new double[tMax];
        for (int i = 0; i < tMax; i++) weights[i] = raw[i] / sum;
    }

    /** w(u) for u in [1, tMax], 0 outside that range. */
    public double weight(int u) {
        if (u < 1 || u > weights.length) return 0.0;
        return weights[u - 1];
    }

    public int tMax() {
        return weights.length;
    }

    public static SerialInterval defaultProfile() {
        // Default per build spec Section 5.2: pathogen-agnostic Gamma(mean=5, sd=2) days.
        return new SerialInterval(5.0, 2.0, 20);
    }
}
