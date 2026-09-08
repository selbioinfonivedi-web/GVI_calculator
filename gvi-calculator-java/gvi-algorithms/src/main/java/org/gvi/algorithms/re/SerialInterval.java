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
 *
 * <h2>Do not "simplify" this to CDF differencing</h2>
 * That accuracy claim is about w(u) itself, and should not be read as meaning the
 * choice of discretization is immaterial downstream -- it is not. Sampling the
 * density at integer days preserves the requested mean (a Gamma(5, 2) profile
 * discretizes to a mean of 5.0002 days). Differencing the CDF instead,
 * {@code F(u) - F(u-1)}, looks more principled -- it integrates the density over
 * each whole day rather than sampling a point -- but assigns each interval's mass
 * to its right endpoint and shifts that mean to 5.50 days.
 * <p>
 * Because Re is estimated by dividing incidence by this profile, a half-day shift
 * propagates straight into Re, and it <em>amplifies with Re</em>. Measured on a
 * simulated renewal process: 0.008 at Re = 1.15, 0.127 at Re = 2.0, 0.209 at
 * Re = 2.5. The bias is therefore invisible when an outbreak is flat and largest
 * exactly when it is growing fastest -- when the number is being relied on most.
 * {@code SerialIntervalTest} pins the mean-preserving property against precisely
 * this substitution.
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
