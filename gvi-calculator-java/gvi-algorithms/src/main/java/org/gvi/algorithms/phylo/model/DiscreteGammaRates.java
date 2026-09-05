package org.gvi.algorithms.phylo.model;

import org.apache.commons.math3.distribution.GammaDistribution;
import org.gvi.core.exception.GviInputException;

/**
 * Discretized Gamma among-site rate heterogeneity (Yang 1994) -- the "+G"
 * in "GTR+Gamma". Real sequence data virtually never evolves at a single
 * uniform rate across all sites (structurally/functionally constrained
 * sites are slower, others faster); modeling that with a single rate, as
 * the plain GTR/JC69/K80 distances elsewhere in this codebase do, biases
 * distance and branch-length estimates downward on divergent sequences.
 * <p>
 * Splits a Gamma(alpha, alpha) distribution (shape=alpha, mean fixed to 1
 * so branch lengths don't need separate rescaling) into {@code k} equal-
 * probability categories and returns each category's mean rate -- the
 * standard "mean of category" discretization (as opposed to the coarser
 * "median of category" method), matching PAML/RAxML/IQ-TREE's convention.
 */
public final class DiscreteGammaRates {

    private DiscreteGammaRates() {
    }

    public static double[] categories(double alpha, int k) {
        if (alpha <= 0) {
            throw new GviInputException("Gamma shape parameter alpha must be positive, got " + alpha);
        }
        if (k < 1) {
            throw new GviInputException("Number of rate categories must be >= 1, got " + k);
        }
        if (k == 1) {
            return new double[]{1.0}; // no heterogeneity: a single category at the mean rate
        }

        double scale = 1.0 / alpha; // shape=alpha, scale=1/alpha -> mean = alpha*scale = 1
        GammaDistribution shapeAlpha = new GammaDistribution(alpha, scale);
        GammaDistribution shapeAlphaPlus1 = new GammaDistribution(alpha + 1.0, scale);

        double[] boundaries = new double[k + 1];
        boundaries[0] = 0.0;
        for (int i = 1; i < k; i++) {
            boundaries[i] = shapeAlpha.inverseCumulativeProbability((double) i / k);
        }
        boundaries[k] = Double.POSITIVE_INFINITY;

        double[] rates = new double[k];
        for (int i = 0; i < k; i++) {
            double cdfLower = boundaries[i] == 0.0 ? 0.0 : shapeAlphaPlus1.cumulativeProbability(boundaries[i]);
            double cdfUpper = Double.isInfinite(boundaries[i + 1]) ? 1.0 : shapeAlphaPlus1.cumulativeProbability(boundaries[i + 1]);
            rates[i] = k * (cdfUpper - cdfLower);
        }
        return rates;
    }
}
