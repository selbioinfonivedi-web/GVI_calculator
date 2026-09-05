package org.gvi.core.model;

/** A generic two-sided interval, e.g. a bootstrap or Bayesian credible interval. {@code level} is e.g. 0.95. */
public record ConfidenceInterval(double lower, double upper, double level) {

    public ConfidenceInterval {
        if (level <= 0 || level >= 1) {
            throw new IllegalArgumentException("confidence level must be in (0,1), got " + level);
        }
    }
}
