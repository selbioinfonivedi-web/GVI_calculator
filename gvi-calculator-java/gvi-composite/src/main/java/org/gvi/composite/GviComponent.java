package org.gvi.composite;

/** One index's contribution to a computed {@link GviResult}. */
public record GviComponent(IndexKey key, double rawValue, double normalizedValue,
                            double configuredWeight, double effectiveWeight, double contribution) {
}
