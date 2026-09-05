package org.gvi.algorithms.gd;

import org.gvi.core.spi.IndexResult;

import java.util.List;

/**
 * Result of a pairwise genetic-distance computation (Section 5.6).
 * {@code distance} is null when the chosen model's log argument went
 * non-positive (saturated at high divergence) -- callers must fall back to
 * {@code pDistance} rather than propagate a NaN.
 */
public record GdResult(String sequenceIdA, String sequenceIdB, GdMethod method,
                        double pDistance, Double distance, boolean saturated,
                        String category, List<String> diagnostics) implements IndexResult {

    @Override
    public String indexName() {
        return "GD";
    }

    @Override
    public double primaryValue() {
        return saturated || distance == null ? pDistance : distance;
    }
}
