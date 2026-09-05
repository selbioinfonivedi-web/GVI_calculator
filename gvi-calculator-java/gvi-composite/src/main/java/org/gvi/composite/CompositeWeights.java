package org.gvi.composite;

import org.gvi.core.exception.GviInputException;

import java.util.EnumMap;
import java.util.Map;

/**
 * Per-index weights for the composite GVI (Section 5.9), Sigma(w) = 1.
 * Defaults are the midpoints of the spec's typical-weight ranges, with the
 * combined "CAI+GC" band (0.03-0.08, midpoint 0.055) split evenly between
 * CAI and GC Content Deviation.
 */
public final class CompositeWeights {

    private static final double TOLERANCE = 1e-6;

    private final Map<IndexKey, Double> weights;

    private CompositeWeights(Map<IndexKey, Double> weights) {
        double sum = weights.values().stream().mapToDouble(Double::doubleValue).sum();
        if (Math.abs(sum - 1.0) > TOLERANCE) {
            throw new GviInputException("Composite weights must sum to 1.0, got " + sum);
        }
        for (double w : weights.values()) {
            if (w < 0) {
                throw new GviInputException("Composite weights must be non-negative");
            }
        }
        this.weights = new EnumMap<>(weights);
    }

    public static CompositeWeights of(Map<IndexKey, Double> weights) {
        return new CompositeWeights(weights);
    }

    public static CompositeWeights defaults() {
        // Midpoints of each index's typical-weight range from Section 5.9's table.
        // Taken independently per index, these midpoints sum to 1.09, not 1.0 (the
        // ranges are individually "typical", not a strict partition) -- so they are
        // rescaled proportionally here to satisfy Sigma(w)=1 while preserving their
        // relative importance exactly as specified. CAI+GC's combined 0.055 midpoint
        // is split evenly between the two sub-components.
        Map<IndexKey, Double> midpoints = new EnumMap<>(IndexKey.class);
        midpoints.put(IndexKey.MU, 0.135);
        midpoints.put(IndexKey.RE, 0.30);
        midpoints.put(IndexKey.PI, 0.10);
        midpoints.put(IndexKey.MB, 0.215);
        midpoints.put(IndexKey.DNDS, 0.15);
        midpoints.put(IndexKey.GD, 0.10);
        midpoints.put(IndexKey.CAI, 0.0275);
        midpoints.put(IndexKey.GC, 0.0275);
        midpoints.put(IndexKey.RI, 0.035);

        double sum = midpoints.values().stream().mapToDouble(Double::doubleValue).sum();
        Map<IndexKey, Double> normalized = new EnumMap<>(IndexKey.class);
        for (var entry : midpoints.entrySet()) {
            normalized.put(entry.getKey(), entry.getValue() / sum);
        }
        return new CompositeWeights(normalized);
    }

    public double get(IndexKey key) {
        return weights.getOrDefault(key, 0.0);
    }

    public Map<IndexKey, Double> asMap() {
        return new EnumMap<>(weights);
    }
}
