package org.gvi.composite;

import org.gvi.core.exception.GviComputationException;
import org.gvi.core.spi.IndexResult;

import java.util.ArrayList;
import java.util.EnumMap;
import java.util.List;
import java.util.Map;

/**
 * Combines whichever of the 9 scalar components are available into the
 * composite GVI (Section 5.9). Missing indices (e.g. no incidence data so
 * Re wasn't computed) have their configured weight redistributed
 * proportionally across the remaining available indices rather than being
 * silently treated as zero -- Section 5.9's explicit requirement.
 */
public final class CompositeGviEngine {

    public GviResult compute(Map<IndexKey, ? extends IndexResult> available, CompositeWeights weights) {
        return compute(available, weights, NormalizationRange.defaults());
    }

    public GviResult compute(Map<IndexKey, ? extends IndexResult> available, CompositeWeights weights,
                              Map<IndexKey, NormalizationRange> ranges) {
        List<IndexKey> excluded = new ArrayList<>();
        Map<IndexKey, Double> configuredWeights = new EnumMap<>(IndexKey.class);
        double availableWeightSum = 0.0;

        for (IndexKey key : IndexKey.values()) {
            double w = weights.get(key);
            if (w <= 0.0) continue; // not part of this weighting scheme at all
            if (available.containsKey(key)) {
                configuredWeights.put(key, w);
                availableWeightSum += w;
            } else {
                excluded.add(key);
            }
        }

        if (availableWeightSum <= 0.0) {
            throw new GviComputationException("Cannot compute composite GVI: none of the weighted indices have a computed result available");
        }

        List<GviComponent> components = new ArrayList<>();
        double gvi = 0.0;
        for (var entry : configuredWeights.entrySet()) {
            IndexKey key = entry.getKey();
            double configuredWeight = entry.getValue();
            double effectiveWeight = configuredWeight / availableWeightSum;
            IndexResult result = available.get(key);
            NormalizationRange range = ranges.getOrDefault(key, NormalizationRange.defaults().get(key));
            double normalized = range.normalize(result.primaryValue());
            double contribution = effectiveWeight * normalized;
            gvi += contribution;
            components.add(new GviComponent(key, result.primaryValue(), normalized, configuredWeight, effectiveWeight, contribution));
        }

        List<String> diagnostics = new ArrayList<>();
        if (!excluded.isEmpty()) {
            diagnostics.add("Excluded (no data supplied): " + excluded.stream().map(IndexKey::label).reduce((a, b) -> a + ", " + b).orElse("")
                    + " -- remaining weights renormalized to sum to 1 (was " + String.format("%.3f", availableWeightSum) + " of full scheme)");
        }
        diagnostics.add(components.size() + " of " + (components.size() + excluded.size()) + " weighted indices contributed to this GVI");

        return new GviResult(gvi, components, excluded, availableWeightSum, diagnostics);
    }

    /**
     * Sensitivity analysis (Section 5.9 / 8): perturb each available index's
     * weight by +/-perturbationFraction (default use 0.2 for +/-20%),
     * renormalizing the rest proportionally, and report the resulting GVI
     * range -- a tornado-chart-ready summary of which index the composite
     * is most sensitive to.
     */
    public List<SensitivityResult> sensitivityAnalysis(Map<IndexKey, ? extends IndexResult> available,
                                                         CompositeWeights baseWeights, Map<IndexKey, NormalizationRange> ranges,
                                                         double perturbationFraction) {
        GviResult base = compute(available, baseWeights, ranges);
        List<SensitivityResult> results = new ArrayList<>();

        for (GviComponent comp : base.components()) {
            IndexKey key = comp.key();
            double baseWeight = baseWeights.get(key);
            double low = recomputeWithPerturbedWeight(available, baseWeights, ranges, key, baseWeight * (1 - perturbationFraction));
            double high = recomputeWithPerturbedWeight(available, baseWeights, ranges, key, baseWeight * (1 + perturbationFraction));
            results.add(new SensitivityResult(key, base.gvi(), Math.min(low, high), Math.max(low, high)));
        }
        return results;
    }

    private double recomputeWithPerturbedWeight(Map<IndexKey, ? extends IndexResult> available, CompositeWeights baseWeights,
                                                  Map<IndexKey, NormalizationRange> ranges, IndexKey perturbedKey, double newWeight) {
        Map<IndexKey, Double> perturbed = baseWeights.asMap();
        double oldWeight = perturbed.getOrDefault(perturbedKey, 0.0);
        double delta = newWeight - oldWeight;
        perturbed.put(perturbedKey, newWeight);

        // redistribute -delta proportionally across all other currently-weighted indices so the scheme still sums to 1
        double othersSum = perturbed.entrySet().stream()
                .filter(e -> e.getKey() != perturbedKey)
                .mapToDouble(Map.Entry::getValue).sum();
        if (othersSum > 0) {
            for (IndexKey k : List.copyOf(perturbed.keySet())) {
                if (k == perturbedKey) continue;
                double share = perturbed.get(k) / othersSum;
                perturbed.put(k, perturbed.get(k) - delta * share);
            }
        }

        CompositeWeights adjusted = CompositeWeights.of(normalizeToUnitSum(perturbed));
        return compute(available, adjusted, ranges).gvi();
    }

    private Map<IndexKey, Double> normalizeToUnitSum(Map<IndexKey, Double> weights) {
        double sum = weights.values().stream().mapToDouble(Double::doubleValue).sum();
        Map<IndexKey, Double> normalized = new EnumMap<>(IndexKey.class);
        for (var e : weights.entrySet()) {
            double v = Math.max(0.0, e.getValue() / sum);
            normalized.put(e.getKey(), v);
        }
        // renormalize once more in case clamping to 0 above shifted the sum away from 1
        double total = normalized.values().stream().mapToDouble(Double::doubleValue).sum();
        Map<IndexKey, Double> finalWeights = new EnumMap<>(IndexKey.class);
        for (var e : normalized.entrySet()) {
            finalWeights.put(e.getKey(), e.getValue() / total);
        }
        return finalWeights;
    }
}
