package org.gvi.cli;

import com.fasterxml.jackson.databind.ObjectMapper;
import org.gvi.composite.CompositeWeights;
import org.gvi.composite.IndexKey;
import org.gvi.core.exception.GviInputException;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.EnumMap;
import java.util.List;
import java.util.Map;

/**
 * Reads a user-supplied composite-weights JSON file for a real run --
 * closes the loop with {@code --calibrate --weights-out}, whose output is
 * exactly the format this reads: a flat JSON object of index label -&gt;
 * weight, e.g. {@code {"Re": 0.30, "MB": 0.22, ...}} (see {@link IndexKey#label()}
 * for the exact 9 labels: mu, Re, pi, MB, dN/dS, GD, CAI, GC_Deviation, RI).
 * <p>
 * Deliberately permissive rather than requiring a complete, already-
 * normalized set: any index NOT mentioned keeps {@link CompositeWeights#defaults()}'s
 * spec-midpoint weight for it, and the full 9-value set (yours plus
 * whatever defaults filled in) is renormalized to sum to 1.0 at the end --
 * so a user who only wants to say "make Re matter more" can supply just
 * {@code {"Re": 0.5}} without having to also re-specify the other 8.
 * An unrecognized label is a warning, not a fatal error -- one typo
 * shouldn't block the whole run when spec defaults are a safe fallback.
 */
public final class CompositeWeightsReader {

    public record Result(CompositeWeights weights, List<String> warnings) {
    }

    private CompositeWeightsReader() {
    }

    public static Result read(Path path) {
        if (!Files.isReadable(path)) {
            throw new GviInputException("Cannot read composite weights file: " + path);
        }
        Map<String, Object> raw;
        try {
            raw = new ObjectMapper().readValue(path.toFile(), Map.class);
        } catch (IOException e) {
            throw new GviInputException("Malformed composite weights JSON '" + path + "': " + e.getMessage(), e);
        }

        List<String> warnings = new ArrayList<>();
        Map<IndexKey, Double> merged = new EnumMap<>(CompositeWeights.defaults().asMap());

        for (var entry : raw.entrySet()) {
            IndexKey key;
            try {
                key = IndexKey.fromLabel(entry.getKey());
            } catch (IllegalArgumentException e) {
                warnings.add("unrecognized index label '" + entry.getKey() + "' ignored (" + e.getMessage() + ")");
                continue;
            }
            double value;
            try {
                value = ((Number) entry.getValue()).doubleValue();
            } catch (ClassCastException | NullPointerException e) {
                warnings.add("weight for '" + entry.getKey() + "' is not a number, ignored (kept the spec-default weight instead)");
                continue;
            }
            if (value < 0) {
                warnings.add("weight for '" + entry.getKey() + "' is negative (" + value + "), ignored (kept the spec-default weight instead)");
                continue;
            }
            merged.put(key, value);
        }

        double sum = merged.values().stream().mapToDouble(Double::doubleValue).sum();
        if (sum <= 0) {
            throw new GviInputException("Composite weights in '" + path + "' sum to " + sum + " after merging with defaults; cannot normalize");
        }
        Map<IndexKey, Double> normalized = new EnumMap<>(IndexKey.class);
        for (var entry : merged.entrySet()) {
            normalized.put(entry.getKey(), entry.getValue() / sum);
        }

        return new Result(CompositeWeights.of(normalized), warnings);
    }
}
