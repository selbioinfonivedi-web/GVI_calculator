package org.gvi.core.model;

import org.gvi.core.exception.GviInputException;
import org.gvi.core.util.GeneticCode;

import java.util.HashMap;
import java.util.Map;

/**
 * Host reference codon-usage table used by CAI (Section 5.8): raw codon
 * counts/frequencies in, per-codon relative synonymous codon usage (RSCU)
 * weight w(codon) precomputed out, where w(codon) = RSCU(codon) /
 * RSCU(max synonymous codon for that amino acid), per Sharp & Li 1987.
 */
public final class CodonUsageTable {

    private final Map<String, Double> weights = new HashMap<>();
    private final String name;

    private CodonUsageTable(String name, Map<String, Double> rawFrequencies) {
        this.name = name;
        computeWeights(rawFrequencies);
    }

    public static CodonUsageTable fromFrequencies(String name, Map<String, Double> codonToFrequency) {
        if (codonToFrequency == null || codonToFrequency.isEmpty()) {
            throw new GviInputException("Codon usage table '" + name + "' is empty");
        }
        return new CodonUsageTable(name, codonToFrequency);
    }

    private void computeWeights(Map<String, Double> raw) {
        for (Map.Entry<Character, java.util.List<String>> aa : GeneticCode.SYNONYMOUS_CODONS.entrySet()) {
            if (aa.getKey() == '*') continue; // stop codons excluded from CAI per convention
            double max = 0.0;
            for (String codon : aa.getValue()) {
                max = Math.max(max, raw.getOrDefault(codon, 0.0));
            }
            for (String codon : aa.getValue()) {
                double freq = raw.getOrDefault(codon, 0.0);
                double w = (max <= 0.0) ? 0.0 : freq / max;
                // Avoid exact-zero weights producing -Infinity in the CAI log-geometric-mean;
                // a codon truly never observed in a whole-genome reference table gets a small floor.
                weights.put(codon, w <= 0.0 ? 0.0001 : w);
            }
        }
    }

    public double weightOf(String codon) {
        Double w = weights.get(codon.toUpperCase());
        if (w == null) {
            throw new GviInputException("Codon '" + codon + "' not found in codon usage table '" + name + "'");
        }
        return w;
    }

    public String getName() {
        return name;
    }
}
