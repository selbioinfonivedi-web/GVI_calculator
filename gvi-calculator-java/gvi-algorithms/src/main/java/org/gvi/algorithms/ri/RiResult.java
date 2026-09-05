package org.gvi.algorithms.ri;

import org.gvi.core.spi.IndexResult;

import java.util.List;

/** Result of Index 7 - Recombination Index (Section 5.7). */
public record RiResult(RiMethod method, double ri, Double windowedPhiStatistic, Double globalPhiStatistic,
                        Double pValue, Integer informativeSitesUsed, Integer permutationsRun,
                        String category, List<String> diagnostics) implements IndexResult {

    @Override
    public String indexName() {
        return "RI";
    }

    @Override
    public double primaryValue() {
        return ri;
    }

    public boolean significant() {
        return pValue != null && pValue < 0.05;
    }
}
