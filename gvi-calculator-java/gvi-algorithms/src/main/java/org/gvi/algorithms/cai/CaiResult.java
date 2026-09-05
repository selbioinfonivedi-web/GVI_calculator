package org.gvi.algorithms.cai;

import org.gvi.core.spi.IndexResult;

import java.util.List;

/** Result of the CAI half of Index 8 (Section 5.8). */
public record CaiResult(String sequenceId, String geneName, double cai, int codonsUsed,
                         int excludedStop, int excludedSingleAminoAcid, int excludedInvalid,
                         String category, List<String> diagnostics) implements IndexResult {

    @Override
    public String indexName() {
        return "CAI";
    }

    @Override
    public double primaryValue() {
        return cai;
    }
}
