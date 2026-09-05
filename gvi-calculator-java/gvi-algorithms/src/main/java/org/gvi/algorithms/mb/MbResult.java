package org.gvi.algorithms.mb;

import org.gvi.core.spi.IndexResult;

import java.util.List;

/** Result of Index 4 - Mutation Burden (Section 5.4). */
public record MbResult(String sequenceId, int mutationBurden, MbSource source,
                        int excludedLowCoverage, int excludedAmbiguousOrUnknownCoverage,
                        String category, List<String> diagnostics) implements IndexResult {

    public enum MbSource { VCF_WITH_COVERAGE_FILTER, ALIGNMENT_NO_COVERAGE_FILTER }

    @Override
    public String indexName() {
        return "MB";
    }

    @Override
    public double primaryValue() {
        return mutationBurden;
    }
}
