package org.gvi.algorithms.pi;

import org.gvi.core.model.ConfidenceInterval;
import org.gvi.core.spi.IndexResult;

import java.util.List;

/** Result of Index 3 - Nucleotide Diversity (Section 5.3). */
public record PiResult(double pi, int sequenceCount, int alignmentLength, long pairsUsed,
                        boolean subsampled, ConfidenceInterval confidenceInterval,
                        String category, List<String> diagnostics) implements IndexResult {

    @Override
    public String indexName() {
        return "pi";
    }

    @Override
    public double primaryValue() {
        return pi;
    }
}
