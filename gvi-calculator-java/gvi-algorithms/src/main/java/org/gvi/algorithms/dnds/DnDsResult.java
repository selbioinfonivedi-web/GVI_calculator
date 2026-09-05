package org.gvi.algorithms.dnds;

import org.gvi.core.spi.IndexResult;

import java.util.List;

/** Result of Index 5 - dN/dS Ratio (Section 5.5), via the Nei-Gojobori (1986) counting method. */
public record DnDsResult(String sequenceId, String geneName, double dN, double dS, double omega,
                          boolean jcCorrectionApplied, double synonymousSites, double nonsynonymousSites,
                          double synonymousDifferences, double nonsynonymousDifferences,
                          int codonsCompared, int excludedCodons, String category, List<String> diagnostics)
        implements IndexResult, DnDsCounts {

    @Override
    public String indexName() {
        return "dN/dS";
    }

    @Override
    public double primaryValue() {
        return omega;
    }
}
