package org.gvi.algorithms.re;

import org.gvi.core.model.ConfidenceInterval;
import org.gvi.core.spi.IndexResult;

import java.time.LocalDate;
import java.util.List;

/** Result of Index 2 - Effective Reproduction Number (Section 5.2). */
public record ReResult(ReMethod method, LocalDate estimateDate, double reMean, ConfidenceInterval interval,
                        List<ReDailyEstimate> series, Double doublingTimeDays,
                        String category, List<String> diagnostics) implements IndexResult {

    @Override
    public String indexName() {
        return "Re";
    }

    @Override
    public double primaryValue() {
        return reMean;
    }
}
