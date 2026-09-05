package org.gvi.core.spi;

import java.util.List;

/**
 * Minimal reusable {@link IndexResult} for ad-hoc/aggregated values that
 * don't warrant their own dedicated record type -- e.g. a composite-facing
 * summary (max/mean) built from several per-gene results.
 */
public record SimpleIndexResult(String indexName, double primaryValue, String category,
                                 List<String> diagnostics) implements IndexResult {
}
