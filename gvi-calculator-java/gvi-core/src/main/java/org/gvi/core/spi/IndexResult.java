package org.gvi.core.spi;

import java.util.List;

/**
 * Contract every one of the 8 GVI component index results implements, so
 * gvi-composite can normalize/weight/combine them uniformly (Section 6:
 * each index is independently computable AND combinable). Individual
 * modules return their own richer record (extra fields, per-gene
 * breakdowns, confidence intervals) that also implements this interface.
 */
public interface IndexResult {

    /** Canonical short name, e.g. "mu", "Re", "pi", "MB", "dN/dS", "GD", "RI", "CAI". */
    String indexName();

    /**
     * The single scalar this index contributes to composite normalization
     * (Section 5.9). What it means is index-specific (a rate, a ratio, a
     * count, a distance) -- gvi-composite maps it to [0,1] using the
     * configured reference range for this index.
     */
    double primaryValue();

    /** Reference-table classification for this value (e.g. "RNA (low-fidelity)"), or "N/A" if not applicable. */
    String category();

    /** Non-fatal caveats about this result (e.g. low temporal signal, distance saturated, estimate from subsampling). */
    List<String> diagnostics();
}
