package org.gvi.algorithms.dnds;

/**
 * The raw Nei-Gojobori site/difference counts behind a dN/dS ratio --
 * implemented by both {@link DnDsResult} (whole-gene/pooled) and
 * {@link SlidingWindowDnDsResult} (per-window), so {@link NeiGojoboriSelectionTest}
 * can test either for statistical significance with one method.
 */
public interface DnDsCounts {
    double synonymousSites();

    double nonsynonymousSites();

    double synonymousDifferences();

    double nonsynonymousDifferences();
}
