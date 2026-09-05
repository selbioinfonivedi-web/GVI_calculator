package org.gvi.algorithms.cai;

import java.util.List;

/** GC content deviation bands from build spec Section 8.4. */
public final class GcReferenceTable {

    public record Band(double upperBound, String status, String interpretation, String evolutionaryAge) {
    }

    public static final List<Band> BANDS = List.of(
            new Band(1.0, "Native", "No HGT signal", "Stable, expected"),
            new Band(3.0, "Minor drift", "Minor recombination / reassortment", "Possible but not recent"),
            new Band(Double.POSITIVE_INFINITY, "Anomaly", "HGT suspected", "Recent foreign acquisition")
    );

    private GcReferenceTable() {
    }

    public static Band classify(double gcDeviationPercent) {
        for (Band b : BANDS) {
            if (gcDeviationPercent <= b.upperBound()) return b;
        }
        return BANDS.get(BANDS.size() - 1);
    }
}
