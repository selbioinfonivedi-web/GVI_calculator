package org.gvi.algorithms.pi;

import java.util.List;

/** Reference thresholds from build spec Section 3.4. */
public final class PiReferenceTable {

    public record Band(double upperBound, String scenario, String interpretation, String implication) {
    }

    public static final List<Band> BANDS = List.of(
            new Band(0.0005, "Acute outbreak", "Very low diversity; single clone", "Originated from single introduction"),
            new Band(0.005, "Established outbreak", "Moderate diversity; chains of transmission", "Multiple lineages co-circulate; endemic risk"),
            new Band(0.02, "Endemic circulation", "High diversity; years of transmission", "Likely antigenically variable; vaccine updates needed"),
            new Band(Double.POSITIVE_INFINITY, "Hyperendemic / highly recombinant", "Very high diversity; multiple species/types", "Complex epidemiology; multiple control targets")
    );

    private PiReferenceTable() {
    }

    public static Band classify(double pi) {
        for (Band b : BANDS) {
            if (pi <= b.upperBound()) return b;
        }
        return BANDS.get(BANDS.size() - 1);
    }
}
