package org.gvi.algorithms.cai;

import java.util.List;

/** CAI bands from build spec Section 8.4. */
public final class CaiReferenceTable {

    public record Band(double upperBound, String status, String interpretation, String evolutionaryAge) {
    }

    public static final List<Band> BANDS = List.of(
            new Band(0.50, "Unadapted", "Recent spillover / zoonotic", "Days-weeks"),
            new Band(0.70, "Intermediate", "Emerging; partial adaptation", "Weeks-months"),
            new Band(0.80, "Adapted", "Established circulation", "Months-years"),
            new Band(Double.POSITIVE_INFINITY, "Optimized", "Long-term human adaptation", "Years-decades")
    );

    private CaiReferenceTable() {
    }

    public static Band classify(double cai) {
        for (Band b : BANDS) {
            if (cai <= b.upperBound()) return b;
        }
        return BANDS.get(BANDS.size() - 1);
    }
}
