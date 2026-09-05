package org.gvi.algorithms.re;

import java.util.List;

/** Reference thresholds from build spec Section 2.4. */
public final class ReReferenceTable {

    public record Band(double upperBound, String category, String transmission, String doublingTime,
                        String casesPer100k, String controlMeasure) {
    }

    public static final List<Band> BANDS = List.of(
            new Band(0.8, "Controlled", "Declining -> Extinction", "n/a", "Decreasing", "Maintenance; surveillance"),
            new Band(1.2, "Endemic", "Stable circulation -> Stable", "n/a", "Constant low", "Routine vaccination"),
            new Band(1.5, "Moderate growth", "Expanding", "Days to weeks", "Exponential", "Targeted NPIs; enhanced surveillance"),
            new Band(2.0, "High growth", "Rapid expansion", "3-7 days", "Rapid exponential", "Strong NPIs; mass vaccination push"),
            new Band(Double.POSITIVE_INFINITY, "Pandemic", "Explosive growth", "2-3 days", "Healthcare overwhelmed", "Full lockdown; emergency measures")
    );

    private ReReferenceTable() {
    }

    public static Band classify(double re) {
        for (Band b : BANDS) {
            if (re <= b.upperBound()) return b;
        }
        return BANDS.get(BANDS.size() - 1);
    }
}
