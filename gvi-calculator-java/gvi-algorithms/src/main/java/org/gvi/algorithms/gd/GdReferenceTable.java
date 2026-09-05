package org.gvi.algorithms.gd;

import java.util.List;

/**
 * Reference thresholds & interpretation table from build spec Section 6.4,
 * defined in Jukes-Cantor distance units. Bundled as data (not scattered
 * magic numbers) so it can be reviewed by the virology domain team and
 * swapped without touching calculation code.
 */
public final class GdReferenceTable {

    public record Band(double upperBound, String timeSinceDivergence, String interpretation,
                        String classification, String confidence) {
    }

    public static final List<Band> BANDS = List.of(
            new Band(0.001, "<1 week", "Nearly identical", "Same transmission event", "Very high"),
            new Band(0.005, "1-4 weeks", "Closely related", "Linked outbreak cluster", "High"),
            new Band(0.015, "1-3 months", "Distinct lineage", "Separate introduction event", "Moderate"),
            new Band(0.05, "3-12 months", "Divergent", "Different geographic source", "Low-moderate"),
            new Band(Double.POSITIVE_INFINITY, ">1 year", "Highly divergent", "Different strain / serotype", "Low (ambiguous without tree)")
    );

    private GdReferenceTable() {
    }

    public static Band classify(double jcDistance) {
        for (Band b : BANDS) {
            if (jcDistance <= b.upperBound()) return b;
        }
        return BANDS.get(BANDS.size() - 1);
    }
}
