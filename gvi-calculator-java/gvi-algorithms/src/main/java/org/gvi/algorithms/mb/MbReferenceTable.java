package org.gvi.algorithms.mb;

import java.util.List;

/** Reference thresholds from build spec Section 4.4. */
public final class MbReferenceTable {

    public record Band(int upperBound, String pathogenExample, String timeline, String characteristics, String controlImplication) {
    }

    public static final List<Band> BANDS = List.of(
            new Band(5, "Within-household transmission", "Days", "Near-identical to reference", "Genetic link confirmable"),
            new Band(15, "Community chain (1-2 weeks)", "1-2 weeks", "Early divergence", "Traceable contact history"),
            new Band(40, "Established outbreak (months)", "1-3 months", "Recognizable variant signature", "Variant of interest (VOI) classification"),
            new Band(100, "Pandemic variant (6-18 months)", "6-18 months", "Multiple functional changes likely", "Variant of concern (VOC); major public health action"),
            new Band(Integer.MAX_VALUE, "Beyond reference table", "> 18 months / possible novel lineage", "Highly diverged", "Full genomic + epidemiological re-investigation warranted")
    );

    private MbReferenceTable() {
    }

    public static Band classify(int mutationBurden) {
        for (Band b : BANDS) {
            if (mutationBurden <= b.upperBound()) return b;
        }
        return BANDS.get(BANDS.size() - 1);
    }
}
