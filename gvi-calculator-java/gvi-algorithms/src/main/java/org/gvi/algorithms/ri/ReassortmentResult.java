package org.gvi.algorithms.ri;

import java.util.List;

/**
 * Result of testing two independently-aligned genome segments (e.g. Bluetongue's Segment 2 vs Segment
 * 10, influenza's HA vs NA) for reassortment -- a distinct biological process from within-locus
 * recombination ({@link RecombinationIndexCalculator}'s PHI test), which only applies within one
 * contiguous alignment and cannot detect segments swapping between co-infecting strains.
 */
public record ReassortmentResult(String segmentALabel, String segmentBLabel, int taxaCompared,
                                  double mantelR, double pValue, double reassortmentIndex,
                                  String category, List<String> diagnostics) {

    public boolean significant() {
        return pValue < 0.05;
    }
}
