package org.gvi.algorithms.dnds;

/**
 * dN/dS pooled over one codon window (1-based, inclusive codon positions
 * within the gene) -- gives site-level resolution instead of one
 * gene-wide ratio, so a small number of positively-selected codons isn't
 * diluted into invisibility by surrounding purifying selection.
 */
public record SlidingWindowDnDsResult(String sequenceId, String geneName, int startCodon, int endCodon,
                                       double dN, double dS, double omega, double synonymousSites,
                                       double nonsynonymousSites, double synonymousDifferences,
                                       double nonsynonymousDifferences, double selectionZScore,
                                       double selectionPValue, NeiGojoboriSelectionTest.Verdict selectionVerdict)
        implements DnDsCounts {
}
