package org.gvi.algorithms.common;

import org.gvi.core.exception.GviComputationException;

import java.util.ArrayList;
import java.util.List;

/**
 * Splits an in-frame CDS nucleotide string into codons, for the two
 * codon-level algorithms (CAI, dN/dS). A trailing partial codon (length not
 * a multiple of 3) is dropped rather than crashing, with the caller
 * expected to surface {@code droppedTrailingBases} in its diagnostics.
 */
public final class CodonUtil {

    public record SplitResult(List<String> codons, int droppedTrailingBases) {
    }

    private CodonUtil() {
    }

    public static SplitResult splitIntoCodons(String cds) {
        if (cds == null || cds.isEmpty()) {
            throw new GviComputationException("Cannot split an empty CDS into codons");
        }
        int usableLength = (cds.length() / 3) * 3;
        List<String> codons = new ArrayList<>(usableLength / 3);
        for (int i = 0; i < usableLength; i += 3) {
            codons.add(cds.substring(i, i + 3));
        }
        return new SplitResult(codons, cds.length() - usableLength);
    }

    /** Extracts the 1-based inclusive [start,end] region from a full sequence, reverse-complementing if strand is '-'. */
    public static String extractCds(String fullSequence, long start, long end, char strand) {
        if (start < 1 || end > fullSequence.length() || end < start) {
            throw new GviComputationException("Gene coordinates [" + start + "," + end
                    + "] fall outside sequence of length " + fullSequence.length());
        }
        String region = fullSequence.substring((int) start - 1, (int) end);
        return strand == '-' ? reverseComplement(region) : region;
    }

    private static String reverseComplement(String seq) {
        StringBuilder sb = new StringBuilder(seq.length());
        for (int i = seq.length() - 1; i >= 0; i--) {
            sb.append(complement(seq.charAt(i)));
        }
        return sb.toString();
    }

    private static char complement(char c) {
        return switch (Character.toUpperCase(c)) {
            case 'A' -> 'T';
            case 'T', 'U' -> 'A';
            case 'C' -> 'G';
            case 'G' -> 'C';
            case 'R' -> 'Y';
            case 'Y' -> 'R';
            case 'W' -> 'W';
            case 'S' -> 'S';
            case 'K' -> 'M';
            case 'M' -> 'K';
            case 'B' -> 'V';
            case 'V' -> 'B';
            case 'D' -> 'H';
            case 'H' -> 'D';
            default -> 'N';
        };
    }
}
