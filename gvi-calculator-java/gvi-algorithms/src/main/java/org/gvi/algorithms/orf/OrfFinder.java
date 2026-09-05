package org.gvi.algorithms.orf;

import org.gvi.core.model.GeneAnnotation;
import org.gvi.core.util.GeneticCode;

import java.util.ArrayList;
import java.util.List;

/**
 * Native open-reading-frame prediction, used as a real fallback when no
 * {@code --gff} annotation is supplied instead of blindly treating the
 * whole sequence as one ORF. Standard 6-frame scan (3 forward, 3
 * reverse-complement): within each frame, the sequence is split into
 * "stop-to-stop" segments (regions between consecutive in-frame stop
 * codons, or from a sequence end to the nearest stop) -- the same
 * convention NCBI's ORFfinder and similar tools use by default -- and each
 * segment's ORF is defined as running from its FIRST ATG through to the
 * terminal stop codon (the longest possible ORF in that segment), which
 * avoids reporting a cascade of nested candidate ORFs from downstream
 * in-frame ATGs within one real gene.
 * <p>
 * Deliberately scoped: this identifies ORF *coordinates* from sequence
 * signal alone (start/stop codon positions), not genuine ab initio gene
 * prediction (no coding-potential model, no splice-site handling, no
 * comparison against known protein families the way Prodigal/GeneMark do)
 * -- appropriate for compact, intron-free viral genomes, where an ORF
 * bounded by in-frame stop codons is usually a real gene, but a much
 * simpler technique than a trained gene-finder.
 */
public final class OrfFinder {

    /** NCBI ORFfinder's own default minimum is 75nt (25 codons); this project's default is a little more permissive to avoid missing small real viral accessory ORFs. */
    public static final int DEFAULT_MIN_ORF_CODONS = 30;

    private final int minOrfCodons;

    public OrfFinder() {
        this(DEFAULT_MIN_ORF_CODONS);
    }

    public OrfFinder(int minOrfCodons) {
        this.minOrfCodons = minOrfCodons;
    }

    public List<GeneAnnotation> findOrfs(String sequence) {
        String upper = sequence.toUpperCase();
        List<GeneAnnotation> orfs = new ArrayList<>();

        for (int frame = 0; frame < 3; frame++) {
            orfs.addAll(scanFrame(upper, frame, '+'));
        }

        String revComp = reverseComplement(upper);
        int length = upper.length();
        for (int frame = 0; frame < 3; frame++) {
            for (GeneAnnotation orf : scanFrame(revComp, frame, '-')) {
                // convert reverse-complement-relative 1-based coordinates back to the original (forward) sequence
                long originalStart = length - orf.end() + 1;
                long originalEnd = length - orf.start() + 1;
                orfs.add(new GeneAnnotation(orf.geneName(), originalStart, originalEnd, '-'));
            }
        }

        orfs.sort((a, b) -> Long.compare(b.length(), a.length())); // longest first, a reasonable default presentation order
        return orfs;
    }

    private List<GeneAnnotation> scanFrame(String seq, int frame, char strand) {
        List<GeneAnnotation> found = new ArrayList<>();
        int n = 0;
        Integer firstAtg = null; // 0-based index of the first ATG seen in the current stop-to-stop segment, if any

        for (int pos = frame; pos + 3 <= seq.length(); pos += 3) {
            String codon = seq.substring(pos, pos + 3);
            if (firstAtg == null && codon.equals("ATG")) {
                firstAtg = pos;
            }
            if (isStop(codon)) {
                if (firstAtg != null) {
                    int orfLengthCodons = (pos + 3 - firstAtg) / 3;
                    if (orfLengthCodons >= minOrfCodons) {
                        n++;
                        found.add(new GeneAnnotation((strand == '+' ? "orf_fwd" : "orf_rev") + frame + "_" + n,
                                firstAtg + 1L, pos + 3L, strand));
                    }
                }
                firstAtg = null; // start scanning a new stop-to-stop segment
            }
        }
        return found;
    }

    private boolean isStop(String codon) {
        if (codon.indexOf('N') >= 0 || codon.indexOf('-') >= 0) return false;
        try {
            return GeneticCode.translate(codon) == GeneticCode.STOP;
        } catch (RuntimeException e) {
            return false; // non-standard/ambiguous codon: treat as neither start nor stop rather than crash the scan
        }
    }

    private String reverseComplement(String seq) {
        StringBuilder sb = new StringBuilder(seq.length());
        for (int i = seq.length() - 1; i >= 0; i--) {
            sb.append(complement(seq.charAt(i)));
        }
        return sb.toString();
    }

    private char complement(char c) {
        return switch (c) {
            case 'A' -> 'T';
            case 'T', 'U' -> 'A';
            case 'C' -> 'G';
            case 'G' -> 'C';
            default -> 'N';
        };
    }
}
