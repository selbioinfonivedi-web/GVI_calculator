package org.gvi.core.util;

import org.gvi.core.model.NucleotideSequence;
import org.gvi.core.model.SequenceAlignment;

import java.util.ArrayList;
import java.util.List;

/**
 * Trims an alignment to the columns every sequence actually covers.
 * <p>
 * <b>The problem this solves.</b> A set of partial sequences with different coverage windows --
 * different amplicons, or the same gene sequenced to different lengths -- aligns correctly, but the
 * shorter ones are gap-padded at the ends. The alignment is not wrong; it simply contains columns
 * where most sequences have no data. Every position-wise index then compares real bases against
 * padding: mutation burden counts the padding as indels, distances are computed over positions only
 * some sequences observed, and the whole set trips the gap-fraction quality gate.
 * <p>
 * Diagnosing this correctly matters, because the instinctive remedy is wrong. Re-running MAFFT does
 * nothing -- the aligner already placed these sequences correctly, and the gaps are a faithful
 * representation of missing coverage, not a misalignment. What the indices need is the shared
 * window: the columns where every sequence has a real base, so that "difference" means an observed
 * difference rather than absent data.
 * <p>
 * On this project's corpus, trimming recovered mutation burden, nucleotide diversity, dN/dS and
 * genetic distance for datasets where the gap fraction ran to 20-29%, retaining roughly 40-47% of
 * columns. It does not rescue an alignment whose sequences are genuinely too divergent to be one
 * gene from one organism -- that is a different failure, and the trimmed result still shows it.
 */
public final class AlignmentTrimmer {

    /** Outcome of a trim, including what it cost, so the caller can report it honestly. */
    public record Result(SequenceAlignment alignment, int originalColumns, int retainedColumns,
                         int longestContiguousBlock, boolean trimmed, String description,
                         int[] keptColumns) {

        public double retainedFraction() {
            return originalColumns == 0 ? 0 : (double) retainedColumns / originalColumns;
        }

        /**
         * Remaps a 1-based coordinate range from the original alignment onto the trimmed one, or
         * empty when the range no longer survives.
         * <p>
         * Trimming renumbers every column, so gene annotations supplied against the ORIGINAL
         * alignment become meaningless afterwards -- a CDS at [1, 1146] on a 456-column trimmed
         * alignment simply fails to load, silently costing the dN/dS and CAI that trimming was
         * meant to rescue. Coordinates must travel with the trim.
         */
        public java.util.Optional<long[]> remap(long start1Based, long end1Based) {
            if (!trimmed) return java.util.Optional.of(new long[]{start1Based, end1Based});
            int newStart = -1;
            int newEnd = -1;
            for (int i = 0; i < keptColumns.length; i++) {
                long original = keptColumns[i] + 1L; // keptColumns is 0-based
                if (original >= start1Based && original <= end1Based) {
                    if (newStart < 0) newStart = i + 1;
                    newEnd = i + 1;
                }
            }
            return newStart < 0 ? java.util.Optional.empty() : java.util.Optional.of(new long[]{newStart, newEnd});
        }
    }

    private AlignmentTrimmer() {
    }

    /**
     * Restricts {@code alignment} to columns where no sequence has a gap. Sequence order, ids and
     * all metadata are preserved; only the sequence strings shorten.
     *
     * @param minRetainedColumns refuse to trim below this many columns -- a 40-column "alignment"
     *                           is not something to compute population genetics on, and silently
     *                           producing one would trade a visible failure for an invisible one
     */
    public static Result trimToFullyCovered(SequenceAlignment alignment, int minRetainedColumns) {
        List<NucleotideSequence> seqs = alignment.getSequences();
        if (seqs.isEmpty()) {
            return new Result(alignment, 0, 0, 0, false, "Alignment is empty; nothing to trim.", new int[0]);
        }

        int length = alignment.length();
        List<Integer> keep = new ArrayList<>();
        int longestBlock = 0;
        int currentBlock = 0;
        for (int col = 0; col < length; col++) {
            boolean covered = true;
            for (NucleotideSequence s : seqs) {
                String str = s.getSequence();
                if (col >= str.length() || IupacUtil.isGap(str.charAt(col))) {
                    covered = false;
                    break;
                }
            }
            if (covered) {
                keep.add(col);
                currentBlock++;
                longestBlock = Math.max(longestBlock, currentBlock);
            } else {
                currentBlock = 0;
            }
        }

        if (keep.size() == length) {
            int[] all = new int[length];
            for (int i = 0; i < length; i++) all[i] = i;
            return new Result(alignment, length, length, length, false,
                    "No trimming needed: every column is covered by every sequence.", all);
        }
        if (keep.size() < minRetainedColumns) {
            return new Result(alignment, length, keep.size(), longestBlock, false, String.format(
                    "Trimming NOT applied: only %d of %d columns (%.0f%%) are covered by every sequence, below the "
                            + "%d-column minimum. The sequences overlap too little to share a usable window -- they are "
                            + "probably different regions or different genes rather than one locus sequenced to differing "
                            + "lengths. Check what was actually pooled into this alignment.",
                    keep.size(), length, 100.0 * keep.size() / length, minRetainedColumns), new int[0]);
        }

        List<NucleotideSequence> trimmed = new ArrayList<>(seqs.size());
        for (NucleotideSequence s : seqs) {
            String str = s.getSequence();
            StringBuilder sb = new StringBuilder(keep.size());
            for (int col : keep) sb.append(str.charAt(col));
            trimmed.add(new NucleotideSequence(s.getId(), sb.toString(),
                    s.getCollectionDate().orElse(null), s.getLocation().orElse(null), s.getHost().orElse(null)));
        }

        String referenceId = alignment.getReference().getId();
        String description = String.format(
                "Alignment trimmed to the %d of %d columns (%.0f%%) covered by all %d sequences; longest contiguous "
                        + "block %d bp. The discarded columns were positions where at least one sequence had no data, so "
                        + "every index now compares observed bases rather than absent coverage. Note this is NOT a "
                        + "re-alignment: the input was correctly aligned, it simply pooled sequences with different "
                        + "coverage windows.",
                keep.size(), length, 100.0 * keep.size() / length, seqs.size(), longestBlock);

        int[] keptArray = new int[keep.size()];
        for (int i = 0; i < keep.size(); i++) keptArray[i] = keep.get(i);
        return new Result(SequenceAlignment.of(trimmed, referenceId), length, keep.size(), longestBlock, true,
                description, keptArray);
    }
}
