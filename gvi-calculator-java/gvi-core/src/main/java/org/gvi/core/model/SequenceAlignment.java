package org.gvi.core.model;

import org.gvi.core.exception.GviInputException;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * A set of aligned nucleotide sequences (equal length, gaps included) plus
 * a designated reference. Sequences of unequal length are rejected up front
 * with a clear message rather than causing an ArrayIndexOutOfBounds deep
 * inside a pairwise-distance loop later.
 */
public final class SequenceAlignment {

    private final List<NucleotideSequence> sequences;
    private final int referenceIndex;
    private final int alignmentLength;

    private SequenceAlignment(List<NucleotideSequence> sequences, int referenceIndex) {
        this.sequences = List.copyOf(sequences);
        this.referenceIndex = referenceIndex;
        this.alignmentLength = sequences.get(0).length();
    }

    /**
     * Builds an alignment, validating equal length across all sequences.
     * The first sequence is treated as the reference unless
     * {@code referenceId} is supplied.
     */
    public static SequenceAlignment of(List<NucleotideSequence> sequences, String referenceId) {
        if (sequences == null || sequences.isEmpty()) {
            throw new GviInputException("Alignment must contain at least one sequence");
        }
        int expectedLength = sequences.get(0).length();
        List<NucleotideSequence> mismatched = new ArrayList<>();
        for (NucleotideSequence seq : sequences) {
            if (seq.length() != expectedLength) {
                mismatched.add(seq);
            }
        }
        if (!mismatched.isEmpty()) {
            String ids = mismatched.stream().map(NucleotideSequence::getId).limit(10)
                    .reduce((a, b) -> a + ", " + b).orElse("");
            throw new GviInputException(
                    "Alignment length mismatch: expected " + expectedLength + " bp (from '"
                            + sequences.get(0).getId() + "'), but these sequences differ: " + ids
                            + (mismatched.size() > 10 ? " (+" + (mismatched.size() - 10) + " more)" : "")
                            + ". This calculator requires a pre-aligned multi-FASTA; align your sequences first.");
        }
        int refIdx = 0;
        if (referenceId != null && !referenceId.isBlank()) {
            refIdx = -1;
            for (int i = 0; i < sequences.size(); i++) {
                if (sequences.get(i).getId().equals(referenceId)) {
                    refIdx = i;
                    break;
                }
            }
            if (refIdx == -1) {
                throw new GviInputException("Reference id '" + referenceId + "' not found in alignment");
            }
        }
        return new SequenceAlignment(sequences, refIdx);
    }

    public static SequenceAlignment of(List<NucleotideSequence> sequences) {
        return of(sequences, null);
    }

    public List<NucleotideSequence> getSequences() {
        return sequences;
    }

    public NucleotideSequence getReference() {
        return sequences.get(referenceIndex);
    }

    public List<NucleotideSequence> getQueries() {
        List<NucleotideSequence> queries = new ArrayList<>(sequences);
        queries.remove(referenceIndex);
        return Collections.unmodifiableList(queries);
    }

    public int size() {
        return sequences.size();
    }

    public int length() {
        return alignmentLength;
    }
}
