package org.gvi.core.model;

import java.time.LocalDate;
import java.util.Optional;

/**
 * A single nucleotide sequence record: FASTA content plus whatever temporal
 * / contextual metadata was joined in from the metadata CSV (Section 4.2 of
 * the build spec). Sequence characters are stored upper-cased; IUPAC
 * ambiguity codes and '-' gaps are preserved as-is (not stripped), since
 * individual algorithms decide how to treat them.
 */
public final class NucleotideSequence {

    private final String id;
    private final String sequence;
    private final LocalDate collectionDate;
    private final String location;
    private final String host;

    public NucleotideSequence(String id, String sequence, LocalDate collectionDate, String location, String host) {
        if (id == null || id.isBlank()) {
            throw new IllegalArgumentException("sequence id must not be blank");
        }
        if (sequence == null) {
            throw new IllegalArgumentException("sequence must not be null");
        }
        this.id = id;
        this.sequence = sequence.toUpperCase();
        this.collectionDate = collectionDate;
        this.location = location;
        this.host = host;
    }

    public NucleotideSequence(String id, String sequence) {
        this(id, sequence, null, null, null);
    }

    public NucleotideSequence withMetadata(LocalDate collectionDate, String location, String host) {
        return new NucleotideSequence(id, sequence, collectionDate, location, host);
    }

    public String getId() {
        return id;
    }

    public String getSequence() {
        return sequence;
    }

    public int length() {
        return sequence.length();
    }

    public Optional<LocalDate> getCollectionDate() {
        return Optional.ofNullable(collectionDate);
    }

    public Optional<String> getLocation() {
        return Optional.ofNullable(location);
    }

    public Optional<String> getHost() {
        return Optional.ofNullable(host);
    }

    @Override
    public String toString() {
        return "NucleotideSequence{" + id + ", len=" + sequence.length() + "}";
    }
}
