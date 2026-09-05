package org.gvi.core.util;

import java.util.Optional;

/**
 * Opportunistically extracts location and host from a FASTA header when it
 * follows the common {@code accession|location|host|date} pipe-delimited
 * convention (e.g. {@code JF416958.1|India|HUMAN|1957},
 * {@code KP821387.1|France|domestic_sheep|2001}) -- a real, widely-used
 * convention, including this project's own KFDV test data.
 * <p>
 * Deliberately conservative, same discipline as {@link FastaHeaderDateParser}:
 * unlike that class (which scans the whole header for a date-shaped
 * substring anywhere), this only fires when the header splits into EXACTLY
 * 4 pipe-delimited fields AND the last field is itself independently
 * recognized as a date/year by {@link FastaHeaderDateParser} -- confirming
 * this really is the accession|location|host|date convention, not some
 * other 4-field header shape that happens to also use pipes. Headers that
 * don't match this specific shape (GISAID-style {@code id|date}, Nextstrain-style
 * {@code id|date|country}, or anything else) return empty rather than
 * guessing which token might be a location vs. a host.
 */
public final class FastaHeaderLocationHostParser {

    public record Result(String location, String host) {
    }

    private FastaHeaderLocationHostParser() {
    }

    public static Optional<Result> extract(String header) {
        if (header == null) return Optional.empty();
        String[] fields = header.split("\\|", -1);
        if (fields.length != 4) return Optional.empty();

        String accession = fields[0].trim();
        String location = fields[1].trim();
        String host = fields[2].trim();
        String dateField = fields[3].trim();

        if (accession.isEmpty() || location.isEmpty() || host.isEmpty()) return Optional.empty();
        if (FastaHeaderDateParser.extract(dateField).isEmpty()) return Optional.empty();

        return Optional.of(new Result(location, host));
    }
}
