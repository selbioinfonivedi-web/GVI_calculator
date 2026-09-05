package org.gvi.core.util;

import java.time.DateTimeException;
import java.time.LocalDate;
import java.util.Optional;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Extracts a collection date directly from a FASTA header, so datasets that
 * already encode dates in their headers (a very common real-world
 * convention -- GISAID EpiCoV-style {@code hCoV-19/USA/CA-1/2020|2020-03-15},
 * Nextstrain-style {@code sample|2021-06-01|USA}, or simply
 * {@code sample_2021-06-01}) don't require a separate metadata CSV just to
 * unlock the temporal indices (mu, Re).
 * <p>
 * Strategy: scan the whole header for the first ISO-8601 {@code YYYY-MM-DD}
 * substring; if none is found, fall back to a standalone 4-digit year
 * (flagged as reduced precision, resolved to the middle of that year --
 * see {@link TemporalUtil#midpointOfYear} for why January 1st would bias
 * every downstream regression toward the start of the year). This is a
 * best-effort heuristic, not a full date-format parser -- ambiguous or
 * malformed dates are simply not matched (returns empty) rather than
 * guessed at.
 */
public final class FastaHeaderDateParser {

    private static final Pattern ISO_DATE = Pattern.compile("(\\d{4})-(\\d{2})-(\\d{2})");
    // Note: NOT \b...\b -- underscore is a "word" character in regex, so \b would fail to match
    // years in the very common "sample_2021_..." header convention (no boundary between '_' and '2').
    // Digit-adjacency lookaround avoids matching a fragment of a longer number instead.
    private static final Pattern YEAR_ONLY = Pattern.compile("(?<!\\d)(19|20)\\d{2}(?!\\d)");

    public record Result(LocalDate date, boolean yearOnly) {
    }

    private FastaHeaderDateParser() {
    }

    public static Optional<Result> extract(String header) {
        if (header == null) return Optional.empty();

        Matcher iso = ISO_DATE.matcher(header);
        while (iso.find()) {
            try {
                LocalDate date = LocalDate.of(
                        Integer.parseInt(iso.group(1)), Integer.parseInt(iso.group(2)), Integer.parseInt(iso.group(3)));
                return Optional.of(new Result(date, false));
            } catch (DateTimeException | NumberFormatException e) {
                // not a real calendar date (e.g. a version-like "2024-99-99"); keep scanning for another match
            }
        }

        Matcher year = YEAR_ONLY.matcher(header);
        if (year.find()) {
            int y = Integer.parseInt(year.group());
            return Optional.of(new Result(TemporalUtil.midpointOfYear(y), true));
        }

        return Optional.empty();
    }
}
