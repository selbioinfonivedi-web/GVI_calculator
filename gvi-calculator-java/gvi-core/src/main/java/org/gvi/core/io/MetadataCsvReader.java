package org.gvi.core.io;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVRecord;
import org.gvi.core.exception.GviInputException;
import org.gvi.core.model.SampleMetadata;
import org.gvi.core.util.TemporalUtil;

import java.io.IOException;
import java.io.Reader;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.time.DateTimeException;
import java.time.LocalDate;
import java.time.format.DateTimeParseException;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Parses the temporal/contextual metadata CSV (required column:
 * sequence_id, collection_date; optional: location, host). Rows with an
 * unparsable date are skipped with a warning, never crash the whole load.
 * Bacterial and parasite deposits routinely carry only year or year-month
 * precision (unlike day-precision viral surveillance dates) -- {@code YYYY}
 * and {@code YYYY-MM} are accepted too, resolved to the midpoint of that
 * period (see {@link TemporalUtil#midpointOfYear}) with a reduced-precision
 * note rather than being rejected outright.
 */
public final class MetadataCsvReader {

    private static final Pattern YEAR_ONLY = Pattern.compile("(\\d{4})");
    private static final Pattern YEAR_MONTH = Pattern.compile("(\\d{4})-(\\d{2})");

    public record Result(List<SampleMetadata> rows, ParseReport report) {
    }

    private MetadataCsvReader() {
    }

    public static Result read(Path path) {
        if (!Files.isReadable(path)) {
            throw new GviInputException("Cannot read metadata CSV: " + path);
        }
        try (Reader r = Files.newBufferedReader(path, StandardCharsets.UTF_8)) {
            return read(r, path.toString());
        } catch (IOException e) {
            throw new GviInputException("I/O error reading metadata CSV: " + path + " (" + e.getMessage() + ")", e);
        }
    }

    public static Result read(Reader reader, String sourceLabel) {
        ParseReport report = new ParseReport();
        List<SampleMetadata> rows = new ArrayList<>();

        CSVFormat format = CSVFormat.DEFAULT.builder()
                .setHeader().setSkipHeaderRecord(true).setTrim(true).setIgnoreSurroundingSpaces(true).build();

        try (CSVParser parser = format.parse(reader)) {
            if (!parser.getHeaderMap().containsKey("sequence_id")) {
                throw new GviInputException("Metadata CSV (" + sourceLabel + ") is missing required column 'sequence_id'");
            }
            boolean hasDate = parser.getHeaderMap().containsKey("collection_date");
            if (!hasDate) {
                report.skip(new ParseWarning(sourceLabel, 1, null,
                        "No 'collection_date' column present -- indices requiring temporal data (mu, Re) will be unavailable"));
            }
            for (CSVRecord rec : parser) {
                long line = rec.getRecordNumber() + 1;
                String id = safeGet(rec, "sequence_id");
                if (id == null || id.isBlank()) {
                    report.skip(new ParseWarning(sourceLabel, line, null, "Blank sequence_id; row skipped"));
                    continue;
                }
                LocalDate date = null;
                if (hasDate) {
                    String dateStr = safeGet(rec, "collection_date");
                    if (dateStr != null && !dateStr.isBlank()) {
                        date = parseDate(dateStr.trim(), sourceLabel, line, id, report);
                    }
                }
                String location = safeGet(rec, "location");
                String host = safeGet(rec, "host");
                rows.add(new SampleMetadata(id, date, location, host));
                report.accept();
            }
        } catch (IOException e) {
            throw new GviInputException("I/O error parsing metadata CSV (" + sourceLabel + "): " + e.getMessage(), e);
        }

        if (rows.isEmpty()) {
            throw new GviInputException("Metadata CSV (" + sourceLabel + ") contained no usable rows");
        }
        return new Result(rows, report);
    }

    /**
     * Full {@code yyyy-MM-dd} first; falls back to {@code yyyy-MM} then bare {@code yyyy} (both resolved to
     * the midpoint of that period, with a note rather than a skip -- the row is still usable, just at
     * reduced precision). Only a string matching none of the three is treated as truly unparsable.
     */
    private static LocalDate parseDate(String dateStr, String sourceLabel, long line, String id, ParseReport report) {
        try {
            return LocalDate.parse(dateStr);
        } catch (DateTimeParseException fullDateFailure) {
            Matcher yearMonth = YEAR_MONTH.matcher(dateStr);
            if (yearMonth.matches()) {
                try {
                    LocalDate date = TemporalUtil.midpointOfYearMonth(
                            Integer.parseInt(yearMonth.group(1)), Integer.parseInt(yearMonth.group(2)));
                    report.note(new ParseWarning(sourceLabel, line, id,
                            "collection_date '" + dateStr + "' has only year-month precision; using " + date
                                    + " (the middle of that month) as a reduced-precision estimate"));
                    return date;
                } catch (DateTimeException invalidMonth) {
                    // e.g. "2019-13" -- fall through to the unparsable case below
                }
            }
            Matcher yearOnly = YEAR_ONLY.matcher(dateStr);
            if (yearOnly.matches()) {
                LocalDate date = TemporalUtil.midpointOfYear(Integer.parseInt(yearOnly.group(1)));
                report.note(new ParseWarning(sourceLabel, line, id,
                        "collection_date '" + dateStr + "' has only year precision; using " + date
                                + " (the middle of that year) as a reduced-precision estimate"));
                return date;
            }
            report.skip(new ParseWarning(sourceLabel, line, id,
                    "Unparsable collection_date '" + dateStr + "' (expected ISO-8601 yyyy-MM-dd, or yyyy-MM / yyyy "
                            + "for reduced precision); date omitted for this row"));
            return null;
        }
    }

    private static String safeGet(CSVRecord rec, String column) {
        try {
            if (!rec.isMapped(column)) return null;
            String v = rec.get(column);
            return (v == null || v.isBlank()) ? null : v;
        } catch (IllegalArgumentException e) {
            return null;
        }
    }
}
