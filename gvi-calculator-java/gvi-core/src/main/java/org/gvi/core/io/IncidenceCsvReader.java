package org.gvi.core.io;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVRecord;
import org.gvi.core.exception.GviInputException;
import org.gvi.core.model.IncidencePoint;

import java.io.IOException;
import java.io.Reader;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.time.LocalDate;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

/** Parses the optional case-incidence time series CSV (date, new_cases) used by the Cori et al. Re estimator. */
public final class IncidenceCsvReader {

    public record Result(List<IncidencePoint> points, ParseReport report) {
    }

    private IncidenceCsvReader() {
    }

    public static Result read(Path path) {
        if (!Files.isReadable(path)) {
            throw new GviInputException("Cannot read incidence CSV: " + path);
        }
        try (Reader r = Files.newBufferedReader(path, StandardCharsets.UTF_8)) {
            return read(r, path.toString());
        } catch (IOException e) {
            throw new GviInputException("I/O error reading incidence CSV: " + path + " (" + e.getMessage() + ")", e);
        }
    }

    public static Result read(Reader reader, String sourceLabel) {
        ParseReport report = new ParseReport();
        List<IncidencePoint> points = new ArrayList<>();
        CSVFormat format = CSVFormat.DEFAULT.builder().setHeader().setSkipHeaderRecord(true).setTrim(true).build();

        try (CSVParser parser = format.parse(reader)) {
            if (!parser.getHeaderMap().containsKey("date") || !parser.getHeaderMap().containsKey("new_cases")) {
                throw new GviInputException("Incidence CSV (" + sourceLabel + ") must have columns 'date' and 'new_cases'");
            }
            for (CSVRecord rec : parser) {
                long line = rec.getRecordNumber() + 1;
                try {
                    LocalDate date = LocalDate.parse(rec.get("date").trim());
                    double cases = Double.parseDouble(rec.get("new_cases").trim());
                    if (cases < 0) {
                        report.skip(new ParseWarning(sourceLabel, line, date.toString(), "Negative new_cases; skipped"));
                        continue;
                    }
                    points.add(new IncidencePoint(date, cases));
                    report.accept();
                } catch (Exception e) {
                    report.skip(new ParseWarning(sourceLabel, line, null, "Unparsable row (" + e.getMessage() + "); skipped"));
                }
            }
        } catch (IOException e) {
            throw new GviInputException("I/O error parsing incidence CSV (" + sourceLabel + "): " + e.getMessage(), e);
        }

        if (points.isEmpty()) {
            throw new GviInputException("Incidence CSV (" + sourceLabel + ") contained no usable rows");
        }
        points.sort(Comparator.comparing(IncidencePoint::date));
        return new Result(points, report);
    }
}
