package org.gvi.cli;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVRecord;
import org.gvi.composite.IndexKey;
import org.gvi.composite.WeightCalibrator;
import org.gvi.core.exception.GviInputException;

import java.io.IOException;
import java.io.Reader;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.EnumMap;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Reads historical calibration data for {@link WeightCalibrator} (Section
 * 8.5): a CSV whose header is a subset of the index labels (mu, Re, pi, MB,
 * dN/dS, GD, CAI, GC_Deviation, RI -- already 0..1-normalized values, one
 * row per historical sample/timepoint) plus a required "target" column
 * (the independently observed outcome to fit against, e.g. a measured
 * growth rate).
 */
public final class CalibrationCsvReader {

    public record Result(List<WeightCalibrator.Observation> observations, List<IndexKey> keysInPlay) {
    }

    private CalibrationCsvReader() {
    }

    public static Result read(Path path) {
        if (!Files.isReadable(path)) {
            throw new GviInputException("Cannot read calibration CSV: " + path);
        }
        try (Reader r = Files.newBufferedReader(path, StandardCharsets.UTF_8)) {
            return read(r, path.toString());
        } catch (IOException e) {
            throw new GviInputException("I/O error reading calibration CSV: " + path + " (" + e.getMessage() + ")", e);
        }
    }

    public static Result read(Reader reader, String sourceLabel) {
        CSVFormat format = CSVFormat.DEFAULT.builder().setHeader().setSkipHeaderRecord(true).setTrim(true).build();
        List<WeightCalibrator.Observation> observations = new ArrayList<>();
        Set<IndexKey> keysInPlay = new LinkedHashSet<>();

        try (CSVParser parser = format.parse(reader)) {
            // Index columns are matched case-insensitively (IndexKey.fromLabel), so the target
            // column is too: a file headed TARGET otherwise reported the column it plainly has as
            // missing. The name as written is kept, because record lookup is case-sensitive.
            String targetColumn = parser.getHeaderNames().stream()
                    .filter(h -> h.equalsIgnoreCase("target"))
                    .findFirst()
                    .orElseThrow(() -> new GviInputException(
                            "Calibration CSV (" + sourceLabel + ") is missing required column 'target'"));
            // Key -> the header exactly as written. IndexKey.fromLabel matches case-insensitively
            // but CSVRecord.get does not, so a file headed "MU" resolved to a key and then failed
            // the row lookup for "mu"; carrying the written spelling is what makes the two agree.
            Map<IndexKey, String> columnKeys = new LinkedHashMap<>();
            for (String header : parser.getHeaderNames()) {
                if (header.equals(targetColumn)) continue;
                try {
                    columnKeys.put(IndexKey.fromLabel(header), header);
                } catch (IllegalArgumentException e) {
                    throw new GviInputException("Calibration CSV (" + sourceLabel + ") has unrecognized column '" + header + "': " + e.getMessage());
                }
            }
            if (columnKeys.isEmpty()) {
                throw new GviInputException("Calibration CSV (" + sourceLabel + ") has no index columns besides 'target'");
            }
            keysInPlay.addAll(columnKeys.keySet());

            for (CSVRecord rec : parser) {
                Map<IndexKey, Double> values = new EnumMap<>(IndexKey.class);
                for (Map.Entry<IndexKey, String> column : columnKeys.entrySet()) {
                    values.put(column.getKey(), Double.parseDouble(rec.get(column.getValue())));
                }
                double target = Double.parseDouble(rec.get(targetColumn));
                observations.add(new WeightCalibrator.Observation(values, target));
            }
        } catch (IOException e) {
            throw new GviInputException("I/O error parsing calibration CSV (" + sourceLabel + "): " + e.getMessage(), e);
        }

        if (observations.isEmpty()) {
            throw new GviInputException("Calibration CSV (" + sourceLabel + ") contained no data rows");
        }
        return new Result(observations, new ArrayList<>(keysInPlay));
    }
}
