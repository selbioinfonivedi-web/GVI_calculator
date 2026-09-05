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
            if (!parser.getHeaderMap().containsKey("target")) {
                throw new GviInputException("Calibration CSV (" + sourceLabel + ") is missing required column 'target'");
            }
            List<IndexKey> columnKeys = new ArrayList<>();
            for (String header : parser.getHeaderNames()) {
                if (header.equalsIgnoreCase("target")) continue;
                try {
                    columnKeys.add(IndexKey.fromLabel(header));
                } catch (IllegalArgumentException e) {
                    throw new GviInputException("Calibration CSV (" + sourceLabel + ") has unrecognized column '" + header + "': " + e.getMessage());
                }
            }
            if (columnKeys.isEmpty()) {
                throw new GviInputException("Calibration CSV (" + sourceLabel + ") has no index columns besides 'target'");
            }
            keysInPlay.addAll(columnKeys);

            for (CSVRecord rec : parser) {
                Map<IndexKey, Double> values = new EnumMap<>(IndexKey.class);
                for (IndexKey key : columnKeys) {
                    values.put(key, Double.parseDouble(rec.get(key.label())));
                }
                double target = Double.parseDouble(rec.get("target"));
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
