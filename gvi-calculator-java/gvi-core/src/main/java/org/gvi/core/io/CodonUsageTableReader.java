package org.gvi.core.io;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVRecord;
import org.gvi.core.exception.GviInputException;
import org.gvi.core.model.CodonUsageTable;

import java.io.IOException;
import java.io.Reader;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.Map;

/** Parses a host codon-usage reference table CSV (columns: codon, frequency) for CAI (Section 5.8). */
public final class CodonUsageTableReader {

    private CodonUsageTableReader() {
    }

    public static CodonUsageTable read(Path path, String tableName) {
        if (!Files.isReadable(path)) {
            throw new GviInputException("Cannot read codon usage table: " + path);
        }
        try (Reader r = Files.newBufferedReader(path, StandardCharsets.UTF_8)) {
            return read(r, tableName);
        } catch (IOException e) {
            throw new GviInputException("I/O error reading codon usage table: " + path + " (" + e.getMessage() + ")", e);
        }
    }

    public static CodonUsageTable read(Reader reader, String tableName) {
        Map<String, Double> freqs = new HashMap<>();
        CSVFormat format = CSVFormat.DEFAULT.builder().setHeader().setSkipHeaderRecord(true).setTrim(true).build();
        try (CSVParser parser = format.parse(reader)) {
            for (CSVRecord rec : parser) {
                String codon = rec.get("codon").toUpperCase().replace('U', 'T');
                double freq = Double.parseDouble(rec.get("frequency"));
                freqs.put(codon, freq);
            }
        } catch (Exception e) {
            throw new GviInputException("Malformed codon usage table '" + tableName + "': " + e.getMessage(), e);
        }
        return CodonUsageTable.fromFrequencies(tableName, freqs);
    }
}
