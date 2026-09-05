package org.gvi.core.io;

import org.gvi.core.exception.GviInputException;
import org.gvi.core.model.VariantRecord;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.GZIPInputStream;

/**
 * Minimal, dependency-free VCF parser (v4.x) covering exactly what Mutation
 * Burden (Section 5.4) needs: chrom/pos/ref/alt and per-site read depth
 * (DP) for the coverage filter. Deliberately not a full VCF library --
 * malformed data lines are skipped with a warning rather than aborting.
 */
public final class VcfReader {

    public record Result(List<VariantRecord> variants, ParseReport report) {
    }

    private VcfReader() {
    }

    public static Result read(Path path) {
        if (!Files.isReadable(path)) {
            throw new GviInputException("Cannot read VCF file: " + path);
        }
        try (InputStream raw = Files.newInputStream(path);
             InputStream in = path.toString().endsWith(".gz") ? new GZIPInputStream(raw) : raw;
             BufferedReader br = new BufferedReader(new java.io.InputStreamReader(in, StandardCharsets.UTF_8))) {
            return read(br, path.toString());
        } catch (IOException e) {
            throw new GviInputException("I/O error reading VCF file: " + path + " (" + e.getMessage() + ")", e);
        }
    }

    public static Result read(BufferedReader br, String sourceLabel) throws IOException {
        ParseReport report = new ParseReport();
        List<VariantRecord> variants = new ArrayList<>();
        String line;
        long lineNo = 0;
        List<String> formatFields = null;
        boolean sawHeader = false;

        while ((line = br.readLine()) != null) {
            lineNo++;
            if (line.isBlank()) continue;
            if (line.startsWith("##")) continue;
            if (line.startsWith("#CHROM")) {
                sawHeader = true;
                continue;
            }
            if (!sawHeader) {
                report.skip(new ParseWarning(sourceLabel, lineNo, null, "Data line before #CHROM header; skipped"));
                continue;
            }
            String[] cols = line.split("\t", -1);
            if (cols.length < 8) {
                report.skip(new ParseWarning(sourceLabel, lineNo, null, "Fewer than 8 mandatory VCF columns; skipped"));
                continue;
            }
            String chrom = cols[0];
            String posStr = cols[1];
            String ref = cols[3];
            String altField = cols[4];
            String info = cols[7];
            long pos;
            try {
                pos = Long.parseLong(posStr);
            } catch (NumberFormatException e) {
                report.skip(new ParseWarning(sourceLabel, lineNo, chrom + ":" + posStr, "Non-numeric POS; skipped"));
                continue;
            }
            if (altField.equals(".") || altField.isBlank()) {
                report.skip(new ParseWarning(sourceLabel, lineNo, chrom + ":" + pos, "No ALT allele (monomorphic site); skipped"));
                continue;
            }

            Integer depth = extractInfoDp(info);
            if (depth == null && cols.length >= 10) {
                formatFields = List.of(cols[8].split(":"));
                depth = extractSampleDp(formatFields, cols[9]);
            }

            for (String alt : altField.split(",")) {
                if (alt.equals("*") || alt.isBlank()) continue; // spanning-deletion placeholder, not a real call
                variants.add(new VariantRecord(chrom, pos, ref, alt, depth, cols.length >= 10 ? "SAMPLE1" : null));
                report.accept();
            }
        }
        return new Result(variants, report);
    }

    private static Integer extractInfoDp(String info) {
        for (String kv : info.split(";")) {
            if (kv.startsWith("DP=")) {
                try {
                    return Integer.parseInt(kv.substring(3).trim());
                } catch (NumberFormatException e) {
                    return null;
                }
            }
        }
        return null;
    }

    private static Integer extractSampleDp(List<String> formatFields, String sampleColumn) {
        int dpIdx = formatFields.indexOf("DP");
        if (dpIdx < 0) return null;
        String[] values = sampleColumn.split(":");
        if (dpIdx >= values.length) return null;
        try {
            return Integer.parseInt(values[dpIdx]);
        } catch (NumberFormatException e) {
            return null;
        }
    }
}
