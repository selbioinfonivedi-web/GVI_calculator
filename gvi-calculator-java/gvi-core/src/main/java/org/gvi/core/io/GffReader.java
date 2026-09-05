package org.gvi.core.io;

import org.gvi.core.exception.GviInputException;
import org.gvi.core.model.GeneAnnotation;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

/**
 * Minimal GFF3 reader that extracts CDS regions (used to pull in-frame
 * codons for dN/dS and CAI, Sections 5.5/5.8). Only the 'CDS' feature type
 * is used; other feature rows are ignored, not errors.
 */
public final class GffReader {

    public record Result(List<GeneAnnotation> genes, ParseReport report) {
    }

    private GffReader() {
    }

    public static Result read(Path path) {
        if (!Files.isReadable(path)) {
            throw new GviInputException("Cannot read GFF3 file: " + path);
        }
        try (BufferedReader br = Files.newBufferedReader(path, StandardCharsets.UTF_8)) {
            return read(br, path.toString());
        } catch (IOException e) {
            throw new GviInputException("I/O error reading GFF3 file: " + path + " (" + e.getMessage() + ")", e);
        }
    }

    public static Result read(BufferedReader br, String sourceLabel) throws IOException {
        ParseReport report = new ParseReport();
        List<GeneAnnotation> genes = new ArrayList<>();
        String line;
        long lineNo = 0;
        while ((line = br.readLine()) != null) {
            lineNo++;
            if (line.isBlank() || line.startsWith("#")) continue;
            String[] cols = line.split("\t", -1);
            if (cols.length < 9) {
                report.skip(new ParseWarning(sourceLabel, lineNo, null, "Fewer than 9 GFF3 columns; skipped"));
                continue;
            }
            if (!"CDS".equalsIgnoreCase(cols[2])) continue;
            try {
                long start = Long.parseLong(cols[3]);
                long end = Long.parseLong(cols[4]);
                char strand = cols[6].isBlank() ? '+' : cols[6].charAt(0);
                if (strand != '+' && strand != '-') strand = '+';
                String name = extractName(cols[8]);
                genes.add(new GeneAnnotation(name, start, end, strand));
                report.accept();
            } catch (Exception e) {
                report.skip(new ParseWarning(sourceLabel, lineNo, null, "Malformed CDS record (" + e.getMessage() + "); skipped"));
            }
        }
        if (genes.isEmpty()) {
            report.skip(new ParseWarning(sourceLabel, 0, null, "No CDS features found; gene-specific dN/dS and CAI will use the whole sequence as one ORF"));
        }
        return new Result(genes, report);
    }

    private static String extractName(String attributes) {
        for (String key : new String[]{"gene=", "Name=", "ID="}) {
            for (String part : attributes.split(";")) {
                if (part.startsWith(key)) {
                    return part.substring(key.length()).trim();
                }
            }
        }
        return "unnamed";
    }
}
