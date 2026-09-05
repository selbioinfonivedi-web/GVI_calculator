package org.gvi.core.io;

import org.gvi.core.exception.GviInputException;
import org.gvi.core.model.NucleotideSequence;
import org.gvi.core.util.FastaHeaderDateParser;
import org.gvi.core.util.FastaHeaderLocationHostParser;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Reader;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.time.LocalDate;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;

/**
 * Streaming FASTA / multi-FASTA parser. Never loads the whole file as one
 * String -- reads line by line so multi-gigabyte alignments don't blow the
 * heap (Section 2). Malformed individual records (blank sequence, duplicate
 * id) are skipped with a {@link ParseWarning} rather than aborting the
 * whole file.
 */
public final class FastaReader {

    private static final Pattern VALID_CHARS = Pattern.compile("[ACGTURYWSKMBDHVN\\-\\.\\?]*", Pattern.CASE_INSENSITIVE);

    public record Result(List<NucleotideSequence> sequences, ParseReport report) {
    }

    private FastaReader() {
    }

    public static Result read(Path path) {
        if (!Files.isReadable(path)) {
            throw new GviInputException("Cannot read FASTA file: " + path);
        }
        try (BufferedReader br = Files.newBufferedReader(path, StandardCharsets.UTF_8)) {
            return read(br, path.toString());
        } catch (IOException e) {
            throw new GviInputException("I/O error reading FASTA file: " + path + " (" + e.getMessage() + ")", e);
        }
    }

    public static Result read(Reader reader, String sourceLabel) {
        ParseReport report = new ParseReport();
        Map<String, NucleotideSequence> byId = new LinkedHashMap<>();
        String currentId = null;
        String currentHeader = null;
        StringBuilder currentSeq = new StringBuilder();
        long lineNo = 0;
        long currentIdLine = 0;

        try (BufferedReader br = (reader instanceof BufferedReader b) ? b : new BufferedReader(reader)) {
            String line;
            while ((line = br.readLine()) != null) {
                lineNo++;
                if (line.isBlank()) continue;
                if (line.charAt(0) == '>') {
                    flush(currentId, currentHeader, currentSeq, currentIdLine, sourceLabel, byId, report);
                    currentHeader = line.substring(1).trim();
                    currentId = firstToken(currentHeader);
                    currentIdLine = lineNo;
                    currentSeq = new StringBuilder();
                } else {
                    if (currentId == null) {
                        report.skip(new ParseWarning(sourceLabel, lineNo, null,
                                "Sequence data before any '>' header line; ignored"));
                        continue;
                    }
                    currentSeq.append(line.trim());
                }
            }
            flush(currentId, currentHeader, currentSeq, currentIdLine, sourceLabel, byId, report);
        } catch (IOException e) {
            throw new GviInputException("I/O error reading FASTA (" + sourceLabel + "): " + e.getMessage(), e);
        }

        if (byId.isEmpty()) {
            throw new GviInputException("No valid sequences found in FASTA source: " + sourceLabel);
        }
        return new Result(new ArrayList<>(byId.values()), report);
    }

    private static void flush(String id, String header, StringBuilder seq, long idLine, String sourceLabel,
                               Map<String, NucleotideSequence> byId, ParseReport report) {
        if (id == null) return;
        if (seq.isEmpty()) {
            report.skip(new ParseWarning(sourceLabel, idLine, id, "Empty sequence; skipped"));
            return;
        }
        if (byId.containsKey(id)) {
            report.skip(new ParseWarning(sourceLabel, idLine, id, "Duplicate sequence id; first occurrence kept, this one skipped"));
            return;
        }
        String seqStr = seq.toString();
        if (!VALID_CHARS.matcher(seqStr).matches()) {
            report.skip(new ParseWarning(sourceLabel, idLine, id,
                    "Contains characters outside the IUPAC nucleotide alphabet; skipped (check for protein sequence or corrupted file)"));
            return;
        }

        // Automatic temporal fallback: if the header itself encodes a collection date (a very common
        // real-world convention -- GISAID/Nextstrain-style "id|2021-03-15|..."), use it so mu/Re work
        // without requiring a separate --metadata CSV. An explicit metadata CSV, if supplied, overrides
        // this later in MetadataJoiner -- deliberate user input wins over a header-parsed guess.
        NucleotideSequence sequence = new NucleotideSequence(id, seqStr);
        var headerDate = FastaHeaderDateParser.extract(header);
        LocalDate date = headerDate.map(FastaHeaderDateParser.Result::date).orElse(null);
        if (headerDate.isPresent() && headerDate.get().yearOnly()) {
            report.note(new ParseWarning(sourceLabel, idLine, id,
                    "Header gave only a year, not a full date; using " + date
                            + " (the middle of that year) as a reduced-precision estimate for temporal indices (mu, Re)"));
        }

        // Same idea, for location/host: fires only for the specific "accession|location|host|date"
        // convention (see FastaHeaderLocationHostParser) -- a real, common convention, not a guess at
        // an arbitrary header's field meanings. Also overridden later by an explicit metadata CSV.
        var headerLocationHost = FastaHeaderLocationHostParser.extract(header);
        String location = headerLocationHost.map(FastaHeaderLocationHostParser.Result::location).orElse(null);
        String host = headerLocationHost.map(FastaHeaderLocationHostParser.Result::host).orElse(null);

        if (date != null || location != null || host != null) {
            sequence = sequence.withMetadata(date, location, host);
        }

        byId.put(id, sequence);
        report.accept();
    }

    private static String firstToken(String header) {
        int sp = header.indexOf(' ');
        int tab = header.indexOf('\t');
        int cut = header.length();
        if (sp >= 0) cut = Math.min(cut, sp);
        if (tab >= 0) cut = Math.min(cut, tab);
        String id = header.substring(0, cut).trim();
        if (id.isEmpty()) {
            throw new GviInputException("FASTA header line has no usable id: '>" + header + "'");
        }
        return id;
    }
}
