package org.gvi.core.io;

import org.gvi.core.model.NucleotideSequence;
import org.gvi.core.model.SampleMetadata;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Joins parsed FASTA sequences with the metadata CSV by sequence_id.
 * Sequences with no matching metadata row are kept (metadata is optional
 * for indices that don't need it) but flagged in the report; per Section
 * 4.3 this must never crash the load.
 * <p>
 * A sequence may already carry a collection date, location, and/or host
 * extracted from its FASTA header (see {@link FastaReader}/
 * {@code FastaHeaderDateParser}/{@code FastaHeaderLocationHostParser})
 * before it ever reaches this join. An explicit metadata CSV value, when
 * present, intentionally overrides the header-derived one for that field
 * (deliberate user input beats a header heuristic) -- but a metadata row
 * with a *blank* cell for a field must not blank out a header-derived value
 * that was already there; it should simply leave it alone. This applies
 * uniformly to date, location, and host.
 */
public final class MetadataJoiner {

    private MetadataJoiner() {
    }

    public record Result(List<NucleotideSequence> sequences, ParseReport report) {
    }

    public static Result join(List<NucleotideSequence> sequences, List<SampleMetadata> metadata) {
        ParseReport report = new ParseReport();
        Map<String, SampleMetadata> byId = new HashMap<>();
        for (SampleMetadata m : metadata) {
            byId.put(m.sequenceId(), m);
        }
        List<NucleotideSequence> joined = new ArrayList<>(sequences.size());
        for (NucleotideSequence seq : sequences) {
            SampleMetadata m = byId.get(seq.getId());
            if (m == null) {
                if (seq.getCollectionDate().isEmpty()) {
                    report.skip(new ParseWarning("metadata-join", 0, seq.getId(),
                            "No metadata row found for this sequence id; temporal fields left empty"));
                } else {
                    report.accept(); // already dated from its FASTA header, nothing to add
                }
                joined.add(seq);
            } else {
                var effectiveDate = m.collectionDate() != null ? m.collectionDate() : seq.getCollectionDate().orElse(null);
                var effectiveLocation = m.location() != null ? m.location() : seq.getLocation().orElse(null);
                var effectiveHost = m.host() != null ? m.host() : seq.getHost().orElse(null);
                joined.add(seq.withMetadata(effectiveDate, effectiveLocation, effectiveHost));
                report.accept();
            }
        }
        return new Result(joined, report);
    }
}
