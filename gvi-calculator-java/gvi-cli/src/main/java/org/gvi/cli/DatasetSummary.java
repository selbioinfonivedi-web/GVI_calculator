package org.gvi.cli;

import java.time.LocalDate;

/**
 * Plain facts about the loaded alignment -- sequence count, alignment length, the reference actually
 * used (resolved from {@code --reference-id}, or the first sequence if that wasn't supplied), and the
 * collection-date span across whichever sequences had one (from metadata or FASTA headers). No new
 * computation -- every field here was already known inside {@link GviPipeline#run} and simply wasn't
 * surfaced to callers before.
 */
public record DatasetSummary(
        int sequenceCount,
        int alignmentLengthBp,
        String referenceId,
        LocalDate earliestCollectionDate,
        LocalDate latestCollectionDate
) {
}
