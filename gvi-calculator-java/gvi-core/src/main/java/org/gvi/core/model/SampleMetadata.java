package org.gvi.core.model;

import java.time.LocalDate;

/** One row of the temporal/contextual metadata CSV (Section 4.2), keyed by sequence_id. */
public record SampleMetadata(String sequenceId, LocalDate collectionDate, String location, String host) {
}
