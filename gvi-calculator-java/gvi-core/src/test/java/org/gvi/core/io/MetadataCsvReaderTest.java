package org.gvi.core.io;

import org.gvi.core.exception.GviInputException;
import org.gvi.core.model.SampleMetadata;
import org.junit.jupiter.api.Test;

import java.io.StringReader;
import java.time.LocalDate;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.assertThatThrownBy;

class MetadataCsvReaderTest {

    @Test
    void parsesValidRows() {
        String csv = "sequence_id,collection_date,location,host\nseq1,2024-01-15,India,Human\nseq2,2024-02-01,,\n";
        MetadataCsvReader.Result result = MetadataCsvReader.read(new StringReader(csv), "test");
        assertThat(result.rows()).hasSize(2);
        SampleMetadata r1 = result.rows().get(0);
        assertThat(r1.sequenceId()).isEqualTo("seq1");
        assertThat(r1.collectionDate()).isEqualTo(LocalDate.of(2024, 1, 15));
        assertThat(r1.location()).isEqualTo("India");
    }

    @Test
    void skipsRowWithUnparsableDateInsteadOfCrashing() {
        String csv = "sequence_id,collection_date\nseq1,not-a-date\nseq2,2024-02-01\n";
        MetadataCsvReader.Result result = MetadataCsvReader.read(new StringReader(csv), "test");
        // row is kept (id still usable) but date is null, with a warning recorded
        assertThat(result.rows()).hasSize(2);
        assertThat(result.rows().get(0).collectionDate()).isNull();
        assertThat(result.report().hasWarnings()).isTrue();
    }

    @Test
    void acceptsYearOnlyDateAtReducedPrecisionInsteadOfSkipping() {
        // Bacterial/parasite GenBank deposits routinely carry only a year -- previously rejected outright
        // (same treatment as "not-a-date"), losing the row's temporal signal entirely. Now resolved to the
        // midpoint of that year (an unbiased single-date estimate) with a note, not a skip.
        String csv = "sequence_id,collection_date\nseq1,2019\n";
        MetadataCsvReader.Result result = MetadataCsvReader.read(new StringReader(csv), "test");
        assertThat(result.rows()).hasSize(1);
        assertThat(result.rows().get(0).collectionDate())
                .isEqualTo(org.gvi.core.util.TemporalUtil.midpointOfYear(2019));
        assertThat(result.report().getSkippedCount()).isZero();
        assertThat(result.report().hasWarnings()).isTrue();
    }

    @Test
    void acceptsYearMonthDateAtReducedPrecisionInsteadOfSkipping() {
        String csv = "sequence_id,collection_date\nseq1,2019-06\n";
        MetadataCsvReader.Result result = MetadataCsvReader.read(new StringReader(csv), "test");
        assertThat(result.rows()).hasSize(1);
        assertThat(result.rows().get(0).collectionDate())
                .isEqualTo(org.gvi.core.util.TemporalUtil.midpointOfYearMonth(2019, 6));
        assertThat(result.report().getSkippedCount()).isZero();
        assertThat(result.report().hasWarnings()).isTrue();
    }

    @Test
    void requiresSequenceIdColumn() {
        String csv = "foo,bar\n1,2\n";
        assertThatThrownBy(() -> MetadataCsvReader.read(new StringReader(csv), "test"))
                .isInstanceOf(GviInputException.class)
                .hasMessageContaining("sequence_id");
    }
}
