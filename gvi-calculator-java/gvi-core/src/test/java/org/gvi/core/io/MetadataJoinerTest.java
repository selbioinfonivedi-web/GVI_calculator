package org.gvi.core.io;

import org.gvi.core.model.NucleotideSequence;
import org.gvi.core.model.SampleMetadata;
import org.junit.jupiter.api.Test;

import java.time.LocalDate;
import java.util.List;

import static org.assertj.core.api.Assertions.assertThat;

class MetadataJoinerTest {

    @Test
    void explicitMetadataDateOverridesHeaderDerivedDate() {
        NucleotideSequence seq = new NucleotideSequence("s1", "ACGT")
                .withMetadata(LocalDate.of(2020, 1, 1), null, null); // simulates a header-derived date
        SampleMetadata meta = new SampleMetadata("s1", LocalDate.of(2021, 6, 15), "India", "Human");

        var result = MetadataJoiner.join(List.of(seq), List.of(meta));

        assertThat(result.sequences().get(0).getCollectionDate()).contains(LocalDate.of(2021, 6, 15));
    }

    @Test
    void metadataRowWithBlankDateDoesNotEraseHeaderDerivedDate() {
        NucleotideSequence seq = new NucleotideSequence("s1", "ACGT")
                .withMetadata(LocalDate.of(2020, 1, 1), null, null); // simulates a header-derived date
        SampleMetadata meta = new SampleMetadata("s1", null, "India", "Human"); // no date column value for this row

        var result = MetadataJoiner.join(List.of(seq), List.of(meta));

        assertThat(result.sequences().get(0).getCollectionDate()).contains(LocalDate.of(2020, 1, 1));
        assertThat(result.sequences().get(0).getLocation()).contains("India"); // other metadata fields still applied
    }

    @Test
    void sequenceWithNoMatchingMetadataRowKeepsItsHeaderDerivedDate() {
        NucleotideSequence seq = new NucleotideSequence("s1", "ACGT")
                .withMetadata(LocalDate.of(2020, 1, 1), null, null);

        var result = MetadataJoiner.join(List.of(seq), List.of());

        assertThat(result.sequences().get(0).getCollectionDate()).contains(LocalDate.of(2020, 1, 1));
        assertThat(result.report().getSkippedCount()).isZero(); // it's not missing data, so not flagged as skipped
    }

    @Test
    void metadataRowWithBlankLocationAndHostDoesNotEraseHeaderDerivedOnes() {
        NucleotideSequence seq = new NucleotideSequence("s1", "ACGT")
                .withMetadata(LocalDate.of(2001, 1, 1), "France", "domestic_sheep"); // simulates header-derived location/host
        SampleMetadata meta = new SampleMetadata("s1", LocalDate.of(2001, 3, 1), null, null); // date given, location/host blank

        var result = MetadataJoiner.join(List.of(seq), List.of(meta));

        assertThat(result.sequences().get(0).getCollectionDate()).contains(LocalDate.of(2001, 3, 1)); // explicit date still wins
        assertThat(result.sequences().get(0).getLocation()).contains("France"); // but blank cells don't erase header values
        assertThat(result.sequences().get(0).getHost()).contains("domestic_sheep");
    }

    @Test
    void explicitMetadataLocationAndHostOverrideHeaderDerivedOnes() {
        NucleotideSequence seq = new NucleotideSequence("s1", "ACGT")
                .withMetadata(LocalDate.of(2001, 1, 1), "France", "domestic_sheep");
        SampleMetadata meta = new SampleMetadata("s1", null, "Spain", "wild_boar");

        var result = MetadataJoiner.join(List.of(seq), List.of(meta));

        assertThat(result.sequences().get(0).getLocation()).contains("Spain");
        assertThat(result.sequences().get(0).getHost()).contains("wild_boar");
    }

    @Test
    void sequenceWithNeitherHeaderDateNorMetadataIsFlaggedButNotDropped() {
        NucleotideSequence seq = new NucleotideSequence("s1", "ACGT");

        var result = MetadataJoiner.join(List.of(seq), List.of());

        assertThat(result.sequences()).hasSize(1);
        assertThat(result.sequences().get(0).getCollectionDate()).isEmpty();
        assertThat(result.report().getSkippedCount()).isEqualTo(1);
    }
}
