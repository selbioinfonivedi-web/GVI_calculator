package org.gvi.core.util;

import org.gvi.core.model.NucleotideSequence;
import org.gvi.core.model.SequenceAlignment;
import org.junit.jupiter.api.Test;

import java.time.LocalDate;
import java.util.List;

import static org.assertj.core.api.Assertions.assertThat;

class AlignmentTrimmerTest {

    @Test
    void trimsColumnsWhereAnySequenceLacksCoverage() {
        // s2 covers only the middle; the flanks are absent data, not observed difference.
        var aln = SequenceAlignment.of(List.of(
                new NucleotideSequence("ref", "AAACCCGGGTTT"),
                new NucleotideSequence("s1", "AAACCCGGGTTT"),
                new NucleotideSequence("s2", "---CCCGGG---")), null);

        var r = AlignmentTrimmer.trimToFullyCovered(aln, 1);

        assertThat(r.trimmed()).isTrue();
        assertThat(r.originalColumns()).isEqualTo(12);
        assertThat(r.retainedColumns()).isEqualTo(6);
        assertThat(r.longestContiguousBlock()).isEqualTo(6);
        assertThat(r.alignment().getSequences()).allSatisfy(s ->
                assertThat(s.getSequence()).isEqualTo("CCCGGG"));
    }

    @Test
    void leavesAFullyCoveredAlignmentAlone() {
        var aln = SequenceAlignment.of(List.of(
                new NucleotideSequence("ref", "ACGTACGT"),
                new NucleotideSequence("s1", "ACGTACGA")), null);

        var r = AlignmentTrimmer.trimToFullyCovered(aln, 1);

        assertThat(r.trimmed()).isFalse();
        assertThat(r.description()).contains("No trimming needed");
        assertThat(r.alignment()).isSameAs(aln);
    }

    @Test
    void refusesToTrimBelowTheMinimumRatherThanManufactureATinyAlignment() {
        var aln = SequenceAlignment.of(List.of(
                new NucleotideSequence("ref", "AAAACCCCGGGG"),
                new NucleotideSequence("s1", "----CCCC----")), null);

        var r = AlignmentTrimmer.trimToFullyCovered(aln, 100);

        assertThat(r.trimmed()).isFalse();
        assertThat(r.alignment()).isSameAs(aln);
        assertThat(r.description()).contains("Trimming NOT applied").contains("overlap too little");
    }

    @Test
    void preservesIdsOrderAndMetadata() {
        var aln = SequenceAlignment.of(List.of(
                new NucleotideSequence("ref", "AAACCCTTT", LocalDate.of(2020, 1, 1), "IN", "cattle"),
                new NucleotideSequence("s1", "---CCC---", LocalDate.of(2021, 6, 2), "KE", "sheep")), null);

        var r = AlignmentTrimmer.trimToFullyCovered(aln, 1);

        assertThat(r.alignment().getSequences()).extracting(NucleotideSequence::getId).containsExactly("ref", "s1");
        var s1 = r.alignment().getSequences().get(1);
        assertThat(s1.getCollectionDate()).contains(LocalDate.of(2021, 6, 2));
        assertThat(s1.getLocation()).contains("KE");
        assertThat(s1.getHost()).contains("sheep");
        assertThat(r.alignment().getReference().getId()).isEqualTo("ref");
    }

    @Test
    void reportsLongestContiguousBlockSeparatelyFromTotalRetained() {
        // Two covered islands: 3 + 2 columns retained, but the longest usable run is only 3.
        var aln = SequenceAlignment.of(List.of(
                new NucleotideSequence("ref", "AAACCCGGGTT"),
                new NucleotideSequence("s1", "AAA---GG-TT")), null);

        var r = AlignmentTrimmer.trimToFullyCovered(aln, 1);

        assertThat(r.retainedColumns()).isEqualTo(7);
        assertThat(r.longestContiguousBlock()).isEqualTo(3);
    }

    @Test
    void remapsCoordinatesOntoTheTrimmedNumbering() {
        // Columns 4-9 (1-based) survive; a gene spanning the whole original must become [1,6].
        var aln = SequenceAlignment.of(List.of(
                new NucleotideSequence("ref", "AAACCCGGGTTT"),
                new NucleotideSequence("s1", "---CCCGGG---")), null);

        var r = AlignmentTrimmer.trimToFullyCovered(aln, 1);

        assertThat(r.remap(1, 12)).hasValueSatisfying(v -> assertThat(v).containsExactly(1L, 6L));
        // A gene wholly inside the discarded flank cannot be remapped.
        assertThat(r.remap(1, 3)).isEmpty();
        // A gene inside the kept window maps to the corresponding offset.
        assertThat(r.remap(4, 6)).hasValueSatisfying(v -> assertThat(v).containsExactly(1L, 3L));
    }

    @Test
    void remapIsIdentityWhenNothingWasTrimmed() {
        var aln = SequenceAlignment.of(List.of(
                new NucleotideSequence("ref", "ACGTACGT"),
                new NucleotideSequence("s1", "ACGTACGA")), null);

        var r = AlignmentTrimmer.trimToFullyCovered(aln, 1);

        assertThat(r.remap(2, 5)).hasValueSatisfying(v -> assertThat(v).containsExactly(2L, 5L));
    }
}
