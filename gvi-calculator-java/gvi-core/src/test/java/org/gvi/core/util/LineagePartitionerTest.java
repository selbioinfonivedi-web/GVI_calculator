package org.gvi.core.util;

import org.gvi.core.model.NucleotideSequence;
import org.gvi.core.model.SequenceAlignment;
import org.junit.jupiter.api.Test;

import java.util.List;

import static org.assertj.core.api.Assertions.assertThat;

class LineagePartitionerTest {

    private static String repeat(String unit, int times) {
        return unit.repeat(times);
    }

    @Test
    void separatesTwoDivergentGroups() {
        // Two clusters ~50% apart; within each, sequences differ by a single base.
        String a = repeat("ACGT", 25);
        String b = repeat("TGCA", 25);
        var aln = SequenceAlignment.of(List.of(
                new NucleotideSequence("a1", a),
                new NucleotideSequence("a2", a.substring(0, 99) + "A"),
                new NucleotideSequence("b1", b),
                new NucleotideSequence("b2", b.substring(0, 99) + "T")), null);

        var r = LineagePartitioner.detect(aln);

        assertThat(r.partitioned()).isTrue();
        assertThat(r.clusters()).hasSize(2);
        assertThat(r.description()).contains("Lineage structure detected").contains("group");
    }

    @Test
    void oneTightGroupIsNotPartitioned() {
        String base = repeat("ACGT", 25);
        var aln = SequenceAlignment.of(List.of(
                new NucleotideSequence("s1", base),
                new NucleotideSequence("s2", base.substring(0, 99) + "A"),
                new NucleotideSequence("s3", base.substring(0, 98) + "AT")), null);

        var r = LineagePartitioner.detect(aln);

        assertThat(r.partitioned()).isFalse();
        assertThat(r.description()).contains("No lineage structure detected");
    }

    @Test
    void everySequenceAppearsInExactlyOneCluster() {
        String a = repeat("ACGT", 25);
        String b = repeat("TGCA", 25);
        var aln = SequenceAlignment.of(List.of(
                new NucleotideSequence("a1", a),
                new NucleotideSequence("b1", b),
                new NucleotideSequence("a2", a)), null);

        var r = LineagePartitioner.detect(aln);

        var all = r.clusters().stream().flatMap(c -> c.sequenceIds().stream()).toList();
        assertThat(all).containsExactlyInAnyOrder("a1", "a2", "b1");
    }

    @Test
    void tooFewSequencesReportsNoStructureRatherThanFailing() {
        var aln = SequenceAlignment.of(List.of(
                new NucleotideSequence("s1", "ACGTACGT"),
                new NucleotideSequence("s2", "ACGTACGA")), null);

        var r = LineagePartitioner.detect(aln);

        assertThat(r.partitioned()).isFalse();
        assertThat(r.description()).contains("Too few sequences");
    }
}
