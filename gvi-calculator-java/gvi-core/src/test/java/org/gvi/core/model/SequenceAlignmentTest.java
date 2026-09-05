package org.gvi.core.model;

import org.gvi.core.exception.GviInputException;
import org.junit.jupiter.api.Test;

import java.util.List;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.assertThatThrownBy;

class SequenceAlignmentTest {

    @Test
    void acceptsEqualLengthSequences() {
        SequenceAlignment aln = SequenceAlignment.of(List.of(
                new NucleotideSequence("ref", "ACGTACGT"),
                new NucleotideSequence("q1", "ACGTACGA")
        ));
        assertThat(aln.size()).isEqualTo(2);
        assertThat(aln.length()).isEqualTo(8);
        assertThat(aln.getReference().getId()).isEqualTo("ref");
        assertThat(aln.getQueries()).hasSize(1);
    }

    @Test
    void rejectsMismatchedLengthsWithActionableMessage() {
        assertThatThrownBy(() -> SequenceAlignment.of(List.of(
                new NucleotideSequence("ref", "ACGTACGT"),
                new NucleotideSequence("q1", "ACGT")
        ))).isInstanceOf(GviInputException.class)
          .hasMessageContaining("q1")
          .hasMessageContaining("align your sequences first");
    }

    @Test
    void selectsNamedReference() {
        SequenceAlignment aln = SequenceAlignment.of(List.of(
                new NucleotideSequence("a", "ACGT"),
                new NucleotideSequence("b", "ACGA")
        ), "b");
        assertThat(aln.getReference().getId()).isEqualTo("b");
    }

    @Test
    void rejectsEmptyAlignment() {
        assertThatThrownBy(() -> SequenceAlignment.of(List.of()))
                .isInstanceOf(GviInputException.class);
    }
}
