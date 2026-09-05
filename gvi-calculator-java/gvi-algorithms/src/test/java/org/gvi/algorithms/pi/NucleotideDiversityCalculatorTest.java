package org.gvi.algorithms.pi;

import org.gvi.core.exception.GviComputationException;
import org.gvi.core.model.NucleotideSequence;
import org.gvi.core.model.SequenceAlignment;
import org.junit.jupiter.api.Test;

import java.util.List;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.assertThatThrownBy;
import static org.assertj.core.api.Assertions.within;

class NucleotideDiversityCalculatorTest {

    // small bootstrap replicate count for fast, deterministic tests
    private final NucleotideDiversityCalculator calc = new NucleotideDiversityCalculator(1234L, 50, 100_000L);

    @Test
    void identicalSequencesGiveZeroDiversity() {
        SequenceAlignment aln = SequenceAlignment.of(List.of(
                new NucleotideSequence("s1", "ACGTACGTAC"),
                new NucleotideSequence("s2", "ACGTACGTAC"),
                new NucleotideSequence("s3", "ACGTACGTAC")
        ));
        PiResult r = calc.compute(aln);
        assertThat(r.pi()).isCloseTo(0.0, within(1e-9));
    }

    @Test
    void matchesHandCalculatedValueForThreeSequences() {
        // s1 vs s2: 1 diff / 10 = 0.1 ; s1 vs s3: 2 diff /10 = 0.2 ; s2 vs s3: 1 diff /10 = 0.1
        // average = (0.1+0.2+0.1)/3 = 0.13333...
        NucleotideSequence s1 = new NucleotideSequence("s1", "AAAAAAAAAA");
        NucleotideSequence s2 = new NucleotideSequence("s2", "AAAAAAAAAG"); // 1 diff from s1
        NucleotideSequence s3 = new NucleotideSequence("s3", "AAAAAAAAGG"); // 2 diff from s1, 1 diff from s2
        SequenceAlignment aln = SequenceAlignment.of(List.of(s1, s2, s3));

        PiResult r = calc.compute(aln);
        assertThat(r.pi()).isCloseTo(0.4 / 3.0, within(1e-9));
        assertThat(r.sequenceCount()).isEqualTo(3);
        assertThat(r.alignmentLength()).isEqualTo(10);
        assertThat(r.confidenceInterval()).isNotNull();
    }

    @Test
    void requiresAtLeastTwoSequences() {
        SequenceAlignment aln = SequenceAlignment.of(List.of(new NucleotideSequence("only", "ACGT")));
        assertThatThrownBy(() -> calc.compute(aln)).isInstanceOf(GviComputationException.class);
    }

    @Test
    void classifiesLowDiversityAsAcuteOutbreakBand() {
        // pi ~ 0.0001 should classify as "Acute outbreak"
        String base = "A".repeat(10000);
        StringBuilder mutant = new StringBuilder(base);
        mutant.setCharAt(0, 'G'); // 1 diff / 10000 = 0.0001
        SequenceAlignment aln = SequenceAlignment.of(List.of(
                new NucleotideSequence("s1", base),
                new NucleotideSequence("s2", mutant.toString())
        ));
        PiResult r = calc.compute(aln);
        assertThat(r.category()).contains("Acute outbreak");
    }
}
