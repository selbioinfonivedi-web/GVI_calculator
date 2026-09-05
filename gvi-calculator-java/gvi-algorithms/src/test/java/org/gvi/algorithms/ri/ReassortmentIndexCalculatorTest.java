package org.gvi.algorithms.ri;

import org.gvi.core.exception.GviComputationException;
import org.gvi.core.model.NucleotideSequence;
import org.gvi.core.model.SequenceAlignment;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.List;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.assertThatThrownBy;

class ReassortmentIndexCalculatorTest {

    private final ReassortmentIndexCalculator calc = new ReassortmentIndexCalculator(500, 42L);

    /** 8 taxa split into two clusters {t0-t3} vs {t4-t7} by a block of cluster-defining SNPs. */
    private SequenceAlignment clusteredSegment(int[] group1, int[] group2) {
        List<NucleotideSequence> seqs = new ArrayList<>();
        for (int i = 0; i < 8; i++) {
            char[] bases = new char[60];
            java.util.Arrays.fill(bases, 'A');
            int taxon = i;
            boolean inGroup1 = java.util.Arrays.stream(group1).anyMatch(g -> g == taxon);
            char allele = inGroup1 ? 'A' : 'G';
            for (int p = 0; p < 30; p++) bases[p] = allele; // cluster-defining block
            seqs.add(new NucleotideSequence("t" + i, new String(bases)));
        }
        return SequenceAlignment.of(seqs);
    }

    @Test
    void concordantSegmentsGiveHighMantelRAndSignificantPValue() {
        int[] groupA1 = {0, 1, 2, 3};
        int[] groupA2 = {4, 5, 6, 7};
        SequenceAlignment segmentA = clusteredSegment(groupA1, groupA2);
        // Same bipartition on segment B (as expected under joint inheritance, no reassortment).
        SequenceAlignment segmentB = clusteredSegment(groupA1, groupA2);

        ReassortmentResult r = calc.compare("seg2", segmentA, "seg10", segmentB);

        assertThat(r.taxaCompared()).isEqualTo(8);
        assertThat(r.mantelR()).isGreaterThan(0.8);
        assertThat(r.pValue()).isLessThan(0.05);
        assertThat(r.significant()).isTrue();
        assertThat(r.reassortmentIndex()).isLessThan(0.2);
        assertThat(r.category()).contains("Concordant");
    }

    @Test
    void discordantSegmentsGiveLowMantelRAndNonSignificantPValue() {
        SequenceAlignment segmentA = clusteredSegment(new int[]{0, 1, 2, 3}, new int[]{4, 5, 6, 7});
        // A DIFFERENT, interleaved bipartition on segment B -- unrelated to segment A's grouping,
        // the signature of this specific taxon set having reassorted between the two segments.
        SequenceAlignment segmentB = clusteredSegment(new int[]{0, 2, 4, 6}, new int[]{1, 3, 5, 7});

        ReassortmentResult r = calc.compare("seg2", segmentA, "seg10", segmentB);

        assertThat(r.taxaCompared()).isEqualTo(8);
        assertThat(r.mantelR()).isLessThan(0.5);
        assertThat(r.significant()).isFalse();
        assertThat(r.reassortmentIndex()).isGreaterThan(0.5);
        assertThat(r.category()).contains("No significant concordance");
    }

    @Test
    void throwsWhenTooFewCommonTaxa() {
        List<NucleotideSequence> few = List.of(
                new NucleotideSequence("t0", "ACGTACGTAC"),
                new NucleotideSequence("t1", "ACGTACGTAG"),
                new NucleotideSequence("t2", "ACGTACGTAT"));
        SequenceAlignment segmentA = SequenceAlignment.of(few);
        SequenceAlignment segmentB = SequenceAlignment.of(few);
        assertThatThrownBy(() -> calc.compare("seg2", segmentA, "seg10", segmentB))
                .isInstanceOf(GviComputationException.class)
                .hasMessageContaining("at least 4");
    }

    @Test
    void throwsWhenNoSharedSequenceIds() {
        List<NucleotideSequence> a = List.of(
                new NucleotideSequence("a0", "ACGTACGTAC"), new NucleotideSequence("a1", "ACGTACGTAG"),
                new NucleotideSequence("a2", "ACGTACGTAT"), new NucleotideSequence("a3", "ACGTACGTAA"));
        List<NucleotideSequence> b = List.of(
                new NucleotideSequence("b0", "ACGTACGTAC"), new NucleotideSequence("b1", "ACGTACGTAG"),
                new NucleotideSequence("b2", "ACGTACGTAT"), new NucleotideSequence("b3", "ACGTACGTAA"));
        SequenceAlignment segmentA = SequenceAlignment.of(a);
        SequenceAlignment segmentB = SequenceAlignment.of(b);
        assertThatThrownBy(() -> calc.compare("seg2", segmentA, "seg10", segmentB))
                .isInstanceOf(GviComputationException.class)
                .hasMessageContaining("found only 0");
    }
}
