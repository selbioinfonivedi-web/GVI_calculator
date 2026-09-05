package org.gvi.algorithms.ri;

import org.gvi.core.exception.GviComputationException;
import org.gvi.core.exception.GviInputException;
import org.gvi.core.model.NucleotideSequence;
import org.gvi.core.model.SequenceAlignment;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.List;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.assertThatThrownBy;
import static org.assertj.core.api.Assertions.within;

class RecombinationIndexCalculatorTest {

    private final RecombinationIndexCalculator calc = new RecombinationIndexCalculator(100, 200, 42L);

    @Test
    void fourGameteTestDetectsCompatibleSites() {
        // only 2 of 4 combinations present (00, 11) -> compatible
        char[] site1 = {'A', 'A', 'G', 'G'};
        char[] site2 = {'C', 'C', 'T', 'T'};
        assertThat(calc.fourGameteIncompatible(site1, site2)).isFalse();
    }

    @Test
    void fourGameteTestDetectsIncompatibleSites() {
        // all 4 combinations present -> incompatible
        char[] site1 = {'A', 'A', 'G', 'G'};
        char[] site2 = {'C', 'T', 'C', 'T'};
        assertThat(calc.fourGameteIncompatible(site1, site2)).isTrue();
    }

    /**
     * Two informative-site "regions" driven by different, cross-cutting bipartitions
     * of the same 24 individuals -- the classic signature of a recombination
     * breakpoint at position 100: within a region, sites are perfectly linked
     * (compatible); across regions, the bipartitions cross-cut so all 4 gametes
     * appear (incompatible). Nearby (same-region) sites should be far more
     * compatible than the genome-wide average.
     */
    @Test
    void detectsStrongRecombinationSignalFromCrossCuttingBipartitions() {
        int[] region1Positions = {0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95};
        int[] region2Positions = {100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 185, 190, 195};

        List<NucleotideSequence> seqs = new ArrayList<>();
        for (int individual = 0; individual < 24; individual++) {
            char[] bases = new char[200];
            java.util.Arrays.fill(bases, 'A');
            char region1Allele = individual < 12 ? 'A' : 'G'; // bipartition #1: first half vs second half
            char region2Allele = individual % 2 == 0 ? 'A' : 'G'; // bipartition #2: even vs odd, cross-cuts #1
            for (int p : region1Positions) bases[p] = region1Allele;
            for (int p : region2Positions) bases[p] = region2Allele;
            seqs.add(new NucleotideSequence("ind" + individual, new String(bases)));
        }
        SequenceAlignment aln = SequenceAlignment.of(seqs);

        RiResult r = calc.computeFromAlignment(aln);
        assertThat(r.informativeSitesUsed()).isEqualTo(40);
        assertThat(r.ri()).isGreaterThan(0.3);
        assertThat(r.pValue()).isLessThan(0.05);
        assertThat(r.significant()).isTrue();
    }

    @Test
    void singleGenomeWideBipartitionGivesZeroRecombinationSignal() {
        // every informative site follows the SAME bipartition everywhere -> every pair compatible, ri exactly 0
        int[] allPositions = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180};
        List<NucleotideSequence> seqs = new ArrayList<>();
        for (int individual = 0; individual < 20; individual++) {
            char[] bases = new char[200];
            java.util.Arrays.fill(bases, 'A');
            char allele = individual < 10 ? 'A' : 'G';
            for (int p : allPositions) bases[p] = allele;
            seqs.add(new NucleotideSequence("ind" + individual, new String(bases)));
        }
        SequenceAlignment aln = SequenceAlignment.of(seqs);

        RiResult r = calc.computeFromAlignment(aln);
        assertThat(r.ri()).isCloseTo(0.0, within(1e-9));
        assertThat(r.significant()).isFalse();
    }

    @Test
    void throwsWhenTooFewInformativeSites() {
        List<NucleotideSequence> seqs = List.of(
                new NucleotideSequence("a", "AAAA"),
                new NucleotideSequence("b", "AAAA")
        );
        SequenceAlignment aln = SequenceAlignment.of(seqs);
        assertThatThrownBy(() -> calc.computeFromAlignment(aln)).isInstanceOf(GviComputationException.class);
    }

    @Test
    void isolateClassificationIsSimpleRatio() {
        RiResult r = calc.fromIsolateClassification(3, 30);
        assertThat(r.ri()).isCloseTo(0.1, within(1e-9));
        assertThat(r.method()).isEqualTo(RiMethod.ISOLATE_CLASSIFICATION);
        assertThat(r.category()).contains("Endemic multi-type");
    }

    @Test
    void isolateClassificationRejectsInvalidCounts() {
        assertThatThrownBy(() -> calc.fromIsolateClassification(5, 3)).isInstanceOf(GviInputException.class);
        assertThatThrownBy(() -> calc.fromIsolateClassification(1, 0)).isInstanceOf(GviInputException.class);
    }
}
