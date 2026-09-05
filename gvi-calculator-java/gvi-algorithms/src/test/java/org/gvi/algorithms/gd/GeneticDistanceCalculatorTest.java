package org.gvi.algorithms.gd;

import org.gvi.core.model.NucleotideSequence;
import org.junit.jupiter.api.Test;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.within;

class GeneticDistanceCalculatorTest {

    private final GeneticDistanceCalculator calc = new GeneticDistanceCalculator();

    @Test
    void hammingDistanceIsSimpleProportionOfDifferences() {
        // 10 sites, 2 differ -> p = 0.2
        NucleotideSequence a = new NucleotideSequence("a", "ACGTACGTAC");
        NucleotideSequence b = new NucleotideSequence("b", "ACGTACGTAG"); // last two differ: C->A, C->G... check below
        GdResult r = calc.compute(a, b, GdMethod.HAMMING);
        // recompute expected diffs manually
        String sa = "ACGTACGTAC", sb = "ACGTACGTAG";
        int diffs = 0;
        for (int i = 0; i < sa.length(); i++) if (sa.charAt(i) != sb.charAt(i)) diffs++;
        double expectedP = diffs / 10.0;
        assertThat(r.pDistance()).isCloseTo(expectedP, within(1e-9));
        assertThat(r.distance()).isCloseTo(expectedP, within(1e-9));
    }

    @Test
    void jukesCantorMatchesHandCalculatedFormula() {
        // p = 0.1 -> d = -3/4 * ln(1 - 4/3*0.1) = -0.75 * ln(0.86666...)
        // Construct 20-site sequences with exactly 2 differences (p=0.1)
        String ref = "A".repeat(20);
        StringBuilder q = new StringBuilder(ref);
        q.setCharAt(0, 'G');
        q.setCharAt(1, 'G');
        NucleotideSequence a = new NucleotideSequence("ref", ref);
        NucleotideSequence b = new NucleotideSequence("q", q.toString());

        GdResult r = calc.compute(a, b, GdMethod.JUKES_CANTOR);
        double p = 0.1;
        double expectedD = -0.75 * Math.log(1 - (4.0 / 3.0) * p);
        assertThat(r.pDistance()).isCloseTo(p, within(1e-9));
        assertThat(r.distance()).isCloseTo(expectedD, within(1e-9));
        assertThat(r.saturated()).isFalse();
    }

    @Test
    void jukesCantorFlagsSaturationInsteadOfCrashingOnNegativeLogArgument() {
        // p >= 0.75 makes (1 - 4/3 p) <= 0 -> undefined ln
        String ref = "AAAA";
        String q = "GGGG"; // p = 1.0
        NucleotideSequence a = new NucleotideSequence("ref", ref);
        NucleotideSequence b = new NucleotideSequence("q", q);

        GdResult r = calc.compute(a, b, GdMethod.JUKES_CANTOR);
        assertThat(r.saturated()).isTrue();
        assertThat(r.distance()).isNull();
        assertThat(r.pDistance()).isCloseTo(1.0, within(1e-9));
        assertThat(r.diagnostics()).anyMatch(s -> s.contains("saturated"));
    }

    @Test
    void kimura2ParameterDistinguishesTransitionsFromTransversions() {
        // 10 sites: 1 transition (A->G), 1 transversion (A->C), rest identical
        String ref = "A".repeat(10);
        StringBuilder q = new StringBuilder(ref);
        q.setCharAt(0, 'G'); // transition
        q.setCharAt(1, 'C'); // transversion
        NucleotideSequence a = new NucleotideSequence("ref", ref);
        NucleotideSequence b = new NucleotideSequence("q", q.toString());

        GdResult r = calc.compute(a, b, GdMethod.KIMURA_2_PARAMETER);
        double bigP = 0.1, bigQ = 0.1;
        double expected = -0.5 * Math.log(1 - 2 * bigP - bigQ) - 0.25 * Math.log(1 - 2 * bigQ);
        assertThat(r.distance()).isCloseTo(expected, within(1e-9));
    }

    @Test
    void pairwiseDeletionExcludesGapsAndAmbiguousSitesFromBothCounts() {
        // Position 0: gap in b -> excluded. Position 1: N in a -> excluded. Positions 2-4 real.
        NucleotideSequence a = new NucleotideSequence("a", "A" + "N" + "ACG");
        NucleotideSequence b = new NucleotideSequence("b", "-" + "A" + "ACA"); // last pos differs
        GdResult r = calc.compute(a, b, GdMethod.HAMMING);
        // comparable sites = positions 2,3,4 = 3; differences = 1 (last G vs A)
        assertThat(r.pDistance()).isCloseTo(1.0 / 3.0, within(1e-9));
    }

    @Test
    void identicalSequencesHaveZeroDistance() {
        NucleotideSequence a = new NucleotideSequence("a", "ACGTACGT");
        NucleotideSequence b = new NucleotideSequence("b", "ACGTACGT");
        GdResult r = calc.compute(a, b, GdMethod.JUKES_CANTOR);
        assertThat(r.distance()).isCloseTo(0.0, within(1e-9));
    }
}
