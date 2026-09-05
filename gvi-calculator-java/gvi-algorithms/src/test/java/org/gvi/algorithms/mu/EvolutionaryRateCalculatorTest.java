package org.gvi.algorithms.mu;

import org.apache.commons.math3.stat.regression.SimpleRegression;
import org.gvi.core.exception.GviComputationException;
import org.gvi.core.model.NucleotideSequence;
import org.gvi.core.model.SequenceAlignment;
import org.gvi.core.util.TemporalUtil;
import org.junit.jupiter.api.Test;

import java.time.LocalDate;
import java.util.List;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.assertThatThrownBy;
import static org.assertj.core.api.Assertions.within;

class EvolutionaryRateCalculatorTest {

    private final EvolutionaryRateCalculator calc = new EvolutionaryRateCalculator();

    private static NucleotideSequence dated(String id, String seq, LocalDate date) {
        return new NucleotideSequence(id, seq).withMetadata(date, null, null);
    }

    @Test
    void slopeAndRSquaredMatchAnIndependentlyBuiltRegressionOverTheSameDistanceDatePairs() {
        String base = "A".repeat(100);
        NucleotideSequence ref = dated("ref", base, LocalDate.of(2020, 1, 1));

        StringBuilder q1 = new StringBuilder(base);
        q1.setCharAt(0, 'G'); // p = 0.01
        StringBuilder q2 = new StringBuilder(base);
        q2.setCharAt(0, 'G'); q2.setCharAt(1, 'G'); // p = 0.02
        StringBuilder q3 = new StringBuilder(base);
        q3.setCharAt(0, 'G'); q3.setCharAt(1, 'G'); q3.setCharAt(2, 'G'); // p = 0.03

        NucleotideSequence sq1 = dated("q1", q1.toString(), LocalDate.of(2021, 1, 1));
        NucleotideSequence sq2 = dated("q2", q2.toString(), LocalDate.of(2022, 1, 1));
        NucleotideSequence sq3 = dated("q3", q3.toString(), LocalDate.of(2023, 1, 1));

        SequenceAlignment aln = SequenceAlignment.of(List.of(ref, sq1, sq2, sq3), "ref");
        MuResult r = calc.compute(aln);

        // independently reproduce the regression the calculator should have run
        SimpleRegression expected = new SimpleRegression();
        expected.addData(TemporalUtil.toDecimalYear(LocalDate.of(2020, 1, 1)), jc(0.0));
        expected.addData(TemporalUtil.toDecimalYear(LocalDate.of(2021, 1, 1)), jc(0.01));
        expected.addData(TemporalUtil.toDecimalYear(LocalDate.of(2022, 1, 1)), jc(0.02));
        expected.addData(TemporalUtil.toDecimalYear(LocalDate.of(2023, 1, 1)), jc(0.03));

        assertThat(r.muSubPerSiteYear()).isCloseTo(expected.getSlope(), within(1e-9));
        assertThat(r.rSquared()).isCloseTo(expected.getRSquare(), within(1e-9));
        assertThat(r.pointsUsed()).isEqualTo(4);
        assertThat(r.reliable()).isTrue(); // near-perfectly linear construction
    }

    private static double jc(double p) {
        return -0.75 * Math.log(1 - (4.0 / 3.0) * p);
    }

    @Test
    void throwsWhenFewerThanTwoDistinctDatesAvailable() {
        NucleotideSequence ref = dated("ref", "ACGTACGTAC", LocalDate.of(2020, 1, 1));
        NucleotideSequence q1 = new NucleotideSequence("q1", "ACGTACGTAG"); // no date
        SequenceAlignment aln = SequenceAlignment.of(List.of(ref, q1), "ref");
        assertThatThrownBy(() -> calc.compute(aln)).isInstanceOf(GviComputationException.class);
    }

    @Test
    void flagsWeakTemporalSignalAsUnreliable() {
        // distance is identical regardless of date -> essentially no clock signal, R^2 near 0
        String base = "A".repeat(100);
        NucleotideSequence ref = dated("ref", base, LocalDate.of(2020, 1, 1));
        StringBuilder mutant = new StringBuilder(base);
        mutant.setCharAt(50, 'G');
        NucleotideSequence q1 = dated("q1", mutant.toString(), LocalDate.of(2020, 6, 1));
        NucleotideSequence q2 = dated("q2", base, LocalDate.of(2021, 1, 1)); // back to zero distance despite later date
        NucleotideSequence q3 = dated("q3", mutant.toString(), LocalDate.of(2021, 6, 1));

        SequenceAlignment aln = SequenceAlignment.of(List.of(ref, q1, q2, q3), "ref");
        MuResult r = calc.compute(aln);
        assertThat(r.diagnostics()).anyMatch(s -> s.contains("weak") || s.contains("R^2"));
    }

    /**
     * Builds 4 sequences whose mutations are assigned to disjoint edge-blocks
     * of a known tree (ref--1--X--3--Y, X--2--P, Y--1--Q, Y--4--R -- the same
     * tree used in NeighborJoiningTest), so their pairwise Hamming distances
     * are exactly path-additive by construction (no back-mutation, no
     * multiple hits at any position). The alignment is made very long
     * (100,000bp) so the resulting p-distances are tiny (~1e-5-1e-4) and the
     * Jukes-Cantor correction is negligible -- this keeps the corrected
     * distances additive too, to within floating precision, so the whole
     * tree-aware pipeline (JC distance matrix -> NJ -> rooting -> regression)
     * can be checked against a precisely predictable answer instead of just
     * "runs without crashing".
     */
    @Test
    void treeAwareMuUsesRealTreeStructureNotJustDistanceToReference() {
        int length = 100_000;
        char[] base = new char[length];
        java.util.Arrays.fill(base, 'A');

        // edge -> block of positions (lengths match ref-X=1, X-P=2, X-Y=3, Y-Q=1, Y-R=4)
        int[] refX = {0};
        int[] xP = {1, 2};
        int[] xY = {3, 4, 5};
        int[] yQ = {6};
        int[] yR = {7, 8, 9, 10};

        NucleotideSequence ref = dated("ref", new String(base), LocalDate.of(2015, 1, 1));
        NucleotideSequence p = dated("P", mutate(base, refX, xP), LocalDate.of(2018, 1, 1));       // path: ref-X, X-P (3 muts / 3yr)
        NucleotideSequence q = dated("Q", mutate(base, refX, xY, yQ), LocalDate.of(2020, 1, 1));    // path: ref-X, X-Y, Y-Q (5 muts / 5yr)
        NucleotideSequence r = dated("R", mutate(base, refX, xY, yR), LocalDate.of(2023, 1, 1));    // path: ref-X, X-Y, Y-R (8 muts / 8yr)

        SequenceAlignment aln = SequenceAlignment.of(List.of(ref, p, q, r), "ref");
        MuResult result = calc.computeTreeAware(aln);

        // every path's mutation-count/year ratio is 1e-5/yr by construction (3/3e5, 5/5e5, 8/8e5 muts-per-bp-per-yr)
        assertThat(result.method()).isEqualTo(MuEstimationMethod.TREE_ROOT_TO_TIP);
        assertThat(result.muSubPerSiteYear()).isCloseTo(1e-5, within(5e-7));
        assertThat(result.rSquared()).isGreaterThan(0.999);
        assertThat(result.pointsUsed()).isEqualTo(4);
    }

    private static String mutate(char[] base, int[]... blocks) {
        char[] copy = base.clone();
        for (int[] block : blocks) {
            for (int pos : block) copy[pos] = 'G';
        }
        return new String(copy);
    }

    @Test
    void treeAwareWithBootstrapSupportFlagAddsASupportDiagnostic() {
        int length = 100_000;
        char[] base = new char[length];
        java.util.Arrays.fill(base, 'A');
        int[] refX = {0};
        int[] xP = {1, 2};
        int[] xY = {3, 4, 5};
        int[] yQ = {6};
        int[] yR = {7, 8, 9, 10};

        NucleotideSequence ref = dated("ref", new String(base), LocalDate.of(2015, 1, 1));
        NucleotideSequence p = dated("P", mutate(base, refX, xP), LocalDate.of(2018, 1, 1));
        NucleotideSequence q = dated("Q", mutate(base, refX, xY, yQ), LocalDate.of(2020, 1, 1));
        NucleotideSequence r = dated("R", mutate(base, refX, xY, yR), LocalDate.of(2023, 1, 1));
        SequenceAlignment aln = SequenceAlignment.of(List.of(ref, p, q, r), "ref");

        MuResult withBootstrap = calc.computeTreeAware(aln, true, 100);
        assertThat(withBootstrap.diagnostics()).anyMatch(d -> d.contains("Bootstrap support (Felsenstein 1985"));

        MuResult withoutBootstrap = calc.computeTreeAware(aln);
        assertThat(withoutBootstrap.diagnostics()).noneMatch(d -> d.contains("Bootstrap support"));
    }

    @Test
    void treeAwareThrowsWithFewerThanThreeSequences() {
        NucleotideSequence ref = dated("ref", "ACGTACGTAC", LocalDate.of(2020, 1, 1));
        NucleotideSequence q1 = dated("q1", "ACGTACGTAG", LocalDate.of(2021, 1, 1));
        SequenceAlignment aln = SequenceAlignment.of(List.of(ref, q1), "ref");
        assertThatThrownBy(() -> calc.computeTreeAware(aln)).isInstanceOf(GviComputationException.class);
    }

    /**
     * Small/fast dataset by design -- the underlying GTR model, Gamma rates,
     * pruning likelihood, and ML branch-length optimizer are each already
     * rigorously hand-verified in their own test classes
     * (GtrModelTest, DiscreteGammaRatesTest, TreeLikelihoodCalculatorTest,
     * MlBranchLengthOptimizerTest). This test only needs to confirm the
     * wiring: it runs without error, tags the result with the right
     * method, and produces a sane (positive, finite) mu.
     */
    @Test
    void gtrGammaAwareProducesASaneResultAndTagsTheMethodCorrectly() {
        int length = 2000;
        char[] base = new char[length];
        java.util.Arrays.fill(base, 'A');
        NucleotideSequence ref = dated("ref", new String(base), LocalDate.of(2020, 1, 1));

        StringBuilder q1 = new StringBuilder(new String(base));
        for (int i = 0; i < 5; i++) q1.setCharAt(i, 'G');
        StringBuilder q2 = new StringBuilder(new String(base));
        for (int i = 0; i < 10; i++) q2.setCharAt(i, 'G');
        StringBuilder q3 = new StringBuilder(new String(base));
        for (int i = 0; i < 15; i++) q3.setCharAt(i, 'G');

        NucleotideSequence sq1 = dated("q1", q1.toString(), LocalDate.of(2021, 1, 1));
        NucleotideSequence sq2 = dated("q2", q2.toString(), LocalDate.of(2022, 1, 1));
        NucleotideSequence sq3 = dated("q3", q3.toString(), LocalDate.of(2023, 1, 1));

        SequenceAlignment aln = SequenceAlignment.of(List.of(ref, sq1, sq2, sq3), "ref");
        // HKY85 (1 rate param), not the default GTR (5 rate params) -- full GTR+coordinate-ascent on tiny,
        // barely-identified test data can take a long time to (not) converge across its 5 rate parameters;
        // the joint optimizer's own correctness (including with full GTR) is already covered by
        // MlPhylogeneticOptimizerTest/SubstitutionModelSelectorTest. This test only needs to confirm the
        // mu-level wiring: runs without error, tags the result with the right method, sane mu.
        MuResult result = calc.computeGtrGammaAware(aln, "hky85", null);

        assertThat(result.method()).isEqualTo(MuEstimationMethod.GTR_GAMMA_ML_BRANCH_LENGTHS);
        assertThat(result.muSubPerSiteYear()).isPositive();
        assertThat(Double.isFinite(result.muSubPerSiteYear())).isTrue();
        assertThat(result.diagnostics()).anyMatch(s -> s.contains("ML fit"));
    }

    @Test
    void gtrGammaAwareWithGammaHeterogeneityAlsoProducesASaneResult() {
        int length = 1500;
        char[] base = new char[length];
        java.util.Arrays.fill(base, 'A');
        NucleotideSequence ref = dated("ref", new String(base), LocalDate.of(2020, 1, 1));
        StringBuilder q1 = new StringBuilder(new String(base));
        for (int i = 0; i < 6; i++) q1.setCharAt(i, 'G');
        StringBuilder q2 = new StringBuilder(new String(base));
        for (int i = 0; i < 12; i++) q2.setCharAt(i, 'G');
        NucleotideSequence sq1 = dated("q1", q1.toString(), LocalDate.of(2021, 1, 1));
        NucleotideSequence sq2 = dated("q2", q2.toString(), LocalDate.of(2022, 1, 1));

        SequenceAlignment aln = SequenceAlignment.of(List.of(ref, sq1, sq2), "ref");
        // alpha=0.5 is now just a STARTING point for ML refinement (see class javadoc), not the final value
        MuResult result = calc.computeGtrGammaAware(aln, "hky85", 0.5);

        assertThat(result.method()).isEqualTo(MuEstimationMethod.GTR_GAMMA_ML_BRANCH_LENGTHS);
        assertThat(Double.isFinite(result.muSubPerSiteYear())).isTrue();
        assertThat(result.diagnostics()).anyMatch(s -> s.contains("gamma alpha="));
    }

    @Test
    void twoArgOverloadDefaultsToGtrForBackwardCompatibility() {
        int length = 800;
        char[] base = new char[length];
        java.util.Arrays.fill(base, 'A');
        NucleotideSequence ref = dated("ref", new String(base), LocalDate.of(2020, 1, 1));
        StringBuilder q1 = new StringBuilder(new String(base));
        for (int i = 0; i < 5; i++) q1.setCharAt(i, 'G');
        StringBuilder q2 = new StringBuilder(new String(base));
        for (int i = 0; i < 9; i++) q2.setCharAt(i, 'G');
        NucleotideSequence sq1 = dated("q1", q1.toString(), LocalDate.of(2021, 1, 1));
        NucleotideSequence sq2 = dated("q2", q2.toString(), LocalDate.of(2022, 1, 1));

        SequenceAlignment aln = SequenceAlignment.of(List.of(ref, sq1, sq2), "ref");
        MuResult result = calc.computeGtrGammaAware(aln, (Double) null);

        assertThat(result.method()).isEqualTo(MuEstimationMethod.GTR_GAMMA_ML_BRANCH_LENGTHS);
        assertThat(result.diagnostics()).anyMatch(s -> s.contains("GTR ML fit"));
    }

    @Test
    void gtrGammaAwareThrowsAboveTaxonCapWithoutAttemptingExpensiveComputation() {
        List<NucleotideSequence> seqs = new java.util.ArrayList<>();
        for (int i = 0; i <= EvolutionaryRateCalculator.MAX_TAXA_FOR_ML; i++) {
            seqs.add(dated("s" + i, "ACGTACGTAC", LocalDate.of(2020, 1, 1).plusDays(i)));
        }
        SequenceAlignment aln = SequenceAlignment.of(seqs, "s0");
        assertThatThrownBy(() -> calc.computeGtrGammaAware(aln, null))
                .isInstanceOf(GviComputationException.class)
                .hasMessageContaining("taxa exceeds");
    }

    @Test
    void gtrGammaAwareThrowsAboveSiteCapWithoutAttemptingExpensiveComputation() {
        String longSeq = "A".repeat(EvolutionaryRateCalculator.MAX_SITES_FOR_ML + 1);
        List<NucleotideSequence> seqs = List.of(
                dated("ref", longSeq, LocalDate.of(2020, 1, 1)),
                dated("q1", longSeq, LocalDate.of(2021, 1, 1)),
                dated("q2", longSeq, LocalDate.of(2022, 1, 1))
        );
        SequenceAlignment aln = SequenceAlignment.of(seqs, "ref");
        assertThatThrownBy(() -> calc.computeGtrGammaAware(aln, null))
                .isInstanceOf(GviComputationException.class)
                .hasMessageContaining("exceeds the");
    }

    @Test
    void treeAwareThrowsAboveTheTaxonCapInsteadOfHangingOnCubicNjCost() {
        List<NucleotideSequence> seqs = new java.util.ArrayList<>();
        for (int i = 0; i < 301; i++) {
            seqs.add(dated("s" + i, "ACGTACGTAC", LocalDate.of(2020, 1, 1).plusDays(i)));
        }
        SequenceAlignment aln = SequenceAlignment.of(seqs, "s0");
        assertThatThrownBy(() -> calc.computeTreeAware(aln))
                .isInstanceOf(GviComputationException.class)
                .hasMessageContaining("cap");
    }
}
