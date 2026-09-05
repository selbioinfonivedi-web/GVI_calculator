package org.gvi.algorithms.mu;

import org.gvi.algorithms.gd.GdMethod;
import org.gvi.algorithms.phylo.PhyloTreeFactory;
import org.gvi.core.model.NucleotideSequence;
import org.gvi.core.model.SequenceAlignment;
import org.gvi.core.util.TemporalUtil;
import org.junit.jupiter.api.Test;

import java.time.LocalDate;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.data.Percentage.withPercentage;

/**
 * Validates the native least-squares dating estimator against data
 * simulated under a KNOWN Jukes-Cantor molecular clock -- the same
 * simulation approach and same 0.002 sub/site/yr true rate as
 * {@link GroundTruthRecoveryTest}, so this estimator's recovery accuracy is
 * directly comparable to the existing root-to-tip regression estimator's
 * (0.5%/0.3% error) on the same underlying ground truth.
 */
class LeastSquaresDatingEstimatorTest {

    private static final long SEED = 20260813L;
    private static final double TRUE_MU = 2.0e-3;
    private static final int GENOME_LENGTH = 5000;
    private static final int TAXA = 30;
    private static final double MAX_YEARS = 10.0;
    private static final char[] BASES = {'A', 'C', 'G', 'T'};

    @Test
    void recoversTheTrueSimulatedClockRateWithinStatisticalTolerance() {
        Random rng = new Random(SEED);
        LocalDate baseDate = LocalDate.of(2020, 1, 1);

        char[] ancestral = new char[GENOME_LENGTH];
        for (int i = 0; i < GENOME_LENGTH; i++) ancestral[i] = BASES[rng.nextInt(4)];

        NucleotideSequence reference = new NucleotideSequence("ref", new String(ancestral)).withMetadata(baseDate, null, null);
        List<NucleotideSequence> all = new ArrayList<>();
        all.add(reference);
        Map<String, Double> tipDates = new HashMap<>();
        tipDates.put("ref", TemporalUtil.toDecimalYear(baseDate));

        for (int t = 0; t < TAXA; t++) {
            double elapsedYears = 0.5 + rng.nextDouble() * (MAX_YEARS - 0.5);
            char[] evolved = simulateJc69(ancestral, TRUE_MU, elapsedYears, rng);
            LocalDate date = baseDate.plusDays(Math.round(elapsedYears * 365.25));
            String id = "taxon" + t;
            all.add(new NucleotideSequence(id, new String(evolved)).withMetadata(date, null, null));
            tipDates.put(id, TemporalUtil.toDecimalYear(date));
        }

        SequenceAlignment alignment = SequenceAlignment.of(all, "ref");
        PhyloTreeFactory.BuiltTree built = PhyloTreeFactory.build(alignment, GdMethod.JUKES_CANTOR);

        LeastSquaresDatingEstimator.Estimate estimate = new LeastSquaresDatingEstimator().estimate(built.tree(), tipDates);
        System.out.println("[LeastSquaresDatingEstimatorTest] true mu=" + TRUE_MU + " recovered=" + estimate.muSubPerSiteYear()
                + " iterations=" + estimate.iterationsUsed() + " weightedSse=" + estimate.weightedSse());

        assertThat(estimate.muSubPerSiteYear()).isCloseTo(TRUE_MU, withPercentage(10));
        // Rooted AT a leaf (reference), not a separate virtual root node: n taxa (incl. reference)
        // + (n-2) internal non-taxon nodes = 2n-2 total, not the 2n-1 a standard internally-rooted binary tree would have.
        assertThat(estimate.nodeDates()).hasSize(2 * (TAXA + 1) - 2);
    }

    /**
     * Same deliberately-additive, noiseless known tree as
     * {@link EvolutionaryRateCalculatorTest#treeAwareMuUsesRealTreeStructureNotJustDistanceToReference}
     * (ref--1--X--3--Y, X--2--P, Y--1--Q, Y--4--R; every path's rate is
     * exactly 1e-5 sub/site/yr by construction) -- a clean, exact
     * correctness check on data with REAL shared internal structure (P
     * shares ancestry with Q/R through X; Q shares more recent ancestry
     * with R through Y), which is exactly the structure this estimator's
     * joint edge-respecting fit is meant to exploit better than a simple
     * root-to-tip regression that only sees (date, cumulative distance)
     * pairs.
     */
    @Test
    void recoversAnExactRateOnANoiselessHierarchicalTree() {
        int length = 100_000;
        char[] base = new char[length];
        java.util.Arrays.fill(base, 'A');
        int[] refX = {0};
        int[] xP = {1, 2};
        int[] xY = {3, 4, 5};
        int[] yQ = {6};
        int[] yR = {7, 8, 9, 10};

        NucleotideSequence ref = new NucleotideSequence("ref", new String(base)).withMetadata(LocalDate.of(2015, 1, 1), null, null);
        NucleotideSequence p = new NucleotideSequence("P", mutate(base, refX, xP)).withMetadata(LocalDate.of(2018, 1, 1), null, null);
        NucleotideSequence q = new NucleotideSequence("Q", mutate(base, refX, xY, yQ)).withMetadata(LocalDate.of(2020, 1, 1), null, null);
        NucleotideSequence r = new NucleotideSequence("R", mutate(base, refX, xY, yR)).withMetadata(LocalDate.of(2023, 1, 1), null, null);

        SequenceAlignment alignment = SequenceAlignment.of(List.of(ref, p, q, r), "ref");
        PhyloTreeFactory.BuiltTree built = PhyloTreeFactory.build(alignment, GdMethod.JUKES_CANTOR);

        Map<String, Double> tipDates = new HashMap<>();
        for (NucleotideSequence seq : List.of(ref, p, q, r)) {
            tipDates.put(seq.getId(), TemporalUtil.toDecimalYear(seq.getCollectionDate().get()));
        }

        LeastSquaresDatingEstimator.Estimate estimate = new LeastSquaresDatingEstimator().estimate(built.tree(), tipDates);
        System.out.println("[LeastSquaresDatingEstimatorTest] hierarchical: recovered=" + estimate.muSubPerSiteYear());

        assertThat(estimate.muSubPerSiteYear()).isCloseTo(1e-5, withPercentage(5));
    }

    private static String mutate(char[] base, int[]... blocks) {
        char[] copy = base.clone();
        for (int[] block : blocks) {
            for (int pos : block) copy[pos] = 'G';
        }
        return new String(copy);
    }

    private static char[] simulateJc69(char[] ancestral, double mu, double elapsedYears, Random rng) {
        double pSame = 0.25 + 0.75 * Math.exp(-(4.0 / 3.0) * mu * elapsedYears);
        char[] evolved = new char[ancestral.length];
        for (int i = 0; i < ancestral.length; i++) {
            if (rng.nextDouble() < pSame) {
                evolved[i] = ancestral[i];
            } else {
                char newBase;
                do {
                    newBase = BASES[rng.nextInt(4)];
                } while (newBase == ancestral[i]);
                evolved[i] = newBase;
            }
        }
        return evolved;
    }
}
