package org.gvi.algorithms.mu;

import org.gvi.algorithms.gd.GdMethod;
import org.gvi.algorithms.phylo.PhyloTreeFactory;
import org.gvi.core.model.NucleotideSequence;
import org.gvi.core.model.SequenceAlignment;
import org.gvi.core.util.TemporalUtil;
import org.junit.jupiter.api.Test;

import java.time.LocalDate;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.data.Percentage.withPercentage;

/**
 * Ground-truth recovery for the relaxed clock -- the two things that must both be true for this to
 * be a real relaxed clock rather than a strict clock with extra machinery: (1) on data with one true
 * rate everywhere, the fitted coefficient of variation must come out near zero, and (2) on data with
 * a genuinely faster branch, it must come out clearly above that near-zero baseline. Neither alone
 * would catch an estimator that always reports {@code sigma=0} (indistinguishable from
 * {@link LeastSquaresDatingEstimator} on test 1) or one that reports noisy nonsense regardless of
 * input (which could still pass a single "coefficient of variation > 0" check).
 * <p>
 * Same deliberately-additive, noiseless topology and dates as {@link
 * LeastSquaresDatingEstimatorTest#recoversAnExactRateOnANoiselessHierarchicalTree} (ref--X--P,
 * X--Y--{Q,R}) so the two estimators' behaviour on the same underlying structure is directly
 * comparable.
 */
class RelaxedClockMlEstimatorTest {

    @Test
    void reportsANearZeroCoefficientOfVariationWhenEveryBranchTrulySharesOneRate() {
        // Every root-to-leaf path accumulates exactly 1e-5 sub/site/yr, precisely as in
        // LeastSquaresDatingEstimatorTest's noiseless fixture -- no true among-branch rate variation.
        int length = 100_000;
        int[] refX = {0};
        int[] xP = {1, 2};
        int[] xY = {3, 4, 5};
        int[] yQ = {6};
        int[] yR = {7, 8, 9, 10};

        RelaxedClockMlEstimator.Estimate estimate = fit(length,
                new LocalDate[]{LocalDate.of(2015, 1, 1), LocalDate.of(2018, 1, 1), LocalDate.of(2020, 1, 1), LocalDate.of(2023, 1, 1)},
                refX, xP, xY, yQ, yR);

        System.out.println("[RelaxedClockMlEstimatorTest] homogeneous: rate0=" + estimate.rate0SubPerSiteYear()
                + " sigma=" + estimate.sigma() + " CV=" + estimate.coefficientOfVariation());

        assertThat(estimate.rate0SubPerSiteYear()).isCloseTo(1e-5, withPercentage(5));
        assertThat(estimate.coefficientOfVariation()).isLessThan(0.05);
    }

    @Test
    void reportsAClearlyHigherCoefficientOfVariationWhenOneBranchTrulyEvolvesFaster() {
        // Same topology and dates, but Q's branch (yQ) carries 8x the mutations of the homogeneous
        // fixture over the same elapsed time -- a real, single fast-evolving lineage, the exact
        // scenario a strict clock cannot represent and a relaxed clock exists to describe.
        int length = 100_000;
        int[] refX = {0};
        int[] xP = {1, 2};
        int[] xY = {3, 4, 5};
        int[] yQ = {6, 7, 8, 9, 10, 11, 12, 13};
        int[] yR = {14, 15, 16, 17};

        RelaxedClockMlEstimator.Estimate homogeneous = fit(length,
                new LocalDate[]{LocalDate.of(2015, 1, 1), LocalDate.of(2018, 1, 1), LocalDate.of(2020, 1, 1), LocalDate.of(2023, 1, 1)},
                new int[]{0}, new int[]{1, 2}, new int[]{3, 4, 5}, new int[]{6}, new int[]{7, 8, 9, 10});
        RelaxedClockMlEstimator.Estimate heterogeneous = fit(length,
                new LocalDate[]{LocalDate.of(2015, 1, 1), LocalDate.of(2018, 1, 1), LocalDate.of(2020, 1, 1), LocalDate.of(2023, 1, 1)},
                refX, xP, xY, yQ, yR);

        System.out.println("[RelaxedClockMlEstimatorTest] heterogeneous: rate0=" + heterogeneous.rate0SubPerSiteYear()
                + " sigma=" + heterogeneous.sigma() + " CV=" + heterogeneous.coefficientOfVariation()
                + " (homogeneous baseline CV=" + homogeneous.coefficientOfVariation() + ")");

        assertThat(heterogeneous.coefficientOfVariation()).isGreaterThan(0.15);
        assertThat(heterogeneous.coefficientOfVariation()).isGreaterThan(homogeneous.coefficientOfVariation() * 3);
        // The fitted typical rate should still land in the right neighbourhood despite one outlier
        // branch -- a lognormal MLE's mean-of-logs is far more robust to one extreme edge than a
        // strict clock's every-edge-weighted-equally least-squares fit would be.
        assertThat(heterogeneous.rate0SubPerSiteYear()).isBetween(3e-6, 3e-5);
    }

    private RelaxedClockMlEstimator.Estimate fit(int length, LocalDate[] dates, int[] refX, int[] xP, int[] xY, int[] yQ, int[] yR) {
        char[] base = new char[length];
        java.util.Arrays.fill(base, 'A');

        NucleotideSequence ref = new NucleotideSequence("ref", new String(base)).withMetadata(dates[0], null, null);
        NucleotideSequence p = new NucleotideSequence("P", mutate(base, refX, xP)).withMetadata(dates[1], null, null);
        NucleotideSequence q = new NucleotideSequence("Q", mutate(base, refX, xY, yQ)).withMetadata(dates[2], null, null);
        NucleotideSequence r = new NucleotideSequence("R", mutate(base, refX, xY, yR)).withMetadata(dates[3], null, null);

        SequenceAlignment alignment = SequenceAlignment.of(List.of(ref, p, q, r), "ref");
        PhyloTreeFactory.BuiltTree built = PhyloTreeFactory.build(alignment, GdMethod.JUKES_CANTOR);

        Map<String, Double> tipDates = new HashMap<>();
        for (NucleotideSequence seq : List.of(ref, p, q, r)) {
            tipDates.put(seq.getId(), TemporalUtil.toDecimalYear(seq.getCollectionDate().get()));
        }
        return new RelaxedClockMlEstimator().estimate(built.tree(), tipDates);
    }

    private static String mutate(char[] base, int[]... blocks) {
        char[] copy = base.clone();
        for (int[] block : blocks) {
            for (int pos : block) copy[pos] = 'G';
        }
        return new String(copy);
    }
}
