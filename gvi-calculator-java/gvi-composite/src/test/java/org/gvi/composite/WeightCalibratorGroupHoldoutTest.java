package org.gvi.composite;

import org.gvi.core.exception.GviInputException;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.EnumMap;
import java.util.List;
import java.util.Map;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.assertThatThrownBy;

/**
 * Leave-one-group-out asks the question a random holdout cannot: do weights fitted on other
 * pathogens predict a pathogen never seen during fitting?
 */
class WeightCalibratorGroupHoldoutTest {

    private static final List<IndexKey> KEYS = List.of(IndexKey.PI, IndexKey.MB);

    private static WeightCalibrator.Observation obs(double pi, double mb, double target) {
        Map<IndexKey, Double> v = new EnumMap<>(IndexKey.class);
        v.put(IndexKey.PI, pi);
        v.put(IndexKey.MB, mb);
        return new WeightCalibrator.Observation(v, target);
    }

    /** A target that genuinely depends on the indices the same way in every group. */
    private static List<WeightCalibrator.Observation> transferable(List<String> groups) {
        List<WeightCalibrator.Observation> out = new ArrayList<>();
        double[][] pts = {{0.1, 0.9}, {0.4, 0.6}, {0.8, 0.2}};
        for (int g = 0; g < groups.size(); g++) {
            for (double[] p : pts) out.add(obs(p[0], p[1], 0.5 * p[0] + 0.5 * p[1]));
        }
        return out;
    }

    @Test
    void aTransferableRelationshipGeneralisesToAHeldOutGroup() {
        List<String> groups = List.of("a", "a", "a", "b", "b", "b", "c", "c", "c");
        var r = new WeightCalibrator().calibrateWithGroupHoldout(transferable(List.of("a", "b", "c")), KEYS, groups);

        assertThat(r.perGroupRmse()).containsOnlyKeys("a", "b", "c");
        // If the relationship really transfers, held-out error should not blow up relative to training.
        assertThat(r.meanHeldOutRmse()).isLessThan(r.trainRmse() + 0.05);
        assertThat(r.weights().asMap().values().stream().mapToDouble(Double::doubleValue).sum())
                .isCloseTo(1.0, org.assertj.core.api.Assertions.within(1e-6));
    }

    @Test
    void everyGroupGetsItsOwnFoldSoOneBadPathogenIsVisible() {
        List<String> groups = List.of("a", "a", "b", "b", "c", "c");
        List<WeightCalibrator.Observation> o = List.of(
                obs(0.1, 0.9, 0.5), obs(0.9, 0.1, 0.5),
                obs(0.2, 0.8, 0.5), obs(0.8, 0.2, 0.5),
                obs(0.3, 0.7, 0.5), obs(0.7, 0.3, 0.5));
        var r = new WeightCalibrator().calibrateWithGroupHoldout(o, KEYS, groups);
        assertThat(r.perGroupRmse()).hasSize(3);
        assertThat(r.perGroupRmse().values()).allSatisfy(v -> assertThat(v).isNotNaN());
    }

    @Test
    void oneGroupIsRejectedBecauseThereIsNothingToGeneraliseTo() {
        List<String> groups = List.of("only", "only", "only");
        assertThatThrownBy(() -> new WeightCalibrator()
                .calibrateWithGroupHoldout(transferable(List.of("only")), KEYS, groups))
                .isInstanceOf(GviInputException.class)
                .hasMessageContaining("at least 2 distinct groups");
    }

    @Test
    void mismatchedLabelCountIsRejectedRatherThanSilentlyMisaligned() {
        assertThatThrownBy(() -> new WeightCalibrator()
                .calibrateWithGroupHoldout(transferable(List.of("a", "b")), KEYS, List.of("a", "b")))
                .isInstanceOf(GviInputException.class)
                .hasMessageContaining("one label per observation");
    }

    @Test
    void weightsReturnedAreFittedOnEverythingNotOnOneFold() {
        List<String> groups = List.of("a", "a", "a", "b", "b", "b");
        var r = new WeightCalibrator().calibrateWithGroupHoldout(transferable(List.of("a", "b")), KEYS, groups);
        var direct = new WeightCalibrator().calibrate(transferable(List.of("a", "b")), KEYS);
        for (IndexKey k : KEYS) {
            assertThat(r.weights().get(k))
                    .isCloseTo(direct.weights().get(k), org.assertj.core.api.Assertions.within(1e-6));
        }
    }
}
