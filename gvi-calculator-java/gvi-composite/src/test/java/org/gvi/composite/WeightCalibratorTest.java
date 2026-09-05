package org.gvi.composite;

import org.gvi.core.exception.GviInputException;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.EnumMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.assertThatThrownBy;
import static org.assertj.core.api.Assertions.within;

class WeightCalibratorTest {

    private final WeightCalibrator calibrator = new WeightCalibrator();

    @Test
    void recoversKnownWeightsFromNoiseFreeSyntheticData() {
        List<IndexKey> keys = List.of(IndexKey.MU, IndexKey.RE, IndexKey.PI);
        double[] trueWeights = {0.6, 0.3, 0.1};

        Random rng = new Random(42);
        List<WeightCalibrator.Observation> observations = new ArrayList<>();
        for (int i = 0; i < 60; i++) {
            Map<IndexKey, Double> values = new EnumMap<>(IndexKey.class);
            double mu = rng.nextDouble(), re = rng.nextDouble(), pi = rng.nextDouble();
            values.put(IndexKey.MU, mu);
            values.put(IndexKey.RE, re);
            values.put(IndexKey.PI, pi);
            double target = trueWeights[0] * mu + trueWeights[1] * re + trueWeights[2] * pi;
            observations.add(new WeightCalibrator.Observation(values, target));
        }

        WeightCalibrator.CalibrationResult result = calibrator.calibrate(observations, keys);

        assertThat(result.weights().get(IndexKey.MU)).isCloseTo(0.6, within(0.02));
        assertThat(result.weights().get(IndexKey.RE)).isCloseTo(0.3, within(0.02));
        assertThat(result.weights().get(IndexKey.PI)).isCloseTo(0.1, within(0.02));
        assertThat(result.trainRmse()).isLessThan(0.01); // near-exact fit on noise-free, correctly-specified data
        assertThat(result.weights().get(IndexKey.MU) + result.weights().get(IndexKey.RE) + result.weights().get(IndexKey.PI))
                .isCloseTo(1.0, within(1e-6));
    }

    @Test
    void holdoutRmseIsAlsoLowWhenTheModelIsCorrectlySpecified() {
        List<IndexKey> keys = List.of(IndexKey.MB, IndexKey.GD);
        double[] trueWeights = {0.7, 0.3};

        Random rng = new Random(7);
        List<WeightCalibrator.Observation> observations = new ArrayList<>();
        for (int i = 0; i < 40; i++) {
            Map<IndexKey, Double> values = new EnumMap<>(IndexKey.class);
            double mb = rng.nextDouble(), gd = rng.nextDouble();
            values.put(IndexKey.MB, mb);
            values.put(IndexKey.GD, gd);
            observations.add(new WeightCalibrator.Observation(values, trueWeights[0] * mb + trueWeights[1] * gd));
        }

        WeightCalibrator.CalibrationResult result = calibrator.calibrateWithHoldout(observations, keys, 0.2, 123L);

        assertThat(result.holdoutRmse()).isNotNull();
        assertThat(result.holdoutRmse()).isLessThan(0.02);
        assertThat(result.observationsUsed()).isEqualTo(32); // 40 - 8 held out (20%)
    }

    @Test
    void rejectsTooFewObservations() {
        assertThatThrownBy(() -> calibrator.calibrate(List.of(), List.of(IndexKey.MU)))
                .isInstanceOf(GviInputException.class);
    }

    @Test
    void rejectsEmptyKeyList() {
        Map<IndexKey, Double> values = new EnumMap<>(IndexKey.class);
        values.put(IndexKey.MU, 0.5);
        List<WeightCalibrator.Observation> obs = List.of(
                new WeightCalibrator.Observation(values, 0.5),
                new WeightCalibrator.Observation(values, 0.5)
        );
        assertThatThrownBy(() -> calibrator.calibrate(obs, List.of())).isInstanceOf(GviInputException.class);
    }
}
