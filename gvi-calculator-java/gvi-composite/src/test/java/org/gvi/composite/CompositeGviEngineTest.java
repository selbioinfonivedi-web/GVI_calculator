package org.gvi.composite;

import org.gvi.core.exception.GviComputationException;
import org.junit.jupiter.api.Test;

import java.util.EnumMap;
import java.util.Map;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.assertThatThrownBy;
import static org.assertj.core.api.Assertions.within;

class CompositeGviEngineTest {

    private final CompositeGviEngine engine = new CompositeGviEngine();

    @Test
    void defaultWeightsSumToOne() {
        CompositeWeights w = CompositeWeights.defaults();
        double sum = 0;
        for (IndexKey k : IndexKey.values()) sum += w.get(k);
        assertThat(sum).isCloseTo(1.0, within(1e-9));
    }

    @Test
    void rejectsWeightsNotSummingToOne() {
        Map<IndexKey, Double> bad = new EnumMap<>(IndexKey.class);
        bad.put(IndexKey.MU, 0.5);
        bad.put(IndexKey.RE, 0.6);
        assertThatThrownBy(() -> CompositeWeights.of(bad)).isInstanceOf(org.gvi.core.exception.GviInputException.class);
    }

    @Test
    void computesWeightedNormalizedSumWithAllIndicesPresent() {
        // Two indices, equal weight 0.5 each, ranges [0,10]. raw=5 -> normalized=0.5 both -> gvi=0.5
        Map<IndexKey, Double> weightMap = new EnumMap<>(IndexKey.class);
        weightMap.put(IndexKey.MU, 0.5);
        weightMap.put(IndexKey.RE, 0.5);
        CompositeWeights weights = CompositeWeights.of(weightMap);

        Map<IndexKey, NormalizationRange> ranges = new EnumMap<>(IndexKey.class);
        ranges.put(IndexKey.MU, new NormalizationRange(0, 10));
        ranges.put(IndexKey.RE, new NormalizationRange(0, 10));

        Map<IndexKey, FakeIndexResult> available = new EnumMap<>(IndexKey.class);
        available.put(IndexKey.MU, new FakeIndexResult("mu", 5.0));
        available.put(IndexKey.RE, new FakeIndexResult("Re", 5.0));

        GviResult r = engine.compute(available, weights, ranges);
        assertThat(r.gvi()).isCloseTo(0.5, within(1e-9));
        assertThat(r.excludedIndices()).isEmpty();
    }

    @Test
    void renormalizesWeightsWhenAnIndexIsMissing() {
        Map<IndexKey, Double> weightMap = new EnumMap<>(IndexKey.class);
        weightMap.put(IndexKey.MU, 0.3);
        weightMap.put(IndexKey.RE, 0.7);
        CompositeWeights weights = CompositeWeights.of(weightMap);

        Map<IndexKey, NormalizationRange> ranges = new EnumMap<>(IndexKey.class);
        ranges.put(IndexKey.MU, new NormalizationRange(0, 10));
        ranges.put(IndexKey.RE, new NormalizationRange(0, 10));

        // only MU available; RE missing entirely
        Map<IndexKey, FakeIndexResult> available = new EnumMap<>(IndexKey.class);
        available.put(IndexKey.MU, new FakeIndexResult("mu", 10.0)); // normalized = 1.0

        GviResult r = engine.compute(available, weights, ranges);
        // MU's weight should be renormalized to 1.0 (it's the only available index), so gvi = 1.0 * 1.0
        assertThat(r.gvi()).isCloseTo(1.0, within(1e-9));
        assertThat(r.excludedIndices()).containsExactly(IndexKey.RE);
        assertThat(r.components().get(0).effectiveWeight()).isCloseTo(1.0, within(1e-9));
    }

    @Test
    void throwsWhenNoWeightedIndexIsAvailable() {
        Map<IndexKey, Double> weightMap = new EnumMap<>(IndexKey.class);
        weightMap.put(IndexKey.MU, 1.0);
        CompositeWeights weights = CompositeWeights.of(weightMap);
        assertThatThrownBy(() -> engine.compute(Map.of(), weights)).isInstanceOf(GviComputationException.class);
    }

    @Test
    void clampsOutOfRangeValuesInsteadOfExceedingZeroOneBounds() {
        Map<IndexKey, Double> weightMap = new EnumMap<>(IndexKey.class);
        weightMap.put(IndexKey.MU, 1.0);
        CompositeWeights weights = CompositeWeights.of(weightMap);
        Map<IndexKey, NormalizationRange> ranges = Map.of(IndexKey.MU, new NormalizationRange(0, 10));

        Map<IndexKey, FakeIndexResult> available = Map.of(IndexKey.MU, new FakeIndexResult("mu", 999.0));
        GviResult r = engine.compute(available, weights, ranges);
        assertThat(r.gvi()).isCloseTo(1.0, within(1e-9));
        assertThat(r.components().get(0).normalizedValue()).isCloseTo(1.0, within(1e-9));
    }

    @Test
    void betaGenomicMatchesSpecExample() {
        // "If GVI = 0.5 and scale_factor = 2.0: beta increases by 100%"
        GviResult r = new GviResult(0.5, java.util.List.of(), java.util.List.of(), 1.0, java.util.List.of());
        assertThat(r.betaGenomic(100.0, 2.0)).isCloseTo(200.0, within(1e-9));
    }

    @Test
    void sensitivityAnalysisIncreasesGviWhenIncreasingWeightOfAHighValueIndex() {
        Map<IndexKey, Double> weightMap = new EnumMap<>(IndexKey.class);
        weightMap.put(IndexKey.MU, 0.5);
        weightMap.put(IndexKey.RE, 0.5);
        CompositeWeights weights = CompositeWeights.of(weightMap);
        Map<IndexKey, NormalizationRange> ranges = new EnumMap<>(IndexKey.class);
        ranges.put(IndexKey.MU, new NormalizationRange(0, 10));
        ranges.put(IndexKey.RE, new NormalizationRange(0, 10));

        Map<IndexKey, FakeIndexResult> available = new EnumMap<>(IndexKey.class);
        available.put(IndexKey.MU, new FakeIndexResult("mu", 10.0)); // normalized = 1.0
        available.put(IndexKey.RE, new FakeIndexResult("Re", 0.0));  // normalized = 0.0

        var results = engine.sensitivityAnalysis(available, weights, ranges, 0.2);
        var muSensitivity = results.stream().filter(s -> s.key() == IndexKey.MU).findFirst().orElseThrow();
        // increasing MU's weight (which has normalized value 1.0) should raise GVI
        assertThat(muSensitivity.gviAtHighWeight()).isGreaterThan(muSensitivity.gviAtLowWeight());
    }

    /**
     * An index clamped at its normalisation ceiling must say so. Clamped and genuinely-maximal both
     * read 1.0 and contribute identically, but a clamped index carries no discriminating information:
     * it would report the same value for a dataset twice as divergent. On this project's corpus pi
     * clamps on 6 of 13 scored datasets and GD on 2, so the reader needs telling.
     */
    @Test
    void anIndexClampedAtItsCeilingIsReportedWithHowFarOverItWent() {
        Map<IndexKey, Double> weightMap = new EnumMap<>(IndexKey.class);
        weightMap.put(IndexKey.PI, 0.5);
        weightMap.put(IndexKey.MU, 0.5);
        CompositeWeights weights = CompositeWeights.of(weightMap);

        Map<IndexKey, NormalizationRange> ranges = new EnumMap<>(IndexKey.class);
        ranges.put(IndexKey.PI, new NormalizationRange(0, 0.02));
        ranges.put(IndexKey.MU, new NormalizationRange(0, 10));

        Map<IndexKey, FakeIndexResult> available = new EnumMap<>(IndexKey.class);
        available.put(IndexKey.PI, new FakeIndexResult("pi", 0.0825));   // 4.1x the ceiling
        available.put(IndexKey.MU, new FakeIndexResult("mu", 5.0));      // mid-range

        GviResult r = engine.compute(available, weights, ranges);

        assertThat(r.diagnostics())
                .anyMatch(d -> d.contains("Normalisation ceiling reached by") && d.contains("pi"));
        assertThat(r.diagnostics()).anyMatch(d -> d.contains("4.1x over"));
        assertThat(r.diagnostics())
                .as("an index inside its range must not be named")
                .noneMatch(d -> d.contains("Normalisation ceiling") && d.contains("mu ("));
    }

    /**
     * A value at the FLOOR is a measurement, not a clamp. RI = 0.0 means the PHI test found no
     * recombination signal, which is a finding. Reporting it as saturation both cried wolf and
     * produced nonsense arithmetic ("0.0 against a 1.0 ceiling, 0.0x over").
     */
    @Test
    void aValueAtTheFloorIsNotReportedAsSaturated() {
        Map<IndexKey, Double> weightMap = new EnumMap<>(IndexKey.class);
        weightMap.put(IndexKey.RI, 1.0);
        CompositeWeights weights = CompositeWeights.of(weightMap);

        Map<IndexKey, NormalizationRange> ranges = new EnumMap<>(IndexKey.class);
        ranges.put(IndexKey.RI, new NormalizationRange(0, 1.0));

        Map<IndexKey, FakeIndexResult> available = new EnumMap<>(IndexKey.class);
        available.put(IndexKey.RI, new FakeIndexResult("RI", 0.0));

        GviResult r = engine.compute(available, weights, ranges);

        assertThat(r.diagnostics()).noneMatch(d -> d.contains("Normalisation ceiling"));
    }

    @Test
    void exceedsCeilingIsTrueOnlyAtOrAboveTheTop() {
        NormalizationRange range = new NormalizationRange(0, 0.02);
        assertThat(range.exceedsCeiling(0.0)).isFalse();
        assertThat(range.exceedsCeiling(0.019)).isFalse();
        assertThat(range.exceedsCeiling(0.02)).isTrue();
        assertThat(range.exceedsCeiling(0.0825)).isTrue();
    }
}
