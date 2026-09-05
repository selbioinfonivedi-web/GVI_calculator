package org.gvi.composite;

import org.junit.jupiter.api.Test;

import java.util.List;

import static org.assertj.core.api.Assertions.assertThat;

/**
 * The composite engine renormalizes whichever indices survived quality gating so their weights sum to 1,
 * which means a GVI built from two indices lands in [0,1] and is numerically indistinguishable from one
 * built from nine. {@code effectiveWeightSum} recorded how much of the scheme actually had data, but
 * nothing consumed it -- so low-coverage scores were published into the same comparison column as
 * well-supported ones.
 */
class GviResultCoverageTest {

    private static GviComponent component(IndexKey key, double effectiveWeight, double normalized) {
        return new GviComponent(key, normalized, normalized, effectiveWeight, effectiveWeight, effectiveWeight * normalized);
    }

    /**
     * The real enterotoxaemia case: seven of nine indices failed quality gating, leaving RI and GC with
     * 5.7% of the scheme between them. RI was 0.0 and held 56% of the renormalized weight, so the reported
     * GVI was 0.053 -- which next to haemorrhagic septicaemia's 0.680 reads as "much lower risk" when the
     * truth is that the data was too broken to score. A low GVI from low coverage is the most dangerous
     * output shape here, because by value alone it is a confident finding of low virulence.
     */
    @Test
    void flagsALowCoverageScoreAsNotComparable() {
        GviResult result = new GviResult(0.0529,
                List.of(component(IndexKey.RI, 0.56, 0.0), component(IndexKey.GC, 0.44, 0.1203)),
                List.of(IndexKey.MU, IndexKey.RE, IndexKey.PI, IndexKey.MB, IndexKey.DNDS, IndexKey.GD, IndexKey.CAI),
                0.0573, List.of());

        assertThat(result.comparable()).isFalse();
        assertThat(result.coverageSummary()).isEqualTo("2 of 9 indices, 6% of the weighting scheme");
    }

    @Test
    void treatsAWellPopulatedScoreAsComparable() {
        GviResult result = new GviResult(0.4010,
                List.of(component(IndexKey.MU, 0.128, 1.0), component(IndexKey.RE, 0.284, 0.359),
                        component(IndexKey.MB, 0.204, 0.011), component(IndexKey.GD, 0.095, 0.049)),
                List.of(IndexKey.RI),
                0.93, List.of());

        assertThat(result.comparable()).isTrue();
        assertThat(result.coverageSummary()).isEqualTo("4 of 5 indices, 93% of the weighting scheme");
    }

    /** The threshold is inclusive, so a score sitting exactly on it is reported as comparable. */
    @Test
    void treatsCoverageExactlyAtTheThresholdAsComparable() {
        GviResult result = new GviResult(0.3,
                List.of(component(IndexKey.MU, 1.0, 0.3)), List.of(),
                GviResult.MIN_COMPARABLE_COVERAGE, List.of());

        assertThat(result.comparable()).isTrue();
    }
}
