package org.gvi.composite;

import org.junit.jupiter.api.Test;

import java.util.List;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.within;

/**
 * beta_genomic must not be driven by a score that already contains Re.
 * <p>
 * In any compartmental model {@code Re = beta * infectious_period}, so feeding a
 * Re-inclusive GVI into a beta multiplier makes the downstream model derive its own
 * Re from the beta that Re itself inflated -- the same signal counted twice.
 */
class BetaGenomicReExclusionTest {

    private static GviComponent component(IndexKey key, double normalized, double effectiveWeight) {
        return new GviComponent(key, normalized, normalized, effectiveWeight, effectiveWeight, effectiveWeight * normalized);
    }

    @Test
    void betaIgnoresTheReComponentAndRenormalizesTheRest() {
        // Re normalized high (0.9) at weight 0.3; two genomic indices at 0.2 each, weight 0.35 apiece.
        List<GviComponent> components = List.of(
                component(IndexKey.RE, 0.9, 0.30),
                component(IndexKey.PI, 0.2, 0.35),
                component(IndexKey.MB, 0.2, 0.35));
        double gvi = 0.30 * 0.9 + 0.35 * 0.2 + 0.35 * 0.2;
        GviResult r = new GviResult(gvi, components, List.of(), 1.0, List.of());

        assertThat(r.gvi()).isCloseTo(0.41, within(1e-9));
        // Without Re: the two remaining indices renormalize from 0.35 each to 0.5 each, both at 0.2.
        assertThat(r.gviExcludingRe()).isCloseTo(0.2, within(1e-9));
        assertThat(r.betaExcludesRe()).isTrue();

        // beta must be built from 0.2, not 0.41.
        assertThat(r.betaGenomic(100.0, 2.0)).isCloseTo(100.0 * (1 + 0.2 * 2.0), within(1e-9));
    }

    @Test
    void betaIsUnchangedWhenReNeverContributed() {
        List<GviComponent> components = List.of(
                component(IndexKey.PI, 0.5, 0.5),
                component(IndexKey.MB, 0.5, 0.5));
        GviResult r = new GviResult(0.5, components, List.of(IndexKey.RE), 1.0, List.of());

        assertThat(r.gviExcludingRe()).isCloseTo(0.5, within(1e-9));
        assertThat(r.betaExcludesRe()).isFalse();
        // Matches the spec example: GVI 0.5 at scale 2.0 doubles beta.
        assertThat(r.betaGenomic(100.0, 2.0)).isCloseTo(200.0, within(1e-9));
    }

    @Test
    void reAsTheOnlyContributorLeavesNoGenomicSignal() {
        List<GviComponent> components = List.of(component(IndexKey.RE, 0.8, 1.0));
        GviResult r = new GviResult(0.8, components, List.of(), 1.0, List.of());

        assertThat(r.gviExcludingRe()).isCloseTo(0.0, within(1e-9));
        assertThat(r.betaGenomic(100.0, 2.0)).isCloseTo(100.0, within(1e-9)); // beta0 unmodified
    }
}
