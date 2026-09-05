package org.gvi.algorithms.mu;

import org.junit.jupiter.api.Test;

import static org.assertj.core.api.Assertions.assertThat;

/**
 * The reference table's terminal DNA band ("No known DNA virus replicates this fast") was computed and
 * printed while the rate it condemned still carried full composite weight -- and because mu normalizes
 * against a 1e-4 DNA ceiling, the more impossible the rate the harder it pushed the score up.
 */
class MuPlausibilityTest {

    /** The real Lumpy Skin Disease rate: 1.85e-4 for a capripoxvirus the table puts around 1e-6. */
    @Test
    void flagsTheLumpySkinDiseaseRateAsImplausibleForADnaGenome() {
        assertThat(MuReferenceTable.isImplausible(1.845e-4, GenomeType.DNA)).isTrue();
        assertThat(MuReferenceTable.classify(1.845e-4, GenomeType.DNA).polymeraseType())
                .isEqualTo("Implausible for a DNA genome");
    }

    @Test
    void acceptsARateInsideTheFastestCredibleDnaBand() {
        // Small ssDNA (circovirus/parvovirus) territory -- fast for DNA, but real.
        assertThat(MuReferenceTable.isImplausible(9e-6, GenomeType.DNA)).isFalse();
    }

    @Test
    void acceptsATypicalCapripoxvirusRate() {
        assertThat(MuReferenceTable.isImplausible(1e-6, GenomeType.DNA)).isFalse();
    }

    /**
     * RNA genomes have no impossibility ceiling here: the top RNA band (VSV, "ultra high error") is an
     * observed regime, not an artifact, so an RNA rate must never be gated by this check.
     */
    @Test
    void neverFlagsAnRnaGenomeHoweverFast() {
        assertThat(MuReferenceTable.isImplausible(1e-2, GenomeType.RNA)).isFalse();
        assertThat(MuReferenceTable.isImplausible(0.2, GenomeType.RNA)).isFalse();
    }

    /** Without a declared genome type there is no scale to judge against, so the check stays out of the way. */
    @Test
    void doesNotGuessWhenGenomeTypeIsUnspecified() {
        assertThat(MuReferenceTable.isImplausible(1.845e-4, GenomeType.UNSPECIFIED)).isFalse();
    }
}
