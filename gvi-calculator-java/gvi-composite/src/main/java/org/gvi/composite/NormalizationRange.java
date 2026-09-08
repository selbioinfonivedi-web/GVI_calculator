package org.gvi.composite;

import java.util.EnumMap;
import java.util.Map;

/**
 * Configurable [min,max] reference range each index's raw value is
 * min-max-normalized against before weighting (Section 5.9). Values are
 * clamped to [0,1] -- an outlier beyond the configured range saturates
 * rather than distorting the composite.
 */
public record NormalizationRange(double min, double max) {

    public NormalizationRange {
        if (max <= min) {
            throw new IllegalArgumentException("max (" + max + ") must be > min (" + min + ")");
        }
    }

    /**
     * True when {@code raw} sits at or beyond the top of the range, so normalisation clamped it there.
     * <p>
     * Only the ceiling counts. A value at the floor of a natural domain is usually a real measurement
     * rather than lost information -- RI = 0.0 means the PHI test found no recombination signal,
     * which is a finding, not a clamp. Reporting it as saturation both cried wolf and produced
     * nonsense arithmetic ("0.0 against a 1.0 ceiling, 0.0x over").
     * <p>
     * A clamped index is indistinguishable, by its normalised value alone, from one that genuinely
     * sits at the extreme of its reference range -- both read 1.0 (or 0.0) and contribute exactly the
     * same. That matters because a clamped index carries no discriminating information: it would
     * report the same value for this dataset and for one twice as divergent. On this project's
     * corpus pi clamps on 6 of 13 scored datasets and GD on 2, so it is the common case rather than
     * an edge one, and the reader should be told rather than left to infer it from a suspiciously
     * round 1.0000.
     */
    public boolean exceedsCeiling(double raw) {
        return raw >= max;
    }

    public double normalize(double raw) {
        double v = (raw - min) / (max - min);
        return Math.max(0.0, Math.min(1.0, v));
    }

    /**
     * Defaults derived from each index's reference-table ceiling (Section
     * 5.1-5.8 reference tables) -- e.g. mu's domain tops out where the
     * spec's highest reference band begins. Fully overridable per
     * deployment/pathogen via the CLI config.
     */
    /**
     * mu's ceiling for a DNA genome. Large dsDNA viruses and bacteria evolve orders of
     * magnitude slower than RNA viruses, so measuring them against the RNA-scaled 1e-2
     * ceiling collapsed every correctly-classified DNA rate to a near-zero contribution
     * regardless of how fast it was <em>for a DNA genome</em>. Mirrors the split that
     * MuReferenceTable already applies to the category labels.
     */
    public static final double MU_MAX_DNA = 1e-4;
    /** mu's ceiling for an RNA genome (and the historical default when genome type is unknown). */
    public static final double MU_MAX_RNA = 1e-2;

    public static Map<IndexKey, NormalizationRange> defaults() {
        return defaults(MU_MAX_RNA);
    }

    /**
     * Ranges with mu scaled for the supplied genome chemistry. Callers pass
     * {@link #MU_MAX_DNA} for a DNA genome and {@link #MU_MAX_RNA} otherwise;
     * every other index is unaffected by genome type.
     */
    public static Map<IndexKey, NormalizationRange> defaults(double muMax) {
        Map<IndexKey, NormalizationRange> m = new EnumMap<>(IndexKey.class);
        m.put(IndexKey.MU, new NormalizationRange(0.0, muMax));
        m.put(IndexKey.RE, new NormalizationRange(0.0, 3.0));
        m.put(IndexKey.PI, new NormalizationRange(0.0, 0.02));
        // MB is normalized as mutation events PER KILOBASE, not as a raw count. A raw count
        // is not comparable between a 500 bp gene and a 3 kb one, and the previous fixed
        // 0-100 raw-count ceiling was already being exceeded by real data (a mean burden of
        // 257 events clamped to 1.0, making any value from 100 upward indistinguishable).
        // 100 events/kb is 10% divergence -- a defensible ceiling for "as different as these
        // reference-based indices can meaningfully describe".
        m.put(IndexKey.MB, new NormalizationRange(0.0, 100.0));
        m.put(IndexKey.DNDS, new NormalizationRange(0.0, 3.0));
        m.put(IndexKey.GD, new NormalizationRange(0.0, 0.05));
        m.put(IndexKey.CAI, new NormalizationRange(0.0, 1.0));
        m.put(IndexKey.GC, new NormalizationRange(0.0, 10.0));
        // RI deliberately departs from the "top of the spec's highest band" convention used above.
        // The spec's RI bands (top band opening at 0.30) were written for "proportion of recombinant
        // isolates", where 0.30 genuinely is a high value. What the composite actually receives is
        // RecombinationIndexCalculator.computeFromAlignment's PHI signal-strength score, which is a
        // *relative* drop in local vs. global incompatibility and is already clamped to its own
        // natural [0,1] domain. Normalizing that against a 0.30 ceiling rescaled it by ~3.33x, so any
        // ri >= 0.30 saturated the composite contribution outright. Min-max against the quantity's own
        // full domain is the identity-preserving choice and is what is used here.
        m.put(IndexKey.RI, new NormalizationRange(0.0, 1.0));
        return m;
    }
}
