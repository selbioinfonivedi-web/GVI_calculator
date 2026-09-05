package org.gvi.algorithms.mu;

import org.gvi.core.spi.IndexResult;

import java.util.List;

/**
 * Result of Index 1 - Evolutionary Rate (Section 5.1), estimated by
 * root-to-tip regression. {@code rSquared} is the temporal-signal strength
 * diagnostic -- callers must treat {@code reliable() == false} results as
 * indicative only, per the build spec's explicit requirement not to present
 * a weak-signal estimate as precise.
 */
public record MuResult(double muSubPerSiteYear, double rSquared, int pointsUsed, double timeSpanYears,
                        MuEstimationMethod method, String category, List<String> diagnostics,
                        boolean temporalSignalConfirmed) implements IndexResult {

    private static final double RELIABLE_R_SQUARED_THRESHOLD = 0.3;

    /**
     * Back-compatible constructor for callers predating the date-randomization test. Treats the
     * temporal signal as confirmed, preserving the previous behaviour exactly.
     */
    public MuResult(double muSubPerSiteYear, double rSquared, int pointsUsed, double timeSpanYears,
                    MuEstimationMethod method, String category, List<String> diagnostics) {
        this(muSubPerSiteYear, rSquared, pointsUsed, timeSpanYears, method, category, diagnostics, true);
    }

    @Override
    public String indexName() {
        return "mu";
    }

    @Override
    public double primaryValue() {
        return muSubPerSiteYear;
    }

    public boolean reliable() {
        return rSquared >= RELIABLE_R_SQUARED_THRESHOLD;
    }
}
