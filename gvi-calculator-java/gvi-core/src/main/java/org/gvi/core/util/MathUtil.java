package org.gvi.core.util;

public final class MathUtil {

    private MathUtil() {
    }

    /** Collapses IEEE-754 negative zero (e.g. -0.75 * ln(1) == -0.0) to positive zero for display/comparison sanity. */
    public static double stripNegativeZero(double v) {
        return v == 0.0 ? 0.0 : v;
    }
}
