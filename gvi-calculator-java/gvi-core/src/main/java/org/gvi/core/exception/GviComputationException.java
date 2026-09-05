package org.gvi.core.exception;

/**
 * Raised when a formula cannot be evaluated meaningfully for the given data
 * (e.g. Jukes-Cantor distance saturated, regression has fewer than 2 distinct
 * timepoints, insufficient informative sites for the PHI test). This signals
 * "the math breaks down here", as opposed to {@link GviInputException} which
 * signals "the input file itself is wrong".
 */
public class GviComputationException extends GviException {

    public GviComputationException(String message) {
        super(message);
    }

    public GviComputationException(String message, Throwable cause) {
        super(message, cause);
    }
}
