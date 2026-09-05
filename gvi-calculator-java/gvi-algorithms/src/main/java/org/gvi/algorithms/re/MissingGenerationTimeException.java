package org.gvi.algorithms.re;

import org.gvi.core.exception.GviInputException;

/**
 * Raised when a {@code --pathogen-id} has no usable generation time in the
 * bundled table -- either no entry at all, or an entry that is deliberately a
 * stub / marked not applicable.
 * <p>
 * This is a hard failure by design. The behaviour it replaces was a silent
 * fallback to a single 5-day default for every organism, which compressed
 * genuinely different growth rates into a Re band of 1.0007-1.0262 across
 * pathogens as different as a DNA arbovirus and a gut bacterium -- reading as
 * "everything is stable" when it actually meant "the generation time was never
 * supplied". Failing loudly is the point.
 */
public class MissingGenerationTimeException extends GviInputException {

    public MissingGenerationTimeException(String message) {
        super(message);
    }
}
