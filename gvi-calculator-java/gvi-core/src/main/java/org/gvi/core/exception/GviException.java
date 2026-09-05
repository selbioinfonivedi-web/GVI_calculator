package org.gvi.core.exception;

/**
 * Base type for every error the GVI application raises deliberately.
 * Never let a raw {@link RuntimeException} (NPE, IndexOutOfBounds, etc.)
 * escape a module boundary uncaught -- wrap it in one of these subtypes
 * with an actionable message instead.
 */
public class GviException extends RuntimeException {

    public GviException(String message) {
        super(message);
    }

    public GviException(String message, Throwable cause) {
        super(message, cause);
    }
}
