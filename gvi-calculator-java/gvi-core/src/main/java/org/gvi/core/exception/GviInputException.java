package org.gvi.core.exception;

/**
 * Raised for malformed, missing, or inconsistent input data (bad FASTA,
 * missing metadata column, unparsable date, mismatched alignment length).
 * These are user-fixable data problems, not bugs -- the message must tell
 * the user exactly what to fix (file, line, field).
 */
public class GviInputException extends GviException {

    public GviInputException(String message) {
        super(message);
    }

    public GviInputException(String message, Throwable cause) {
        super(message, cause);
    }
}
