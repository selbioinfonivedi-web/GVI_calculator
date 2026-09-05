package org.gvi.core.io;

/**
 * A single non-fatal problem encountered while parsing an input file.
 * The parser records these and keeps going instead of aborting the whole
 * batch for one bad record (see build spec Section 2 / 4.3).
 */
public record ParseWarning(String source, long lineNumber, String recordId, String reason) {

    @Override
    public String toString() {
        String loc = lineNumber > 0 ? (":" + lineNumber) : "";
        String id = (recordId == null || recordId.isBlank()) ? "" : (" [" + recordId + "]") ;
        return source + loc + id + " - " + reason;
    }
}
