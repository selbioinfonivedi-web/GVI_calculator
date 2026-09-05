package org.gvi.core.io;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Accumulates {@link ParseWarning}s across a parse run and the counts of
 * accepted vs. skipped records, so every report can tell the user exactly
 * what was excluded and why -- never a silent drop.
 */
public final class ParseReport {

    private final List<ParseWarning> warnings = new ArrayList<>();
    private int acceptedCount = 0;
    private int skippedCount = 0;

    public void accept() {
        acceptedCount++;
    }

    public void skip(ParseWarning warning) {
        skippedCount++;
        warnings.add(warning);
    }

    /** A caveat about an accepted record (e.g. reduced date precision) -- unlike {@link #skip}, the record is still kept. */
    public void note(ParseWarning warning) {
        warnings.add(warning);
    }

    public List<ParseWarning> getWarnings() {
        return Collections.unmodifiableList(warnings);
    }

    public int getAcceptedCount() {
        return acceptedCount;
    }

    public int getSkippedCount() {
        return skippedCount;
    }

    public boolean hasWarnings() {
        return !warnings.isEmpty();
    }

    public String summary() {
        return String.format("accepted=%d, skipped=%d, warnings=%d", acceptedCount, skippedCount, warnings.size());
    }
}
