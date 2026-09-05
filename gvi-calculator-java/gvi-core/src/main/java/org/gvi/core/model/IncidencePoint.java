package org.gvi.core.model;

import java.time.LocalDate;

/** One (date, new case count) observation feeding the Cori et al. Re estimator (Section 5.2a). */
public record IncidencePoint(LocalDate date, double newCases) {

    public IncidencePoint {
        if (newCases < 0) {
            throw new IllegalArgumentException("newCases must be >= 0 on " + date + ", got " + newCases);
        }
    }
}
