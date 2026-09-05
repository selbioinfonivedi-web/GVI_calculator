package org.gvi.algorithms.re;

import org.gvi.core.model.ConfidenceInterval;

import java.time.LocalDate;

/** One point of the Re(t) time series (Cori et al. sliding-window estimate). */
public record ReDailyEstimate(LocalDate date, double mean, ConfidenceInterval credibleInterval) {
}
