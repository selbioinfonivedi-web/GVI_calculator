package org.gvi.core.util;

import java.time.LocalDate;
import java.time.temporal.ChronoUnit;

/** Shared date <-> decimal-year conversion for temporal-signal regressions (mu, Re). */
public final class TemporalUtil {

    private TemporalUtil() {
    }

    /** e.g. 2024-07-02 -> ~2024.5 */
    public static double toDecimalYear(LocalDate date) {
        LocalDate startOfYear = LocalDate.of(date.getYear(), 1, 1);
        LocalDate startOfNextYear = LocalDate.of(date.getYear() + 1, 1, 1);
        long dayOfYear = ChronoUnit.DAYS.between(startOfYear, date);
        long daysInYear = ChronoUnit.DAYS.between(startOfYear, startOfNextYear);
        return date.getYear() + (double) dayOfYear / daysInYear;
    }

    public static double yearsBetween(LocalDate earlier, LocalDate later) {
        return toDecimalYear(later) - toDecimalYear(earlier);
    }

    /**
     * Best single-date estimate for a bare-year collection date (e.g. GenBank/GFF metadata that only
     * records "2019", common for bacterial and parasite deposits, or a FASTA header with only a
     * standalone year). Resolving to January 1st -- this project's previous behavior -- silently
     * introduces about half a year of systematic bias toward the start of the year into every downstream
     * root-to-tip regression; the middle of the year is the unbiased choice when nothing more specific is
     * known.
     */
    public static LocalDate midpointOfYear(int year) {
        LocalDate startOfYear = LocalDate.of(year, 1, 1);
        LocalDate startOfNextYear = LocalDate.of(year + 1, 1, 1);
        long daysInYear = ChronoUnit.DAYS.between(startOfYear, startOfNextYear);
        return startOfYear.plusDays(daysInYear / 2);
    }

    /** Same idea as {@link #midpointOfYear}, for a collection date known only to year-and-month precision (e.g. "2019-06"). */
    public static LocalDate midpointOfYearMonth(int year, int month) {
        LocalDate startOfMonth = LocalDate.of(year, month, 1);
        return startOfMonth.plusDays(startOfMonth.lengthOfMonth() / 2);
    }
}
