package org.gvi.core.io;

import org.gvi.core.exception.GviInputException;
import org.gvi.core.model.IncidencePoint;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import java.io.IOException;
import java.io.StringReader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.time.LocalDate;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.assertThatThrownBy;

/**
 * The incidence CSV is the input for the Cori et al. estimator -- the one Re path that recovers a
 * known value to within 0.02, against tree-shape inference which does far worse. It was at 0%
 * coverage, and no corpus dataset supplies such a file, so the accurate path was entirely
 * unexercised from the outside.
 * <p>
 * Row-level tolerance is the delicate part: a bad row is skipped rather than aborting the series,
 * which is right for a surveillance export with the odd blank cell, but it means a file that is
 * mostly unparsable would otherwise produce an Re from whatever few rows survived. Every skip
 * therefore has to reach the report.
 */
class IncidenceCsvReaderTest {

    @TempDir
    Path tempDir;

    private static IncidenceCsvReader.Result parse(String csv) {
        return IncidenceCsvReader.read(new StringReader(csv), "test.csv");
    }

    @Test
    void aWellFormedSeriesIsReadInFull() {
        var result = parse("""
                date,new_cases
                2021-01-01,5
                2021-01-02,8
                2021-01-03,13
                """);

        assertThat(result.points()).hasSize(3);
        assertThat(result.points()).extracting(IncidencePoint::newCases).containsExactly(5.0, 8.0, 13.0);
        assertThat(result.points().get(0).date()).isEqualTo(LocalDate.of(2021, 1, 1));
        assertThat(result.report().getAcceptedCount()).isEqualTo(3);
        assertThat(result.report().hasWarnings()).isFalse();
    }

    /**
     * Cori's estimator convolves the series against the serial-interval distribution by position, so
     * an out-of-order export would silently mis-align the whole convolution.
     */
    @Test
    void rowsAreSortedByDateRegardlessOfFileOrder() {
        var result = parse("""
                date,new_cases
                2021-03-01,30
                2021-01-01,10
                2021-02-01,20
                """);

        assertThat(result.points()).extracting(IncidencePoint::date)
                .containsExactly(LocalDate.of(2021, 1, 1), LocalDate.of(2021, 2, 1), LocalDate.of(2021, 3, 1));
        assertThat(result.points()).extracting(IncidencePoint::newCases).containsExactly(10.0, 20.0, 30.0);
    }

    @Test
    void fractionalCountsAreAcceptedBecauseSmoothedSurveillanceDataIsNotIntegral() {
        assertThat(parse("date,new_cases\n2021-01-01,5.5\n").points().get(0).newCases()).isEqualTo(5.5);
    }

    @Test
    void whitespaceAroundValuesIsTolerated() {
        var result = parse("date,new_cases\n  2021-01-01  ,  7  \n");
        assertThat(result.points().get(0).newCases()).isEqualTo(7.0);
        assertThat(result.points().get(0).date()).isEqualTo(LocalDate.of(2021, 1, 1));
    }

    @Test
    void aZeroCaseDayIsRealDataAndIsKept() {
        var result = parse("date,new_cases\n2021-01-01,0\n2021-01-02,4\n");
        assertThat(result.points()).hasSize(2);
        assertThat(result.points().get(0).newCases()).isZero();
    }

    // ── rejections that must be loud ─────────────────────────────────────────
    @Test
    void missingColumnsAreRejectedByName() {
        assertThatThrownBy(() -> parse("date,cases\n2021-01-01,5\n"))
                .isInstanceOf(GviInputException.class)
                .hasMessageContaining("test.csv")
                .hasMessageContaining("new_cases");

        assertThatThrownBy(() -> parse("day,new_cases\n2021-01-01,5\n"))
                .isInstanceOf(GviInputException.class)
                .hasMessageContaining("date");
    }

    /** An Re fitted to zero observations would be meaningless, so this must fail rather than return. */
    @Test
    void aFileWhoseEveryRowIsUnusableIsRejectedRatherThanReturningAnEmptySeries() {
        assertThatThrownBy(() -> parse("date,new_cases\nnot-a-date,5\n2021-01-02,not-a-number\n"))
                .isInstanceOf(GviInputException.class)
                .hasMessageContaining("no usable rows");
    }

    @Test
    void aHeaderWithNoRowsIsRejected() {
        assertThatThrownBy(() -> parse("date,new_cases\n"))
                .isInstanceOf(GviInputException.class)
                .hasMessageContaining("no usable rows");
    }

    // ── skips that must be reported, not silent ──────────────────────────────
    @Test
    void anUnparsableRowIsSkippedAndReportedWithItsLineNumber() {
        var result = parse("""
                date,new_cases
                2021-01-01,5
                garbage,row
                2021-01-03,13
                """);

        assertThat(result.points()).hasSize(2);
        assertThat(result.report().getSkippedCount()).isEqualTo(1);
        assertThat(result.report().getWarnings()).anyMatch(w -> w.reason().contains("Unparsable row"));
    }

    /**
     * A negative case count is a data-entry error, not a real observation, and would drag the
     * estimated Re downward if summed into the renewal equation.
     */
    @Test
    void aNegativeCountIsSkippedAndSaysSo() {
        var result = parse("""
                date,new_cases
                2021-01-01,5
                2021-01-02,-3
                2021-01-03,13
                """);

        assertThat(result.points()).hasSize(2);
        assertThat(result.report().getSkippedCount()).isEqualTo(1);
        assertThat(result.report().getWarnings())
                .anyMatch(w -> w.reason().contains("Negative new_cases")
                        && "2021-01-02".equals(w.recordId()));
    }

    @Test
    void theAcceptedAndSkippedCountsTogetherAccountForEveryDataRow() {
        var result = parse("""
                date,new_cases
                2021-01-01,5
                2021-01-02,-3
                bad,bad
                2021-01-04,9
                """);

        assertThat(result.report().getAcceptedCount() + result.report().getSkippedCount()).isEqualTo(4);
        assertThat(result.report().summary()).contains("accepted=2").contains("skipped=2");
    }

    // ── the Path overload ────────────────────────────────────────────────────
    @Test
    void readingFromDiskMatchesReadingFromAReader() throws IOException {
        Path csv = tempDir.resolve("inc.csv");
        Files.writeString(csv, "date,new_cases\n2021-01-01,5\n2021-01-02,8\n");

        var result = IncidenceCsvReader.read(csv);
        assertThat(result.points()).hasSize(2);
        assertThat(result.points().get(1).newCases()).isEqualTo(8.0);
    }

    @Test
    void anAbsentFileIsAnInputErrorNamingThePath() {
        Path missing = tempDir.resolve("nope.csv");
        assertThatThrownBy(() -> IncidenceCsvReader.read(missing))
                .isInstanceOf(GviInputException.class)
                .hasMessageContaining("Cannot read incidence CSV")
                .hasMessageContaining("nope.csv");
    }
}
