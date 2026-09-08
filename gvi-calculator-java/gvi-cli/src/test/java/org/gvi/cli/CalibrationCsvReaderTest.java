package org.gvi.cli;

import org.gvi.composite.IndexKey;
import org.gvi.core.exception.GviInputException;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import java.io.IOException;
import java.io.StringReader;
import java.nio.file.Files;
import java.nio.file.Path;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.assertThatThrownBy;

/**
 * The calibration reader is the one place a user's own numbers enter the weighting scheme, so every
 * way a file can be wrong has to produce a message that names the file and the problem rather than a
 * NumberFormatException or a silently short fit. It was at 0% coverage.
 */
class CalibrationCsvReaderTest {

    @TempDir
    Path tempDir;

    private static CalibrationCsvReader.Result read(String csv) {
        return CalibrationCsvReader.read(new StringReader(csv), "test.csv");
    }

    @Test
    void aWellFormedFileYieldsOneObservationPerRowAndTheColumnsInFileOrder() {
        var result = read("""
                mu,pi,target
                0.10,0.20,1.5
                0.30,0.40,2.5
                0.50,0.60,3.5
                """);

        assertThat(result.observations()).hasSize(3);
        assertThat(result.keysInPlay()).containsExactly(IndexKey.MU, IndexKey.PI);

        var first = result.observations().get(0);
        assertThat(first.observedTarget()).isEqualTo(1.5);
        assertThat(first.normalizedValues()).containsEntry(IndexKey.MU, 0.10).containsEntry(IndexKey.PI, 0.20);
    }

    /** Labels are the public contract; the reader must accept them as written in the reports. */
    @Test
    void columnLabelsMatchTheIndexLabelsIncludingTheAwkwardOnes() {
        var result = read("""
                dN/dS,GC_Deviation,target
                0.2,0.3,1.0
                """);
        assertThat(result.keysInPlay()).containsExactly(IndexKey.DNDS, IndexKey.GC);
    }

    @Test
    void headerMatchingIsCaseInsensitive() {
        assertThat(read("MU,TARGET\n0.1,1.0\n").keysInPlay()).containsExactly(IndexKey.MU);
    }

    /**
     * Without a target there is nothing to fit against, and the fit would otherwise silently use
     * whichever column happened to sort last.
     */
    @Test
    void aMissingTargetColumnIsRejectedByName() {
        assertThatThrownBy(() -> read("mu,pi\n0.1,0.2\n"))
                .isInstanceOf(GviInputException.class)
                .hasMessageContaining("test.csv")
                .hasMessageContaining("target");
    }

    /**
     * A typo'd column must not be dropped: silently fitting on eight of nine columns produces
     * weights that look fine and are wrong.
     */
    @Test
    void anUnrecognizedColumnIsRejectedAndNamed() {
        assertThatThrownBy(() -> read("mu,virulence,target\n0.1,0.2,1.0\n"))
                .isInstanceOf(GviInputException.class)
                .hasMessageContaining("unrecognized column 'virulence'");
    }

    @Test
    void aFileOfNothingButTargetHasNoIndicesToFit() {
        assertThatThrownBy(() -> read("target\n1.0\n"))
                .isInstanceOf(GviInputException.class)
                .hasMessageContaining("no index columns");
    }

    @Test
    void aHeaderWithNoDataRowsIsRejectedRatherThanFittedOnNothing() {
        assertThatThrownBy(() -> read("mu,pi,target\n"))
                .isInstanceOf(GviInputException.class)
                .hasMessageContaining("no data rows");
    }

    @Test
    void anUnreadablePathIsReportedAsSuchAndNotAsAParseFailure() {
        assertThatThrownBy(() -> CalibrationCsvReader.read(tempDir.resolve("absent.csv")))
                .isInstanceOf(GviInputException.class)
                .hasMessageContaining("Cannot read calibration CSV");
    }

    @Test
    void theSamePathReadFromDiskGivesTheSameResultAsFromAReader() throws IOException {
        Path csv = tempDir.resolve("calib.csv");
        Files.writeString(csv, "mu,pi,target\n0.1,0.2,1.0\n0.3,0.4,2.0\n");

        var fromDisk = CalibrationCsvReader.read(csv);
        assertThat(fromDisk.observations()).hasSize(2);
        assertThat(fromDisk.keysInPlay()).containsExactly(IndexKey.MU, IndexKey.PI);
    }
}
