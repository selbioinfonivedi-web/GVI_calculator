package org.gvi.core.io;

import org.gvi.core.exception.GviInputException;
import org.gvi.core.model.GeneAnnotation;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.StringReader;
import java.nio.file.Files;
import java.nio.file.Path;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.assertThatThrownBy;

/**
 * {@link GffReader} decides which stretches of the alignment are coding, and therefore what dN/dS
 * and CAI are computed over at all. It was at 0% coverage. The failure mode that matters is silent:
 * a file whose CDS rows are all skipped yields an empty gene list, and the pipeline then falls back
 * to treating the whole sequence as one ORF -- a different, usually worse, number reported without
 * comment. So these check not only that good rows parse, but that every skip is recorded.
 */
class GffReaderTest {

    @TempDir
    Path tempDir;

    private static GffReader.Result parse(String gff) throws IOException {
        return GffReader.read(new BufferedReader(new StringReader(gff)), "test.gff3");
    }

    /** Tabs, not spaces: GFF3 is tab-delimited and a space-separated file is a common hand-edit. */
    private static String row(String... cols) {
        return String.join("\t", cols);
    }

    @Test
    void aCdsRowBecomesAGeneWithOneBasedInclusiveCoordinates() throws IOException {
        var result = parse("##gff-version 3\n"
                + row("chr1", ".", "CDS", "10", "99", ".", "+", "0", "ID=cds1;gene=VP1") + "\n");

        assertThat(result.genes()).hasSize(1);
        GeneAnnotation g = result.genes().get(0);
        assertThat(g.geneName()).isEqualTo("VP1");
        assertThat(g.start()).isEqualTo(10);
        assertThat(g.end()).isEqualTo(99);
        assertThat(g.strand()).isEqualTo('+');
        assertThat(g.length()).as("GFF3 coordinates are inclusive at both ends").isEqualTo(90);
        assertThat(result.report().getAcceptedCount()).isEqualTo(1);
    }

    @Test
    void commentsBlankLinesAndNonCdsFeaturesAreIgnoredWithoutBeingCountedAsErrors() throws IOException {
        var result = parse("##gff-version 3\n"
                + "# a comment\n"
                + "\n"
                + row("chr1", ".", "gene", "1", "500", ".", "+", ".", "ID=g1") + "\n"
                + row("chr1", ".", "exon", "1", "500", ".", "+", ".", "ID=e1") + "\n"
                + row("chr1", ".", "CDS", "1", "300", ".", "+", "0", "ID=c1") + "\n");

        assertThat(result.genes()).hasSize(1);
        assertThat(result.report().getSkippedCount())
                .as("a gene or exon row is not an error, just not a CDS")
                .isZero();
    }

    @Test
    void theFeatureTypeIsMatchedCaseInsensitively() throws IOException {
        assertThat(parse(row("c", ".", "cds", "1", "9", ".", "+", "0", "ID=x")).genes()).hasSize(1);
        assertThat(parse(row("c", ".", "Cds", "1", "9", ".", "+", "0", "ID=x")).genes()).hasSize(1);
    }

    /** The name is what a gene-specific reference table is later looked up by, so precedence matters. */
    @Test
    void theGeneNameIsTakenFromGeneThenNameThenIdAndFallsBackToUnnamed() throws IOException {
        assertThat(parse(row("c", ".", "CDS", "1", "9", ".", "+", "0", "ID=i;Name=n;gene=g")).genes()
                .get(0).geneName()).isEqualTo("g");
        assertThat(parse(row("c", ".", "CDS", "1", "9", ".", "+", "0", "ID=i;Name=n")).genes()
                .get(0).geneName()).isEqualTo("n");
        assertThat(parse(row("c", ".", "CDS", "1", "9", ".", "+", "0", "ID=i")).genes()
                .get(0).geneName()).isEqualTo("i");
        assertThat(parse(row("c", ".", "CDS", "1", "9", ".", "+", "0", ".")).genes()
                .get(0).geneName()).isEqualTo("unnamed");
    }

    @Test
    void theMinusStrandIsPreservedAndAnUnknownStrandDefaultsToPlus() throws IOException {
        assertThat(parse(row("c", ".", "CDS", "1", "9", ".", "-", "0", "ID=x")).genes()
                .get(0).strand()).isEqualTo('-');
        assertThat(parse(row("c", ".", "CDS", "1", "9", ".", ".", "0", "ID=x")).genes()
                .get(0).strand()).isEqualTo('+');
        assertThat(parse(row("c", ".", "CDS", "1", "9", ".", "", "0", "ID=x")).genes()
                .get(0).strand()).isEqualTo('+');
    }

    // ── the silent-loss paths ────────────────────────────────────────────────
    /**
     * A space-separated file looks fine to a human and parses to one column here. Every row must be
     * reported, because the alternative is an empty gene list and a whole-sequence ORF fallback that
     * nothing announces.
     */
    @Test
    void aSpaceSeparatedFileIsReportedRowByRowRatherThanSilentlyYieldingNoGenes() throws IOException {
        var result = parse("chr1 . CDS 1 300 . + 0 ID=c1\nchr1 . CDS 400 700 . + 0 ID=c2\n");

        assertThat(result.genes()).isEmpty();
        assertThat(result.report().getSkippedCount()).isGreaterThanOrEqualTo(2);
        assertThat(result.report().getWarnings())
                .anyMatch(w -> w.reason().contains("Fewer than 9 GFF3 columns"));
    }

    @Test
    void aNonNumericCoordinateIsSkippedAndReportedWithItsLineNumber() throws IOException {
        var result = parse(row("c", ".", "CDS", "one", "300", ".", "+", "0", "ID=c1") + "\n"
                + row("c", ".", "CDS", "400", "700", ".", "+", "0", "ID=c2") + "\n");

        assertThat(result.genes()).hasSize(1);
        assertThat(result.report().getSkippedCount()).isEqualTo(1);
        assertThat(result.report().getWarnings()).anyMatch(w ->
                w.reason().contains("Malformed CDS record") && w.lineNumber() == 1);
    }

    /** end &lt; start is rejected by GeneAnnotation; the reader must catch it, not propagate it. */
    @Test
    void aReversedCoordinatePairIsSkippedNotThrown() throws IOException {
        var result = parse(row("c", ".", "CDS", "700", "400", ".", "+", "0", "ID=c1") + "\n");

        assertThat(result.genes()).isEmpty();
        assertThat(result.report().getWarnings()).anyMatch(w -> w.reason().contains("Malformed CDS record"));
    }

    /**
     * An annotation with no CDS at all is the case the pipeline silently falls back on, so the reader
     * states it explicitly and names the consequence.
     */
    @Test
    void aFileWithNoCdsRecordsSaysWhatWillHappenInstead() throws IOException {
        var result = parse("##gff-version 3\n" + row("c", ".", "gene", "1", "500", ".", "+", ".", "ID=g1") + "\n");

        assertThat(result.genes()).isEmpty();
        assertThat(result.report().getWarnings()).anyMatch(w ->
                w.reason().contains("No CDS features found") && w.reason().contains("whole sequence as one ORF"));
    }

    @Test
    void anEmptyFileIsReportedTheSameWayAsOneWithNoCds() throws IOException {
        assertThat(parse("").report().getWarnings())
                .anyMatch(w -> w.reason().contains("No CDS features found"));
    }

    // ── the Path overload ────────────────────────────────────────────────────
    @Test
    void readingFromDiskGivesTheSameAnswerAsReadingFromAReader() throws IOException {
        Path gff = tempDir.resolve("a.gff3");
        Files.writeString(gff, "##gff-version 3\n"
                + row("chr1", ".", "CDS", "10", "99", ".", "+", "0", "gene=VP1") + "\n");

        var result = GffReader.read(gff);
        assertThat(result.genes()).hasSize(1);
        assertThat(result.genes().get(0).geneName()).isEqualTo("VP1");
    }

    @Test
    void anAbsentFileIsAnInputErrorNamingThePath() {
        Path missing = tempDir.resolve("nope.gff3");
        assertThatThrownBy(() -> GffReader.read(missing))
                .isInstanceOf(GviInputException.class)
                .hasMessageContaining("Cannot read GFF3 file")
                .hasMessageContaining("nope.gff3");
    }

    @Test
    void multipleCdsRecordsAreKeptInFileOrder() throws IOException {
        var result = parse(row("c", ".", "CDS", "1", "99", ".", "+", "0", "gene=A") + "\n"
                + row("c", ".", "CDS", "100", "199", ".", "-", "0", "gene=B") + "\n"
                + row("c", ".", "CDS", "200", "299", ".", "+", "0", "gene=C") + "\n");

        assertThat(result.genes()).extracting(GeneAnnotation::geneName).containsExactly("A", "B", "C");
        assertThat(result.report().getAcceptedCount()).isEqualTo(3);
    }
}
