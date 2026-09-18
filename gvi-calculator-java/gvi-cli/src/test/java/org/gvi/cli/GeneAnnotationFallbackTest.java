package org.gvi.cli;

import org.gvi.algorithms.gd.GdMethod;
import org.gvi.algorithms.mu.GenomeType;
import org.gvi.core.model.OrganismClass;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Map;
import java.util.Set;

import static org.assertj.core.api.Assertions.assertThat;

/**
 * When gene coordinates cannot be loaded, the pipeline guesses them with a 6-frame ORF scan and
 * still reports a dN/dS. That is a defensible fallback -- but only while the reason it happened is
 * stated accurately, because the user's next action depends on it.
 * <p>
 * It was not accurate. {@code predictOrfs} is reached three ways -- no {@code --gff}, an unreadable
 * one, or one that parsed to no CDS record -- and hardcoded "No --gff supplied" for all three. A
 * GFF3 saved with spaces instead of tabs (every row fails the nine-column check) therefore produced
 * a guessed dN/dS accompanied by an instruction to supply the flag the user had just supplied,
 * pointing away from the actual fault in their file.
 */
class GeneAnnotationFallbackTest {

    @TempDir
    Path tempDir;

    private static final String FASTA = """
            >reference 2015-01-01
            ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAGATGGCCATTGTAATGGGCCGCATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG
            >q1 2018-01-01
            ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAGATGGCCATTGTAATGGCCCGCATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG
            >q2 2021-01-01
            ATGGCCATTGTAATGGCCCGCTGAAAGGGTGCCCGATACATGGCCATTGTAATGGGCCGCATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG
            """;

    private Path fasta() throws IOException {
        Path p = tempDir.resolve("in.fasta");
        Files.writeString(p, FASTA);
        return p;
    }

    private PipelineResult run(Path gff) throws IOException {
        return new GviPipeline().run(new PipelineConfig(
                fasta(), null, gff, null, null, null, null, "reference",
                Set.of("dnds", "pi"),
                OrganismClass.VIRUS, null, true, 5.0,
                GdMethod.JUKES_CANTOR, 1.0, 2.0, false, null, "gtr", false, false, false, false, false,
                null, GenomeType.RNA, Map.of(), null, false, 100, false, null));
    }

    private static String warningsOf(PipelineResult r) {
        return String.join("\n", r.warnings());
    }

    @Test
    void withNoGffAtAllTheWarningSaysSo() throws IOException {
        String w = warningsOf(run(null));
        assertThat(w).contains("No --gff supplied");
    }

    /**
     * The case that motivated this: a space-separated GFF3. It looks correct to a human, parses to
     * zero CDS rows, and must not be reported as an absent annotation.
     */
    @Test
    void aSpaceSeparatedGffIsNotReportedAsAMissingGff() throws IOException {
        Path gff = tempDir.resolve("spaces.gff3");
        Files.writeString(gff, "##gff-version 3\n"
                + "ref manual CDS 1 60 . + 0 ID=g1;gene=VP1\n");

        String w = warningsOf(run(gff));

        assertThat(w).as("the file was supplied; saying otherwise sends the user to the wrong fix")
                .doesNotContain("No --gff supplied");
        assertThat(w).contains("The supplied --gff parsed but contained no usable CDS record");
        assertThat(w).as("and it should name the actual likely cause")
                .contains("spaces instead of tabs");
    }

    @Test
    void aGffWithOnlyNonCdsFeaturesIsAlsoReportedAsSuppliedButUnusable() throws IOException {
        Path gff = tempDir.resolve("nocds.gff3");
        Files.writeString(gff, "##gff-version 3\n"
                + "ref\tmanual\tgene\t1\t60\t.\t+\t.\tID=g1\n");

        String w = warningsOf(run(gff));
        assertThat(w).doesNotContain("No --gff supplied");
        assertThat(w).contains("contained no usable CDS record");
    }

    @Test
    void anAbsentGffPathIsReportedAsUnreadableRatherThanAsNotSupplied() throws IOException {
        String w = warningsOf(run(tempDir.resolve("missing.gff3")));

        assertThat(w).doesNotContain("No --gff supplied");
        assertThat(w).containsAnyOf("could not be read", "Cannot read GFF3 file");
    }

    /**
     * The dN/dS premise advisory carried the same false claim, and it is the one attached directly
     * to the number a reader is deciding whether to trust.
     */
    @Test
    void theDnDsPremiseAdvisoryDoesNotClaimNoGffWasSuppliedWhenOneWas() throws IOException {
        Path gff = tempDir.resolve("spaces.gff3");
        Files.writeString(gff, "ref manual CDS 1 60 . + 0 ID=g1\n");

        String w = warningsOf(run(gff));
        if (w.contains("dN/dS/CAI premise")) {
            assertThat(w).doesNotContain("no --gff was supplied");
            assertThat(w).contains("no usable gene annotation was available");
        }
    }

    /** A good GFF3 must still be used, and none of these warnings raised. */
    @Test
    void aValidGffIsUsedAndRaisesNoFallbackWarning() throws IOException {
        Path gff = tempDir.resolve("good.gff3");
        Files.writeString(gff, "##gff-version 3\n"
                + "ref\tmanual\tCDS\t1\t60\t.\t+\t0\tID=g1;gene=VP1\n");

        String w = warningsOf(run(gff));
        assertThat(w).doesNotContain("No --gff supplied");
        assertThat(w).doesNotContain("contained no usable CDS record");
        assertThat(w).doesNotContain("could not be read");
    }
}
