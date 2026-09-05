package org.gvi.cli;

import org.gvi.algorithms.gd.GdMethod;
import org.gvi.algorithms.mu.GenomeType;
import org.gvi.composite.IndexKey;
import org.gvi.composite.SensitivityResult;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Map;
import java.util.Set;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.assertThatThrownBy;

class GviPipelineTest {

    @TempDir
    Path tempDir;

    /**
     * Whole-file GVI on indices that need nothing but the FASTA itself (pi,
     * GD, MB, RI -- no dates/GFF/codon-usage/incidence required) must also
     * carry a sensitivity-analysis breakdown for every index that
     * contributed a nonzero default weight (Section 5.9/8).
     */
    @Test
    void wholeFileGviCarriesSensitivityAnalysisForContributingIndices() throws IOException {
        Path fasta = tempDir.resolve("samples.fasta");
        Files.writeString(fasta, """
                >reference
                ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG
                >query1
                ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG
                >query2
                ATGGCCATTGTAATGGCCCGCTGAAAGGGTGCCCGATAC
                >query3
                ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCACGATAG
                """);

        PipelineConfig config = new PipelineConfig(
                fasta, null, null, null, null, null, null, "reference",
                Set.of("pi", "gd", "mb", "ri"),
                5.0, GdMethod.JUKES_CANTOR, 1.0, 2.0,
                false, null, "gtr", false, false, false, false, null, GenomeType.UNSPECIFIED, Map.of(), null);

        PipelineResult result = new GviPipeline().run(config);

        assertThat(result.datasetGvi()).isNotNull();
        assertThat(result.sensitivity()).isNotEmpty();
        Set<IndexKey> sensitivityKeys = result.sensitivity().stream().map(SensitivityResult::key).collect(java.util.stream.Collectors.toSet());
        Set<IndexKey> componentKeys = result.datasetGvi().components().stream().map(c -> c.key()).collect(java.util.stream.Collectors.toSet());
        assertThat(sensitivityKeys).isEqualTo(componentKeys);

        for (SensitivityResult s : result.sensitivity()) {
            assertThat(s.baseGvi()).isEqualTo(result.datasetGvi().gvi());
            assertThat(s.gviAtHighWeight()).isNotEqualTo(s.gviAtLowWeight());
        }
    }

    /**
     * {@link DatasetSummary} is plumbed straight out of the already-built {@code SequenceAlignment} --
     * no new computation, just surfacing facts (sequence count, alignment length, resolved reference,
     * collection-date span from FASTA-header dates) that {@link GviPipeline#run} already knew but
     * previously discarded.
     */
    @Test
    void wholeFileResultCarriesADatasetSummaryWithRealFacts() throws IOException {
        Path fasta = tempDir.resolve("dated.fasta");
        // Space before the year (not a pipe) so firstToken() -- which splits on whitespace only --
        // resolves the id to plain "reference"/"query1"/"query2"; FastaHeaderDateParser scans the
        // whole header line regardless, so the embedded year is still found.
        Files.writeString(fasta, """
                >reference 2015
                ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG
                >query1 2011
                ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG
                >query2 2013
                ATGGCCATTGTAATGGCCCGCTGAAAGGGTGCCCGATAC
                """);

        PipelineConfig config = new PipelineConfig(
                fasta, null, null, null, null, null, null, "reference",
                Set.of("pi"),
                5.0, GdMethod.JUKES_CANTOR, 1.0, 2.0,
                false, null, "gtr", false, false, false, false, null, GenomeType.UNSPECIFIED, Map.of(), null);

        PipelineResult result = new GviPipeline().run(config);

        assertThat(result.datasetSummary()).isNotNull();
        assertThat(result.datasetSummary().sequenceCount()).isEqualTo(3);
        assertThat(result.datasetSummary().alignmentLengthBp()).isEqualTo(39);
        assertThat(result.datasetSummary().referenceId()).isEqualTo("reference");
        assertThat(result.datasetSummary().earliestCollectionDate().getYear()).isEqualTo(2011);
        assertThat(result.datasetSummary().latestCollectionDate().getYear()).isEqualTo(2015);
    }

    /**
     * A sequence sampled LONGER after the reference that is nonetheless MORE similar to it than one
     * sampled shortly after (deliberately inverted here) drives a negative root-to-tip slope --
     * {@link GviPipeline}'s data-quality pass must flag this rather than silently reporting a negative
     * "rate," since a negative substitutions-per-site-per-year is never a real molecular-clock estimate.
     */
    @Test
    void negativeMuTriggersADataQualityAdvisory() throws IOException {
        Path fasta = tempDir.resolve("negative_mu.fasta");
        Files.writeString(fasta, """
                >reference 2000
                ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG
                >query1 2001
                %s
                >query2 2020
                ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG
                """.formatted("T".repeat(39)));

        PipelineConfig config = new PipelineConfig(
                fasta, null, null, null, null, null, null, "reference",
                Set.of("mu"),
                5.0, GdMethod.JUKES_CANTOR, 1.0, 2.0,
                false, null, "gtr", false, false, false, false, null, GenomeType.UNSPECIFIED, Map.of(), null);

        PipelineResult result = new GviPipeline().run(config);

        assertThat(result.warnings()).anyMatch(w -> w.contains("Data quality advisory -- mu"));
    }

    /** A gap-heavy "alignment" (padding instead of real alignment) must be flagged, not silently accepted. */
    @Test
    void highGapFractionAlignmentTriggersADataQualityAdvisory() throws IOException {
        Path fasta = tempDir.resolve("gappy.fasta");
        Files.writeString(fasta, """
                >reference
                ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG
                >query1
                ATGGCCATTGTAATGGGCC%s
                """.formatted("-".repeat(20)));

        PipelineConfig config = new PipelineConfig(
                fasta, null, null, null, null, null, null, "reference",
                Set.of("pi"),
                5.0, GdMethod.JUKES_CANTOR, 1.0, 2.0,
                false, null, "gtr", false, false, false, false, null, GenomeType.UNSPECIFIED, Map.of(), null);

        PipelineResult result = new GviPipeline().run(config);

        assertThat(result.warnings()).anyMatch(w -> w.contains("Data quality advisory -- alignment"));
    }

    /** Mutation burden far exceeding what's plausible for closely-related strains must be flagged. */
    @Test
    void highMutationBurdenFractionTriggersADataQualityAdvisory() throws IOException {
        Path fasta = tempDir.resolve("high_mb.fasta");
        Files.writeString(fasta, """
                >reference
                ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG
                >query1
                %sATGGGCCGCTGAAAGGGTGCCCGATAG
                """.formatted("T".repeat(12)));

        PipelineConfig config = new PipelineConfig(
                fasta, null, null, null, null, null, null, "reference",
                Set.of("mb"),
                5.0, GdMethod.JUKES_CANTOR, 1.0, 2.0,
                false, null, "gtr", false, false, false, false, null, GenomeType.UNSPECIFIED, Map.of(), null);

        PipelineResult result = new GviPipeline().run(config);

        assertThat(result.warnings()).anyMatch(w -> w.contains("Data quality advisory -- mutation burden"));
    }

    /**
     * A misspelled {@code --codon-usage-species} used to be caught and turned into a "CAI skipped"
     * note. That is the wrong disposal for a caller error: the composite renormalizes the surviving
     * weights and the run prints a confident GVI computed from one index fewer, so a typo produced a
     * plausible wrong number rather than a visible failure. It cost the corpus driver CAI on all 7
     * datasets that could compute it. The unknown name must abort the run.
     */
    @Test
    void unknownCodonUsageSpeciesAbortsRatherThanSilentlyDroppingCai() throws IOException {
        Path fasta = tempDir.resolve("samples.fasta");
        Files.writeString(fasta, """
                >reference
                ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG
                >query1
                ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG
                >query2
                ATGGCCATTGTAATGGCCCGCTGAAAGGGTGCCCGATAC
                """);

        PipelineConfig config = new PipelineConfig(
                fasta, null, null, null, "bos_taurus", null, null, "reference",
                Set.of("cai"),
                5.0, GdMethod.JUKES_CANTOR, 1.0, 2.0,
                false, null, "gtr", false, false, false, false, null, GenomeType.UNSPECIFIED, Map.of(), null);

        assertThatThrownBy(() -> new GviPipeline().run(config))
                .isInstanceOf(org.gvi.core.exception.GviInputException.class)
                .hasMessageContaining("bos_taurus")
                .hasMessageContaining("cattle")
                .hasMessageContaining("aborted");
    }

    /** The correctly-spelled neighbour of the name above must still work. */
    @Test
    void knownCodonUsageSpeciesStillComputesCai() throws IOException {
        Path fasta = tempDir.resolve("samples.fasta");
        Files.writeString(fasta, """
                >reference
                ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG
                >query1
                ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG
                >query2
                ATGGCCATTGTAATGGCCCGCTGAAAGGGTGCCCGATAC
                """);

        PipelineConfig config = new PipelineConfig(
                fasta, null, null, null, "cattle", null, null, "reference",
                Set.of("cai"),
                5.0, GdMethod.JUKES_CANTOR, 1.0, 2.0,
                false, null, "gtr", false, false, false, false, null, GenomeType.UNSPECIFIED, Map.of(), null);

        PipelineResult result = new GviPipeline().run(config);
        assertThat(result.skipped()).noneMatch(s -> s.startsWith("CAI:"));
    }
}
