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

    /**
     * A pathogen whose generation-time entry is still a stub must skip Re and let every other
     * index score, not abort the run.
     * <p>
     * It used to abort. Resolution happened in {@code GviCli} before the config was built, so
     * {@code --pathogen-id fmd} exited with code 2 and produced no GVI at all -- discarding the
     * eight indices that need no generation time because of one that does. Every other gate in
     * this pipeline excludes its own index and continues; a missing generation time is absent
     * data, not a caller error, and is now disposed of the same way.
     */
    @Test
    void aStubGenerationTimeSkipsReRatherThanAbortingTheRun() throws IOException {
        Path fasta = tempDir.resolve("samples.fasta");
        Files.writeString(fasta, """
                >reference 2015-01-01
                ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG
                >query1 2018-01-01
                ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG
                >query2 2021-01-01
                ATGGCCATTGTAATGGCCCGCTGAAAGGGTGCCCGATAC
                >query3 2023-01-01
                ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCACGATAG
                """);

        // Specifically an UNKNOWN entry -- a generation time nobody has filled in yet. The table
        // also carries FALSE entries, where Re is not merely unmeasured but meaningless: anthrax
        // is acquired from environmental spores rather than from a preceding case, so it has no
        // serial interval at all. Both must skip Re rather than abort, but they are different
        // findings and the pipeline words them differently.
        String stubId = org.gvi.algorithms.re.GenerationTimeTable.bundled().allEntries().stream()
                .filter(e -> e.reApplicable() == org.gvi.algorithms.re.GenerationTimeTable.Applicability.UNKNOWN)
                .map(org.gvi.algorithms.re.GenerationTimeTable.Entry::pathogenId)
                .findFirst().orElseThrow();

        PipelineConfig config = new PipelineConfig(
                fasta, null, null, null, null, null, null, "reference",
                Set.of("pi", "gd", "mb", "re"),
                org.gvi.core.model.OrganismClass.VIRUS, stubId, false, 5.0,
                GdMethod.JUKES_CANTOR, 1.0, 2.0, false, null, "gtr", false, false, false, true, false,
                null, GenomeType.RNA, Map.of(), null, false, 100, false, null);

        PipelineResult result = new GviPipeline().run(config);

        assertThat(result.datasetGvi()).as("the run must still produce a composite").isNotNull();
        assertThat(result.skipped())
                .as("Re should be skipped, naming the stub entry as the reason")
                .anyMatch(sk -> sk.startsWith("Re:") && sk.contains("stub"));
        assertThat(result.datasetGvi().components())
                .as("the indices that need no generation time must still score")
                .isNotEmpty();
    }

    /**
     * The companion case: a pathogen for which Re is not merely unmeasured but meaningless. The
     * table marks these {@code re_applicable: false} -- anthrax is acquired from environmental
     * spores, not from a preceding case, so there is no transmission chain and no serial interval.
     * The run must still complete and say so, rather than reporting a placeholder Re.
     */
    @Test
    void aPathogenWithNoTransmissionChainSkipsReWithThatReason() throws IOException {
        Path fasta = tempDir.resolve("nochain.fasta");
        Files.writeString(fasta, """
                >reference 2015-01-01
                ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG
                >query1 2018-01-01
                ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG
                >query2 2021-01-01
                ATGGCCATTGTAATGGCCCGCTGAAAGGGTGCCCGATAC
                >query3 2023-01-01
                ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCACGATAG
                """);

        String notApplicableId = org.gvi.algorithms.re.GenerationTimeTable.bundled().allEntries().stream()
                .filter(e -> e.reApplicable() == org.gvi.algorithms.re.GenerationTimeTable.Applicability.FALSE)
                .map(org.gvi.algorithms.re.GenerationTimeTable.Entry::pathogenId)
                .findFirst().orElseThrow();

        PipelineConfig config = new PipelineConfig(
                fasta, null, null, null, null, null, null, "reference",
                Set.of("pi", "gd", "mb", "re"),
                org.gvi.core.model.OrganismClass.BACTERIUM, notApplicableId, false, 5.0,
                GdMethod.JUKES_CANTOR, 1.0, 2.0, false, null, "gtr", false, false, false, true, false,
                null, GenomeType.DNA, Map.of(), null, false, 100, false, null);

        PipelineResult result = new GviPipeline().run(config);

        assertThat(result.datasetGvi()).isNotNull();
        assertThat(result.skipped())
                .anyMatch(sk -> sk.startsWith("Re:") && sk.contains("not applicable"));
    }
}
