package org.gvi.web;

import org.gvi.algorithms.gd.GdMethod;
import org.gvi.algorithms.mu.GenomeType;
import org.gvi.cli.GviPipeline;
import org.gvi.cli.PipelineConfig;
import org.gvi.cli.PipelineResult;
import org.gvi.core.exception.GviInputException;
import org.gvi.core.model.OrganismClass;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;
import java.util.Set;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.assertThatThrownBy;
import static org.assertj.core.api.Assertions.within;

class AnalysisServiceTest {

    private final AnalysisService service = new AnalysisService();

    /**
     * Eight sequences with a handful of substitutions -- enough for the position-wise indices to run
     * without tripping the gap-fraction gate.
     */
    private static String alignment() {
        String base = "ATGGCACGTTCAGGCACCTTAGCACGTGGCACCTTAGCACGTTCAGGCACCTTAGCACGTTCAGGCACCTTAGCACGTTCAGGCTAA";
        StringBuilder fasta = new StringBuilder();
        for (int i = 0; i < 8; i++) {
            StringBuilder seq = new StringBuilder(base);
            // Vary a few third-codon positions per sequence so there is real, frame-aware divergence.
            for (int m = 0; m <= i; m++) {
                int pos = 5 + m * 9;
                if (pos < seq.length()) seq.setCharAt(pos, seq.charAt(pos) == 'A' ? 'G' : 'A');
            }
            fasta.append(">seq").append(i).append('\n').append(seq).append('\n');
        }
        return fasta.toString();
    }

    /**
     * The web front end must be a different presentation of the same computation, not a second
     * implementation of it. The JavaFX UI regressed exactly here: it called the back-compatible
     * {@link PipelineConfig} constructor, which defaulted organism class and genome type to UNSPECIFIED,
     * so identical input scored differently depending on which front end ran it. This pins the web
     * service's output against a direct pipeline run configured the same way.
     */
    @Test
    void producesTheSameCompositeScoreAsADirectPipelineRun(@TempDir Path tmp) throws Exception {
        AnalyzeRequest request = new AnalyzeRequest();
        request.fasta = alignment();
        request.organismClass = "virus";
        request.genomeType = "rna";

        AnalyzeResponse web = service.analyze(request);

        Path fasta = Files.writeString(tmp.resolve("in.fasta"), alignment());
        PipelineConfig config = new PipelineConfig(
                fasta, null, null, null, null, null, null, null, Set.of("all"),
                OrganismClass.VIRUS, null, false, 5.0, GdMethod.JUKES_CANTOR, 1.0, 2.0,
                false, null, "gtr", false, false, true, false, null,
                GenomeType.RNA, Map.of(), null, false, 100, false, null);
        PipelineResult direct = new GviPipeline().run(config);

        assertThat(web.gvi).isNotNull();
        assertThat(direct.datasetGvi()).isNotNull();
        assertThat(web.gvi).isCloseTo(direct.datasetGvi().gvi(), within(1e-12));
        assertThat(web.comparable).isEqualTo(direct.datasetGvi().comparable());
        assertThat(web.coverageSummary).isEqualTo(direct.datasetGvi().coverageSummary());
    }

    /**
     * An index that quality gating removed must still appear in the response. The exclusion reason
     * refers to that index's value, so hiding the value would leave the explanation pointing at a
     * number the analyst cannot see.
     */
    @Test
    void reportsExcludedIndicesAlongsideTheOnesThatScored() throws Exception {
        AnalyzeRequest request = new AnalyzeRequest();
        request.fasta = alignment();
        request.organismClass = "virus";
        request.genomeType = "rna";

        AnalyzeResponse response = service.analyze(request);

        assertThat(response.indices).isNotEmpty();
        for (AnalyzeResponse.Exclusion exclusion : response.exclusions) {
            assertThat(response.indices).containsKey(exclusion.key);
            assertThat(response.indices.get(exclusion.key).scored).isFalse();
            assertThat(exclusion.reason).isNotBlank();
        }
        List<String> scoredKeys = response.components.stream().map(c -> c.key).toList();
        for (String key : scoredKeys) {
            assertThat(response.indices.get(key).scored).isTrue();
        }
    }

    @Test
    void rejectsAMissingAlignmentWithAnActionableMessage() {
        AnalyzeRequest request = new AnalyzeRequest();
        request.fasta = "   ";

        assertThatThrownBy(() -> service.analyze(request))
                .isInstanceOf(GviInputException.class)
                .hasMessageContaining("aligned multi-FASTA");
    }

    @Test
    void rejectsInputOverTheSizeCapRatherThanAttemptingIt() {
        AnalyzeRequest request = new AnalyzeRequest();
        request.fasta = ">x\n" + "A".repeat(AnalysisService.MAX_INPUT_CHARS + 1);

        assertThatThrownBy(() -> service.analyze(request))
                .isInstanceOf(GviInputException.class)
                .hasMessageContaining("over the");
    }

    @Test
    void rejectsAnUnknownGenomeTypeInsteadOfSilentlyFallingBack() {
        AnalyzeRequest request = new AnalyzeRequest();
        request.fasta = alignment();
        request.genomeType = "protein";

        assertThatThrownBy(() -> service.analyze(request))
                .isInstanceOf(GviInputException.class)
                .hasMessageContaining("expected rna or dna");
    }

    /** The temp directory each analysis stages its inputs into must not survive the call. */
    @Test
    void cleansUpItsTemporaryInputFiles() throws Exception {
        Path tmpRoot = Path.of(System.getProperty("java.io.tmpdir"));
        long before = countWorkDirs(tmpRoot);

        AnalyzeRequest request = new AnalyzeRequest();
        request.fasta = alignment();
        service.analyze(request);

        assertThat(countWorkDirs(tmpRoot)).isEqualTo(before);
    }

    private static long countWorkDirs(Path root) throws Exception {
        try (var paths = Files.list(root)) {
            return paths.filter(p -> p.getFileName().toString().startsWith("gvi-web-")).count();
        }
    }
}
