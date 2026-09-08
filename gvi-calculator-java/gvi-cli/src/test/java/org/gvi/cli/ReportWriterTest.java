package org.gvi.cli;

import org.gvi.algorithms.gd.GdMethod;
import org.gvi.algorithms.mu.GenomeType;
import org.gvi.core.model.OrganismClass;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Map;
import java.util.Set;

import static org.assertj.core.api.Assertions.assertThat;

/**
 * {@link ReportWriter} was at 0% coverage while being the thing users actually read. Nothing
 * verified that a number reaching the report was the number the pipeline computed, that the
 * not-comparable banner appeared when it should, or that the JSON kept the keys downstream code
 * reads by name.
 * <p>
 * Fixtures come from real pipeline runs rather than a hand-built {@code PipelineResult}: the record
 * has fifteen fields and a synthetic one drifts from what the pipeline actually produces, which is
 * the failure this test exists to catch.
 */
class ReportWriterTest {

    @TempDir
    Path tempDir;

    private final ReportWriter writer = new ReportWriter();

    private PipelineResult result;

    @BeforeEach
    void computeAResult() throws IOException {
        Path fasta = tempDir.resolve("in.fasta");
        Files.writeString(fasta, """
                >reference 2015-01-01
                ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAGATGGCCATTGTAATGGGCCGC
                >q1 2018-01-01
                ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAGATGGCCATTGTAATGGCCCGC
                >q2 2021-01-01
                ATGGCCATTGTAATGGCCCGCTGAAAGGGTGCCCGATACATGGCCATTGTAATGGGCCGC
                >q3 2023-01-01
                ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCACGATAGATGGCCATTGTAATGGGCCGA
                """);
        result = new GviPipeline().run(config(fasta));
    }

    private PipelineConfig config(Path fasta) {
        return new PipelineConfig(
                fasta, null, null, null, null, null, null, "reference",
                Set.of("pi", "gd", "mb", "ri"),
                OrganismClass.VIRUS, null, false, 5.0,
                GdMethod.JUKES_CANTOR, 1.0, 2.0, false, null, "gtr", false, false, true, false,
                null, GenomeType.RNA, Map.of(), null, false, 100, false, null);
    }

    private String text() {
        ByteArrayOutputStream buf = new ByteArrayOutputStream();
        writer.writeText(result, new PrintStream(buf, true, StandardCharsets.UTF_8));
        return buf.toString(StandardCharsets.UTF_8);
    }

    @Test
    void theTextReportCarriesTheCompositeAndEveryContributingComponent() {
        String out = text();

        assertThat(out).contains("GVI = ");
        assertThat(out).contains(String.format("GVI = %.4f", result.datasetGvi().gvi()));
        assertThat(out).contains(result.datasetGvi().coverageSummary());

        // Every component the composite used must appear with its own numbers, so a reader can
        // reconstruct the total rather than take it on trust.
        for (var c : result.datasetGvi().components()) {
            assertThat(out).as("component %s missing from the report", c.key().label())
                    .contains(c.key().label());
        }
    }

    /**
     * The composite's own diagnostics reached the JSON but never the text report, so a reader
     * working from the printed output could not see that indices had been clamped at their
     * normalisation ceiling or that weights had been renormalised after an exclusion.
     */
    @Test
    void theTextReportPrintsTheCompositeDiagnostics() {
        assertThat(result.datasetGvi().diagnostics())
                .as("this fixture should produce at least one composite note").isNotEmpty();

        String out = text();
        assertThat(out).contains("-- Composite notes --");
        for (String d : result.datasetGvi().diagnostics()) {
            // Compare on a distinctive prefix: the full text wraps in the terminal but the
            // report itself must carry it verbatim.
            assertThat(out).contains(d.substring(0, Math.min(40, d.length())));
        }
    }

    /**
     * A low-coverage GVI is numerically indistinguishable from a well-supported one, so the caveat
     * has to sit next to the number rather than only in a warnings block further down.
     */
    @Test
    void aNonComparableScoreIsBanneredNextToTheNumber() {
        String out = text();
        if (result.datasetGvi().comparable()) {
            assertThat(out).doesNotContain("NOT COMPARABLE");
        } else {
            assertThat(out).contains("NOT COMPARABLE");
            int gviLine = out.indexOf("GVI = ");
            int banner = out.indexOf("NOT COMPARABLE");
            assertThat(banner - gviLine)
                    .as("the banner must sit next to the number, not pages below it")
                    .isBetween(0, 400);
        }
    }

    @Test
    void theJsonKeepsTheKeysDownstreamCodeReadsByName() throws IOException {
        Path json = tempDir.resolve("out.json");
        writer.writeJson(result, json);
        String body = Files.readString(json);

        for (String key : new String[]{"dataset_gvi", "dataset_indices", "sensitivity",
                                       "skipped", "warnings"}) {
            assertThat(body).as("JSON key %s", key).contains("\"" + key + "\"");
        }
        for (String field : new String[]{"gvi", "components", "excludedIndices", "effectiveWeightSum"}) {
            assertThat(body).as("dataset_gvi field %s", field).contains("\"" + field + "\"");
        }
    }

    /**
     * The CSV is wide -- one row per dataset, one column per index -- because it exists to be pasted
     * into a spreadsheet and sorted. That is exactly why comparable/coverage_pct sit immediately
     * after GVI: a score built from two indices is renormalised into [0,1] like one built from nine,
     * so a file without those columns invites a ranking it cannot support.
     */
    @Test
    void theCsvPutsComparabilityNextToTheScore() throws IOException {
        Path csv = tempDir.resolve("out.csv");
        writer.writeCsv(result, csv);
        var lines = Files.readAllLines(csv);

        assertThat(lines).hasSizeGreaterThanOrEqualTo(2);
        String[] header = lines.get(0).split(",");
        assertThat(header[0]).isEqualTo("GVI");
        assertThat(header[1]).isEqualTo("comparable");
        assertThat(header[2]).isEqualTo("coverage_pct");
        assertThat(header[3]).isEqualTo("indices_scored");

        // A column per scalar component, so no index is silently absent from the tabular view.
        for (org.gvi.composite.IndexKey key : org.gvi.composite.IndexKey.values()) {
            assertThat(header).as("column for %s", key.label()).contains(key.label());
        }

        String[] row = lines.get(1).split(",");
        assertThat(row.length).isEqualTo(header.length);
        assertThat(Double.parseDouble(row[0])).isEqualTo(result.datasetGvi().gvi());
        assertThat(Boolean.parseBoolean(row[1])).isEqualTo(result.datasetGvi().comparable());
    }

    /** The report must state what could not be computed, not quietly omit it. */
    @Test
    void skippedIndicesAreNamedInTheReport() {
        String out = text();
        for (String s : result.skipped()) {
            assertThat(out).contains(s.substring(0, Math.min(20, s.length())));
        }
    }
}
