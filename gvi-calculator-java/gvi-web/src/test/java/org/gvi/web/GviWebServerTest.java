package org.gvi.web;

import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.net.URI;
import java.net.http.HttpClient;
import java.net.http.HttpRequest;
import java.net.http.HttpResponse;

import static org.assertj.core.api.Assertions.assertThat;

/** End-to-end checks over the real HTTP surface, on an ephemeral loopback port. */
class GviWebServerTest {

    private final ObjectMapper mapper = new ObjectMapper();
    private final HttpClient client = HttpClient.newHttpClient();
    private GviWebServer server;
    private String base;

    @BeforeEach
    void startServer() throws Exception {
        server = new GviWebServer();
        server.start("127.0.0.1", 0);
        base = "http://127.0.0.1:" + server.port();
    }

    @AfterEach
    void stopServer() {
        if (server != null) server.stop();
    }

    private HttpResponse<String> post(String path, String body) throws Exception {
        return client.send(HttpRequest.newBuilder(URI.create(base + path))
                .header("Content-Type", "application/json")
                .POST(HttpRequest.BodyPublishers.ofString(body)).build(), HttpResponse.BodyHandlers.ofString());
    }

    private HttpResponse<String> get(String path) throws Exception {
        return client.send(HttpRequest.newBuilder(URI.create(base + path)).GET().build(),
                HttpResponse.BodyHandlers.ofString());
    }

    @Test
    void servesTheApplicationPage() throws Exception {
        HttpResponse<String> response = get("/");

        assertThat(response.statusCode()).isEqualTo(200);
        assertThat(response.body()).contains("Genomic Virulence Index Calculator");
    }

    @Test
    void reportsHealth() throws Exception {
        assertThat(mapper.readTree(get("/api/health").body()).get("status").asText()).isEqualTo("ok");
    }

    /** Bad input is the user's problem to fix, so it must come back as a 400 with the reason, not a 500. */
    @Test
    void returnsBadRequestWithTheReasonWhenTheAlignmentIsMissing() throws Exception {
        HttpResponse<String> response = post("/api/analyze", "{\"fasta\":\"\"}");

        assertThat(response.statusCode()).isEqualTo(400);
        assertThat(mapper.readTree(response.body()).get("error").asText()).contains("aligned multi-FASTA");
    }

    @Test
    void returnsBadRequestForMalformedJson() throws Exception {
        HttpResponse<String> response = post("/api/analyze", "{not json");

        assertThat(response.statusCode()).isEqualTo(400);
        assertThat(mapper.readTree(response.body()).get("error").asText()).contains("not valid JSON");
    }

    @Test
    void computesAScoreOverHttp() throws Exception {
        StringBuilder fasta = new StringBuilder();
        String base = "ATGGCACGTTCAGGCACCTTAGCACGTGGCACCTTAGCACGTTCAGGCACCTTAGCACGTTCAGGCTAA";
        for (int i = 0; i < 6; i++) {
            StringBuilder seq = new StringBuilder(base);
            seq.setCharAt(5 + i * 6, 'T');
            fasta.append(">seq").append(i).append("\\n").append(seq).append("\\n");
        }
        String body = mapper.writeValueAsString(java.util.Map.of(
                "fasta", fasta.toString().replace("\\n", "\n"),
                "organismClass", "virus", "genomeType", "rna"));

        HttpResponse<String> response = post("/api/analyze", body);

        assertThat(response.statusCode()).isEqualTo(200);
        JsonNode json = mapper.readTree(response.body());
        assertThat(json.has("gvi")).isTrue();
        assertThat(json.has("comparable")).isTrue();
        assertThat(json.get("reportText").asText()).contains("Genomic Virulence Index");
    }

    @Test
    void rejectsPathTraversalOnTheStaticHandler() throws Exception {
        assertThat(get("/../pom.xml").statusCode()).isIn(400, 404);
    }

    @Test
    void rejectsTheWrongHttpMethod() throws Exception {
        assertThat(get("/api/analyze").statusCode()).isEqualTo(405);
    }

    /**
     * The codon-table dropdown is populated from this endpoint rather than a list hardcoded in the
     * page, so the browser cannot offer a name the loader will reject. A misspelled species is not a
     * typo the pipeline absorbs -- CAI is dropped, the surviving weights are renormalized, and a
     * confident GVI computed from one index fewer comes back. Keeping the two in sync by hand is
     * exactly how that bug reached the corpus driver, so the list must come from the loader's own map.
     */
    @Test
    void servesTheBundledCodonSpeciesSoTheUiCannotOfferAnUnknownOne() throws Exception {
        HttpResponse<String> res = get("/api/codon-species");
        assertThat(res.statusCode()).isEqualTo(200);
        JsonNode species = mapper.readTree(res.body()).get("species");
        assertThat(species.isArray()).isTrue();
        assertThat(species).isNotEmpty();

        java.util.Set<String> served = new java.util.LinkedHashSet<>();
        species.forEach(n -> served.add(n.asText()));
        assertThat(served).isEqualTo(new java.util.TreeSet<>(
                org.gvi.algorithms.cai.BundledCodonUsageTables.availableSpecies()));

        // Every served name must actually load, or the dropdown is offering a dead option.
        for (String name : served) {
            assertThat(org.gvi.algorithms.cai.BundledCodonUsageTables.load(name)).isNotNull();
        }
    }

    /** The front end fetches app.js and styles.css by absolute path; neither may 404. */
    @Test
    void servesTheFrontEndAssetsThePageReferences() throws Exception {
        for (String asset : new String[]{"/app.js", "/styles.css"}) {
            HttpResponse<String> res = get(asset);
            assertThat(res.statusCode()).as(asset).isEqualTo(200);
            assertThat(res.body()).as(asset).isNotEmpty();
        }
        String page = get("/").body();
        assertThat(page).contains("/app.js").contains("/styles.css");
    }

    /**
     * Every field the page renders from must be present in the response. These are read by name in
     * app.js; a rename on the Java side would blank a panel in the browser with no error anywhere.
     */
    @Test
    void analyzeResponseCarriesEveryFieldTheFrontEndRenders() throws Exception {
        String fasta = ">a\nATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG\n"
                     + ">b\nATGGCCATTGTAATGGCCCGCTGAAAGGGTGCCCGATAC\n"
                     + ">c\nATGGCCATTGTAATGGGCCGCTGAAAGGGTGCACGATAG\n";
        HttpResponse<String> res = post("/api/analyze", mapper.writeValueAsString(
                java.util.Map.of("fasta", fasta, "organismClass", "virus", "genomeType", "rna")));
        assertThat(res.statusCode()).isEqualTo(200);
        JsonNode body = mapper.readTree(res.body());

        assertThat(body.has("gvi")).isTrue();
        assertThat(body.has("comparable")).isTrue();
        assertThat(body.has("coverageSummary")).isTrue();
        assertThat(body.has("effectiveWeightSum")).isTrue();
        assertThat(body.has("components")).isTrue();
        assertThat(body.has("exclusions")).isTrue();
        assertThat(body.has("warnings")).isTrue();
        assertThat(body.has("indices")).isTrue();

        JsonNode first = body.get("components").get(0);
        for (String field : new String[]{"key", "rawValue", "normalizedValue", "effectiveWeight", "contribution"}) {
            assertThat(first.has(field)).as("component field " + field).isTrue();
        }
    }

    /**
     * The pathogen dropdown is populated from this endpoint rather than a list hardcoded in the
     * page. Same reasoning as the codon-species endpoint, and it matters more here: the generation
     * time scales the whole birth-death process behind Re, which carries the largest weight in the
     * composite.
     */
    @Test
    void servesTheBundledGenerationTimeTable() throws Exception {
        HttpResponse<String> res = get("/api/pathogens");
        assertThat(res.statusCode()).isEqualTo(200);
        JsonNode rows = mapper.readTree(res.body()).get("pathogens");
        assertThat(rows.isArray()).isTrue();
        assertThat(rows).hasSize(org.gvi.algorithms.re.GenerationTimeTable.bundled().size());

        boolean sawUsable = false, sawUnusable = false;
        for (JsonNode row : rows) {
            for (String field : new String[]{"pathogenId", "displayName", "tDays", "confidence",
                                             "reApplicable", "vectorBorne", "note"}) {
                assertThat(row.has(field)).as("pathogen row field " + field).isTrue();
            }
            // An entry with a usable value must actually carry one, and one without must say so
            // rather than quietly presenting a null the UI could read as zero.
            if ("true".equals(row.get("reApplicable").asText())) {
                assertThat(row.get("tDays").isNull()).isFalse();
                assertThat(row.get("tDays").asDouble()).isGreaterThan(0.0);
                sawUsable = true;
            } else {
                assertThat(row.get("tDays").isNull()).isTrue();
                sawUnusable = true;
            }
        }
        assertThat(sawUsable).as("the bundled table should contain at least one usable entry").isTrue();
        assertThat(sawUnusable).as("entries without a value must still be listed, marked unusable").isTrue();
    }

    /**
     * A supplied generation time must reach the pipeline and be reported back. Before this was
     * plumbed the web layer hardcoded 5 days for every request, so a browser user had no way to
     * set the parameter Re is most sensitive to.
     */
    @Test
    void aSuppliedGenerationTimeIsUsedAndItsSourceReported() throws Exception {
        String fasta = ">a 2015-01-01\nATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG\n"
                     + ">b 2018-01-01\nATGGCCATTGTAATGGCCCGCTGAAAGGGTGCCCGATAC\n"
                     + ">c 2021-01-01\nATGGCCATTGTAATGGGCCGCTGAAAGGGTGCACGATAG\n";

        JsonNode dflt = mapper.readTree(post("/api/analyze", mapper.writeValueAsString(
                java.util.Map.of("fasta", fasta, "organismClass", "virus", "genomeType", "rna"))).body());
        assertThat(dflt.at("/summary/generationTimeDays").asDouble()).isEqualTo(5.0);
        assertThat(dflt.at("/summary/generationTimeSource").asText()).isEqualTo("default");

        JsonNode supplied = mapper.readTree(post("/api/analyze", mapper.writeValueAsString(
                java.util.Map.of("fasta", fasta, "organismClass", "virus", "genomeType", "rna",
                                 "generationTimeDays", 14.0))).body());
        assertThat(supplied.at("/summary/generationTimeDays").asDouble()).isEqualTo(14.0);
        assertThat(supplied.at("/summary/generationTimeSource").asText()).isEqualTo("supplied");
    }

    /**
     * Where an estimator produces a confidence interval it must reach the browser. Re's
     * profile-likelihood interval is the reason: on realistic tree sizes it is wide, and a UI
     * showing only the point estimate would hide exactly the part that matters.
     */
    @Test
    void anIndexIntervalIsCarriedThroughToTheResponseShape() throws Exception {
        String fasta = ">a 2015-01-01\nATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG\n"
                     + ">b 2018-01-01\nATGGCCATTGTAATGGCCCGCTGAAAGGGTGCCCGATAC\n";
        HttpResponse<String> res = post("/api/analyze", mapper.writeValueAsString(
                java.util.Map.of("fasta", fasta, "organismClass", "virus", "genomeType", "rna")));
        assertThat(res.statusCode()).isEqualTo(200);
        JsonNode indices = mapper.readTree(res.body()).get("indices");

        // Every index view carries the field, whether or not this run produced an interval --
        // the browser reads it by name and must never see it simply absent.
        indices.fields().forEachRemaining(e -> {
            assertThat(e.getValue().has("interval")).as("interval field on " + e.getKey()).isTrue();
            JsonNode iv = e.getValue().get("interval");
            if (!iv.isNull()) {
                assertThat(iv.get("lower").asDouble()).isLessThanOrEqualTo(iv.get("upper").asDouble());
                assertThat(iv.get("level").asDouble()).isBetween(0.0, 1.0);
            }
        });
    }

    /**
     * A pathogen id must resolve the generation time through the bundled table, and the response
     * must report the value the run actually used.
     * <p>
     * It did not. The resolution lived in {@code GviCli}, which runs before a PipelineConfig is
     * built, so the web layer -- which constructs the config directly -- never consulted the table.
     * A request naming blue_tongue ran on the 5-day default while the response announced
     * "bundled table (blue_tongue)": a provenance claim for a lookup that never happened. The
     * resolution now lives on GenerationTimeTable, shared by both front ends.
     */
    @Test
    void aPathogenIdResolvesTheGenerationTimeAndTheResponseSaysWhatWasUsed() throws Exception {
        String fasta = ">a 2015-01-01\nATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG\n"
                     + ">b 2018-01-01\nATGGCCATTGTAATGGCCCGCTGAAAGGGTGCCCGATAC\n"
                     + ">c 2021-01-01\nATGGCCATTGTAATGGGCCGCTGAAAGGGTGCACGATAG\n";

        var table = org.gvi.algorithms.re.GenerationTimeTable.bundled();
        String usableId = table.usableIds().stream().findFirst().orElseThrow();
        double expected = table.requireGenerationTimeDays(usableId);

        JsonNode body = mapper.readTree(post("/api/analyze", mapper.writeValueAsString(
                java.util.Map.of("fasta", fasta, "organismClass", "virus", "genomeType", "rna",
                                 "pathogenId", usableId))).body());

        assertThat(body.at("/summary/generationTimeDays").asDouble())
                .as("a pathogen id must resolve to the table's value, not the default")
                .isEqualTo(expected);
        assertThat(body.at("/summary/generationTimeSource").asText()).contains(usableId);
        // The claim and the computation must agree: reporting the table as the source while
        // running on the default is worse than either alone.
        assertThat(body.at("/summary/generationTimeDays").asDouble())
                .isNotEqualTo(AnalysisService.DEFAULT_GENERATION_TIME_DAYS);
    }

    /** An explicit value still wins over a pathogen id, matching the command line's precedence. */
    @Test
    void anExplicitGenerationTimeOverridesThePathogenId() throws Exception {
        String fasta = ">a 2015-01-01\nATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG\n"
                     + ">b 2018-01-01\nATGGCCATTGTAATGGCCCGCTGAAAGGGTGCCCGATAC\n";
        var table = org.gvi.algorithms.re.GenerationTimeTable.bundled();
        String usableId = table.usableIds().stream().findFirst().orElseThrow();

        JsonNode body = mapper.readTree(post("/api/analyze", mapper.writeValueAsString(
                java.util.Map.of("fasta", fasta, "organismClass", "virus", "genomeType", "rna",
                                 "pathogenId", usableId, "generationTimeDays", 99.0))).body());
        assertThat(body.at("/summary/generationTimeDays").asDouble()).isEqualTo(99.0);
        assertThat(body.at("/summary/generationTimeSource").asText()).isEqualTo("supplied");
    }
}
