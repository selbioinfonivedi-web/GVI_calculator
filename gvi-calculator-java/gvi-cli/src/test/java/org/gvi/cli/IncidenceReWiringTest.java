package org.gvi.cli;

import org.gvi.algorithms.gd.GdMethod;
import org.gvi.algorithms.mu.GenomeType;
import org.gvi.algorithms.re.ReMethod;
import org.gvi.composite.IndexKey;
import org.gvi.core.model.OrganismClass;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.time.LocalDate;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.within;

/**
 * Case incidence is the accurate path for Re: the Cori et al. (2013) estimator recovers a known
 * value to within about 0.02, where tree-shape inference does far worse. The estimator itself is
 * well covered by {@code ReGroundTruthRecoveryTest}, but nothing exercised the <em>wiring</em> --
 * that a real CSV on disk reaches it, that its presence actually switches the pipeline away from
 * the birth-death fit, and that the resulting Re lands in the composite rather than being gated
 * out. No corpus dataset supplies incidence data, so this path shipped unexercised end to end.
 * <p>
 * Re carries the largest single weight in the scheme, so "the good estimator exists but nothing
 * routes to it" is an expensive failure to leave untested.
 */
class IncidenceReWiringTest {

    @TempDir
    Path tempDir;

    private static final long SEED = 20260813L;

    /** Same alignment shape the other pipeline tests use: enough signal to build a tree from. */
    private static final String FASTA = """
            >reference 2015-01-01
            ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAGATGGCCATTGTAATGGGCCGC
            >q1 2018-01-01
            ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAGATGGCCATTGTAATGGCCCGC
            >q2 2021-01-01
            ATGGCCATTGTAATGGCCCGCTGAAAGGGTGCCCGATACATGGCCATTGTAATGGGCCGC
            >q3 2023-01-01
            ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCACGATAGATGGCCATTGTAATGGGCCGA
            """;

    private Path fasta() throws IOException {
        Path p = tempDir.resolve("in.fasta");
        Files.writeString(p, FASTA);
        return p;
    }

    /**
     * A renewal process with a known Re and Poisson case counts -- the generative model the Cori
     * estimator is derived to invert, so recovering rTrue from the file proves the whole chain.
     */
    private Path incidenceCsv(double rTrue) throws IOException {
        int tMax = 20;
        double mean = 5.0, sd = 2.0;
        double[] w = serialIntervalWeights(mean, sd, tMax);

        int totalDays = 70;
        double[] incidence = new double[totalDays];
        for (int i = 0; i < tMax; i++) incidence[i] = 150.0;
        Random rng = new Random(SEED);
        for (int t = tMax; t < totalDays; t++) {
            double lambda = 0;
            for (int u = 1; u <= tMax; u++) lambda += incidence[t - u] * w[u];
            incidence[t] = poissonSample(rTrue * lambda, rng);
        }

        StringBuilder sb = new StringBuilder("date,new_cases\n");
        LocalDate start = LocalDate.of(2024, 1, 1);
        for (int i = 0; i < totalDays; i++) {
            sb.append(start.plusDays(i)).append(',').append((long) incidence[i]).append('\n');
        }
        Path p = tempDir.resolve("incidence.csv");
        Files.writeString(p, sb);
        return p;
    }

    /** Discretised gamma, matching SerialInterval's own convention. */
    private static double[] serialIntervalWeights(double mean, double sd, int tMax) {
        double shape = (mean / sd) * (mean / sd);
        double scale = sd * sd / mean;
        var gamma = new org.apache.commons.math3.distribution.GammaDistribution(shape, scale);
        double[] w = new double[tMax + 1];
        double sum = 0;
        for (int u = 1; u <= tMax; u++) {
            w[u] = gamma.cumulativeProbability(u) - gamma.cumulativeProbability(u - 1);
            sum += w[u];
        }
        for (int u = 1; u <= tMax; u++) w[u] /= sum;
        return w;
    }

    private static double poissonSample(double lambda, Random rng) {
        if (lambda <= 0) return 0;
        if (lambda > 500) return Math.max(0, Math.round(lambda + Math.sqrt(lambda) * rng.nextGaussian()));
        double l = Math.exp(-lambda), p = 1.0;
        int k = 0;
        do { k++; p *= rng.nextDouble(); } while (p > l);
        return k - 1;
    }

    private PipelineConfig config(Path fasta, Path incidence) {
        return new PipelineConfig(
                fasta, null, null, null, null, incidence, null, "reference",
                Set.of("re", "pi", "gd"),
                OrganismClass.VIRUS, null, true, 5.0,
                GdMethod.JUKES_CANTOR, 1.0, 2.0, false, null, "gtr", false, false, true, false,
                null, GenomeType.RNA, Map.of(), null, false, 100, false, null);
    }

    // ── the wiring ───────────────────────────────────────────────────────────
    @Test
    void anIncidenceFileOnDiskReachesTheCoriEstimatorAndRecoversTheTrueRe() throws IOException {
        double rTrue = 1.15;
        PipelineResult result = new GviPipeline().run(config(fasta(), incidenceCsv(rTrue)));

        var re = result.populationIndices().get(IndexKey.RE);
        assertThat(re).as("Re should be scored when incidence is supplied; skipped: %s", result.skipped())
                .isNotNull();

        assertThat(((org.gvi.algorithms.re.ReResult) re).method())
                .as("supplying real case counts must select Cori, not the birth-death or LTT fallback")
                .isEqualTo(ReMethod.CORI_INCIDENCE);

        assertThat(re.primaryValue())
                .as("the whole chain -- CSV on disk, reader, estimator -- should recover the simulated Re")
                .isCloseTo(rTrue, within(0.10));
    }

    /**
     * The contrast that makes the previous test mean something: the same alignment with no incidence
     * file takes a different estimator entirely. If both routes produced the same method, the first
     * assertion would pass without the file being read at all.
     */
    @Test
    void thesameAlignmentWithoutIncidenceUsesADifferentEstimator() throws IOException {
        PipelineResult withCases = new GviPipeline().run(config(fasta(), incidenceCsv(1.15)));
        PipelineResult without = new GviPipeline().run(config(fasta(), null));

        var a = (org.gvi.algorithms.re.ReResult) withCases.populationIndices().get(IndexKey.RE);
        var b = (org.gvi.algorithms.re.ReResult) without.populationIndices().get(IndexKey.RE);

        assertThat(a.method()).isEqualTo(ReMethod.CORI_INCIDENCE);
        if (b != null) {
            assertThat(b.method())
                    .as("without case counts Re must come from tree shape, not Cori")
                    .isNotEqualTo(ReMethod.CORI_INCIDENCE);
        }
    }

    /** Re carries the largest weight in the scheme, so it has to actually enter the composite. */
    @Test
    void anIncidenceDerivedReIsNotGatedOutOfTheComposite() throws IOException {
        PipelineResult result = new GviPipeline().run(config(fasta(), incidenceCsv(1.15)));

        assertThat(result.datasetGvi().excludedIndices())
                .as("Re came from real case counts; it must not be excluded")
                .doesNotContain(IndexKey.RE);
        assertThat(result.datasetGvi().components())
                .anyMatch(c -> c.key() == IndexKey.RE && c.effectiveWeight() > 0);
    }

    /**
     * The default-generation-time warning exists for tree-derived Re. Cori works from the serial
     * interval in the case counts themselves, so firing it here would tell the user to fix
     * something that is not affecting their number.
     */
    @Test
    void theGenerationTimeWarningIsNotRaisedWhenReCameFromCaseCounts() throws IOException {
        PipelineResult result = new GviPipeline().run(config(fasta(), incidenceCsv(1.15)));

        assertThat(result.warnings())
                .as("warnings: %s", result.warnings())
                .noneMatch(w -> w.contains("left at its default (5 days"));
    }

    /** Parse warnings from the incidence file must survive into the run's own warnings. */
    @Test
    void skippedIncidenceRowsAreReportedInTheRunsWarnings() throws IOException {
        Path csv = tempDir.resolve("dirty.csv");
        StringBuilder sb = new StringBuilder("date,new_cases\n");
        LocalDate start = LocalDate.of(2024, 1, 1);
        for (int i = 0; i < 40; i++) {
            sb.append(start.plusDays(i)).append(',').append(100 + i).append('\n');
        }
        sb.append("2024-03-01,-5\n");           // negative: skipped
        sb.append("not-a-date,7\n");            // unparsable: skipped
        Files.writeString(csv, sb);

        PipelineResult result = new GviPipeline().run(config(fasta(), csv));

        assertThat(result.warnings()).as("warnings: %s", result.warnings())
                .anyMatch(w -> w.toLowerCase().contains("incidence"));
    }

    /**
     * A named-but-absent file is the user's error. The pipeline skips Re outright rather than
     * quietly falling back to tree shape -- which matters, because a silent fallback would hand back
     * a materially worse Re under the impression it came from the case counts that were asked for.
     */
    @Test
    void anAbsentIncidenceFileSkipsReRatherThanSilentlyFallingBackToTreeShape() throws IOException {
        PipelineResult result = new GviPipeline().run(config(fasta(), tempDir.resolve("nope.csv")));

        assertThat(result.skipped())
                .as("skipped=%s", result.skipped())
                .anyMatch(s -> s.startsWith("Re:") && s.contains("Cannot read incidence CSV")
                        && s.contains("nope.csv"));
        assertThat(result.populationIndices().get(IndexKey.RE))
                .as("Re must not be scored from another estimator behind the user's back")
                .isNull();
    }
}
