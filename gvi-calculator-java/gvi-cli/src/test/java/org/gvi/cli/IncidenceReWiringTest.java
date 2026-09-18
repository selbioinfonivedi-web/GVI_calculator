package org.gvi.cli;

import org.gvi.algorithms.gd.GdMethod;
import org.gvi.algorithms.mu.GenomeType;
import org.gvi.algorithms.re.ReMethod;
import org.gvi.algorithms.re.SerialInterval;
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
     * <p>
     * The infectiousness profile comes from {@link SerialInterval} itself rather than a
     * reimplementation here. That is deliberate, and was learned the hard way: an earlier version of
     * this test discretised the gamma by CDF differencing, {@code F(u) - F(u-1)}, which looks
     * equivalent and is not. It assigns each interval's mass to its right endpoint and shifts the
     * effective mean from 5.00 to 5.50 days, so the estimator was being asked to invert data from a
     * half-day-slower epidemic than the one it had been told about. That cost 0.008 at Re = 1.15 and
     * 0.127 at Re = 2.0 -- invisible at the one point originally tested, disqualifying across the
     * range. Sharing the profile keeps the oracle honest rather than weakening it: the forward
     * simulation is still independent of the inversion under test, and the serial interval is a
     * parameter both sides are entitled to agree on.
     */
    private Path incidenceCsv(double rTrue) throws IOException {
        SerialInterval si = new SerialInterval(5.0, 2.0, 20);
        int tMax = si.tMax();

        int totalDays = 70;
        double[] incidence = new double[totalDays];
        for (int i = 0; i < tMax; i++) incidence[i] = 150.0;
        Random rng = new Random(SEED);
        for (int t = tMax; t < totalDays; t++) {
            double lambda = 0;
            for (int u = 1; u <= tMax; u++) lambda += incidence[t - u] * si.weight(u);
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
                GdMethod.JUKES_CANTOR, 1.0, 2.0, false, null, "gtr", false, false, false, true, false,
                null, GenomeType.RNA, Map.of(), null, false, 100, false, null);
    }

    // ── the wiring ───────────────────────────────────────────────────────────
    /**
     * Recovery is checked across the range rather than at one point, and the tolerance is calibrated
     * rather than guessed.
     * <p>
     * A single point with a generous band is a weak oracle twice over: an estimator that returned a
     * constant near that point would pass, and a band wide enough to absorb Poisson noise also
     * absorbs a subtly wrong estimator. Measured on this fixture, the error at a fixed seed is
     * 0.001-0.04 and the spread across twelve seeds is 0.0496 wide, so 0.06 clears the noise with
     * headroom. For scale, a serial interval wrong by a factor of two moves the estimate by 0.136 --
     * comfortably outside this band, where the 0.10 originally used here would have swallowed it.
     */
    @Test
    void recoversTheTrueReAcrossTheRangeNotJustAtOnePoint() throws IOException {
        double[] rTrue = {0.80, 1.15, 1.50, 2.00};
        double[] recovered = new double[rTrue.length];

        for (int i = 0; i < rTrue.length; i++) {
            PipelineResult result = new GviPipeline().run(config(fasta(), incidenceCsv(rTrue[i])));
            var re = result.populationIndices().get(IndexKey.RE);
            assertThat(re).as("Re should be scored at rTrue=%.2f; skipped: %s", rTrue[i], result.skipped())
                    .isNotNull();
            assertThat(((org.gvi.algorithms.re.ReResult) re).method()).isEqualTo(ReMethod.CORI_INCIDENCE);

            recovered[i] = re.primaryValue();
            assertThat(recovered[i]).as("recovering rTrue=%.2f", rTrue[i]).isCloseTo(rTrue[i], within(0.06));
        }

        // A constant-returning stub passes any single-point check; it cannot pass this one.
        assertThat(recovered).as("the estimate must track the truth, not sit at a fixed value")
                .isSorted();

        // The epidemiologically load-bearing property: growing and shrinking must not be confused.
        assertThat(recovered[0]).as("rTrue=0.80 must read as a shrinking outbreak").isLessThan(1.0);
        assertThat(recovered[2]).as("rTrue=1.50 must read as a growing one").isGreaterThan(1.0);
    }

    /**
     * That the case counts <em>decided</em> the answer, not merely labelled it.
     * <p>
     * The weak form of this test -- "with a file, method is CORI" -- passes even if Cori won by
     * walkover because the birth-death fit failed on the fixture. It does not: BDSKY succeeds here
     * and returns about 0.996, while Cori on the same alignment returns about 1.148. So the two
     * routes disagree materially, and the incidence run must land on Cori's answer rather than
     * BDSKY's. A refactor that read the file and then ignored it would still fail this.
     */
    @Test
    void theCaseCountsDecideTheAnswerAndNotJustTheLabel() throws IOException {
        PipelineResult without = new GviPipeline().run(config(fasta(), null));
        var treeRe = (org.gvi.algorithms.re.ReResult) without.populationIndices().get(IndexKey.RE);

        assertThat(treeRe).as("the fixture must be one the tree-shape route can actually fit, "
                + "otherwise Cori wins by default and this test proves nothing").isNotNull();
        assertThat(treeRe.method())
                .as("without case counts Re must come from tree shape")
                .isNotEqualTo(ReMethod.CORI_INCIDENCE);

        PipelineResult withCases = new GviPipeline().run(config(fasta(), incidenceCsv(1.15)));
        var coriRe = (org.gvi.algorithms.re.ReResult) withCases.populationIndices().get(IndexKey.RE);

        assertThat(coriRe.method()).isEqualTo(ReMethod.CORI_INCIDENCE);
        assertThat(coriRe.primaryValue())
                .as("the two routes must disagree here, or the comparison below is vacuous")
                .isNotCloseTo(treeRe.primaryValue(), within(0.05));
        assertThat(coriRe.primaryValue())
                .as("the answer must be the one the case counts imply, not the tree's")
                .isCloseTo(1.15, within(0.06));
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
