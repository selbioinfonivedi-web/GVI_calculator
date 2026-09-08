package org.gvi.cli;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import picocli.CommandLine;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;

import static org.assertj.core.api.Assertions.assertThat;

/**
 * The exit codes are a documented contract -- 0 success, 1 bad input/usage, 2 computation could not
 * complete, 3 self-test failed, 9 internal bug -- and scripts branch on them. Nothing verified them,
 * and {@link GviCli} was at 0% coverage, so a refactor that turned a bad path into an exit 9 (a bug
 * report) instead of an exit 1 (a user error) would have gone unnoticed.
 * <p>
 * These drive the real command line rather than calling methods, because the class under test is the
 * argument parsing and error mapping, not the pipeline.
 */
class GviCliTest {

    @TempDir
    Path tempDir;

    private String out = "";
    private String err = "";

    /** Runs the CLI as a script would, capturing both streams. Returns the exit code. */
    private int run(String... args) {
        ByteArrayOutputStream o = new ByteArrayOutputStream();
        ByteArrayOutputStream e = new ByteArrayOutputStream();
        PrintStream oldOut = System.out;
        PrintStream oldErr = System.err;
        try {
            System.setOut(new PrintStream(o, true, StandardCharsets.UTF_8));
            System.setErr(new PrintStream(e, true, StandardCharsets.UTF_8));
            CommandLine cmd = new CommandLine(new GviCli());
            cmd.setOut(new PrintWriter(o, true));
            cmd.setErr(new PrintWriter(e, true));
            return cmd.execute(args);
        } finally {
            System.setOut(oldOut);
            System.setErr(oldErr);
            out = o.toString(StandardCharsets.UTF_8);
            err = e.toString(StandardCharsets.UTF_8);
        }
    }

    private Path fastaWith(String body) throws IOException {
        Path p = tempDir.resolve("in.fasta");
        Files.writeString(p, body);
        return p;
    }

    private static final String VALID_FASTA = """
            >reference 2015-01-01
            ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAGATGGCCATTGTAATGGGCCGC
            >q1 2018-01-01
            ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAGATGGCCATTGTAATGGCCCGC
            >q2 2021-01-01
            ATGGCCATTGTAATGGCCCGCTGAAAGGGTGCCCGATACATGGCCATTGTAATGGGCCGC
            """;

    // ── the success path ─────────────────────────────────────────────────────
    @Test
    void aRunThatProducesAScoreExitsZeroAndPrintsIt() throws IOException {
        int code = run("--fasta", fastaWith(VALID_FASTA).toString(),
                       "--indices", "pi,gd,mb,ri", "--organism-class", "virus");

        assertThat(code).as("stderr was: %s", err).isZero();
        assertThat(out).contains("GVI = ");
    }

    @Test
    void theJsonOutputPathIsWrittenAndParseable() throws IOException {
        Path json = tempDir.resolve("out.json");
        int code = run("--fasta", fastaWith(VALID_FASTA).toString(),
                       "--indices", "pi,gd", "--organism-class", "virus",
                       "--json", json.toString());

        assertThat(code).isZero();
        assertThat(json).exists();
        assertThat(Files.readString(json)).contains("\"dataset_gvi\"");
    }

    /** --self-test is what an installer runs to prove the binary works; it must pass here. */
    @Test
    void theSelfTestPassesInsideThisBuildAndExitsZero() {
        int code = run("--self-test");

        assertThat(code).as("a failing self-test exits 3; output was: %s", out).isZero();
        assertThat(out).contains("Self-Test");
    }

    // ── exit 1: the user's input or usage is wrong ───────────────────────────
    @Test
    void noFastaAndNoSelfTestIsAUsageErrorNotACrash() {
        int code = run("--organism-class", "virus");

        assertThat(code).isEqualTo(1);
        assertThat(err).contains("--fasta is required");
    }

    @Test
    void aMissingFastaFileIsAnInputErrorNotAnInternalOne() {
        int code = run("--fasta", tempDir.resolve("absent.fasta").toString());

        assertThat(code).as("a file the user named but did not provide is their error, not a bug")
                .isEqualTo(1);
        assertThat(err).startsWith("Error: ");
        assertThat(err).doesNotContain("Unexpected internal error");
    }

    @Test
    void aFileThatIsNotFastaIsAnInputError() throws IOException {
        int code = run("--fasta", fastaWith("this is not a FASTA file at all\n").toString());

        assertThat(code).isEqualTo(1);
        assertThat(err).doesNotContain("Unexpected internal error");
    }

    /** An unparseable option is picocli's own usage failure, which must not read as a crash. */
    @Test
    void anUnknownOptionIsAUsageErrorWithUsageText() {
        int code = run("--not-an-option");

        assertThat(code).isNotZero();
        assertThat(code).as("a typo must not report itself as an internal bug").isNotEqualTo(9);
        assertThat(err).containsIgnoringCase("unknown option");
    }

    /** No stack trace ever reaches the user: that is the whole point of the catch in call(). */
    @Test
    void noErrorPathPrintsARawStackTrace() {
        run("--fasta", tempDir.resolve("absent.fasta").toString());

        assertThat(err).doesNotContain("\tat org.gvi");
        assertThat(err).doesNotContain("Exception in thread");
    }

    // ── calibration, which has its own entrypoint before --fasta is checked ──
    @Test
    void calibrationRunsWithoutAFastaBecauseItFitsWeightsNotSequences() throws IOException {
        Path csv = tempDir.resolve("calib.csv");
        Files.writeString(csv, """
                mu,pi,GD,target
                0.10,0.20,0.30,0.55
                0.20,0.30,0.35,0.70
                0.35,0.45,0.50,0.95
                0.50,0.55,0.60,1.20
                0.65,0.70,0.75,1.45
                """);

        int code = run("--calibrate", csv.toString());

        assertThat(code).as("stderr was: %s", err).isZero();
        assertThat(out).contains("Calibration");
    }

    @Test
    void aBrokenCalibrationFileIsAnInputErrorNamingTheProblem() throws IOException {
        Path csv = tempDir.resolve("bad.csv");
        Files.writeString(csv, "mu,pi\n0.1,0.2\n");

        int code = run("--calibrate", csv.toString());

        assertThat(code).isEqualTo(1);
        assertThat(err).contains("target");
    }

    // ── --help / --version, which scripts and packagers call ────────────────
    @Test
    void helpAndVersionExitZero() {
        assertThat(run("--help")).isZero();
        assertThat(out).contains("--fasta");
        assertThat(run("--version")).isZero();
        assertThat(out).contains("GVI Calculator");
    }
}
