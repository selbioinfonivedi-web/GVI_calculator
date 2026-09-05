package org.gvi.cli;

import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.SerializationFeature;
import org.gvi.algorithms.gd.GdMethod;
import org.gvi.composite.IndexKey;
import org.gvi.composite.WeightCalibrator;
import org.gvi.core.exception.GviException;
import org.gvi.core.util.GlobalExceptionHandler;
import org.gvi.selftest.SelfTestCase;
import org.gvi.selftest.SelfTestReport;
import org.gvi.selftest.SelfTestRunner;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import picocli.CommandLine;
import picocli.CommandLine.Command;
import picocli.CommandLine.Option;

import java.io.IOException;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Callable;

/**
 * Headless CLI entrypoint (Section 6/7 of the build spec). Every error path
 * is caught here and turned into a clear message and a non-zero exit code
 * -- nothing propagates as a raw stack trace to the user, and nothing calls
 * System.exit() from deeper in the call stack.
 * <p>
 * Exit codes: 0 = success, 1 = bad input/usage, 2 = computation could not
 * complete, 3 = self-test failed, 9 = unexpected internal error (a bug).
 */
@Command(name = "gvi-calculator", mixinStandardHelpOptions = true, version = "GVI Calculator 0.1.0",
        description = "Standalone calculator for the 8 Genomic Virulence Index component indices and the composite GVI.")
public final class GviCli implements Callable<Integer> {

    private static final Logger log = LoggerFactory.getLogger(GviCli.class);

    @Option(names = "--self-test", description = "Run the bundled self-diagnostic suite and exit (validates this installation offline).")
    boolean selfTest;

    @Option(names = "--fasta", description = "Multi-FASTA alignment (required unless --self-test).")
    Path fasta;

    @Option(names = "--metadata", description = "Metadata CSV (sequence_id, collection_date, location, host).")
    Path metadata;

    @Option(names = "--gff", description = "GFF3 gene annotation (CDS features) for dN/dS and CAI.")
    Path gff;

    @Option(names = "--gff-out", description = "Write the native ORF-scan's predicted gene coordinates to this path as a real "
            + "GFF3 file -- only applies when no --gff was supplied (or it failed to parse) and the pipeline had to guess gene "
            + "boundaries itself for dN/dS/CAI. Lets you inspect what was guessed, hand-correct it, or feed it back in via "
            + "--gff on a later run instead of re-guessing every time. No-op if --gff was supplied and used successfully.")
    Path gffOut;

    @Option(names = "--codon-usage", description = "Host codon usage table CSV (codon,frequency) for CAI.")
    Path codonUsage;

    @Option(names = "--codon-usage-species", description = "Use a bundled, real, precomputed host codon usage table "
            + "instead of --codon-usage -- fetched directly from the Kazusa Codon Usage Database (Nakamura, Gojobori & "
            + "Ikemura 2000), not synthesized. Available: human, mouse, pig, wild_boar, cattle, buffalo, sheep, goat, "
            + "horse, ecoli, aedes_aegypti (see BundledCodonUsageTables "
            + "for full provenance -- CDS/codon counts, fetch date). Ignored if --codon-usage is also supplied (an "
            + "explicit table always wins). A bundled table is a real default, not a substitute for a dataset-specific "
            + "reference built from the actual host's own highly-expressed genes -- prefer --codon-usage whenever you "
            + "have real host data.")
    String codonUsageSpecies;

    @Option(names = "--incidence", description = "Case-incidence CSV (date,new_cases) for the Re estimator.")
    Path incidence;

    @Option(names = "--reference-gc", description = "Host/reference genome GC%% for GC Content Deviation.")
    Double referenceGc;

    @Option(names = "--reference-id", description = "Sequence id to use as reference (default: first sequence in the FASTA).")
    String referenceId;

    @Option(names = "--indices", split = ",", description = "Which indices to compute: mu,re,pi,mb,dnds,gd,cai,gc,ri,all (default: all).")
    String[] indices = {"all"};

    @Option(names = "--gd-method", description = "hamming|jukes_cantor|kimura_2_parameter (default: jukes_cantor)")
    String gdMethod = "jukes_cantor";

    @Option(names = "--genome-type", description = "rna|dna -- this pathogen's genome chemistry, used to classify the computed mu "
            + "against the right reference table. Without it, mu is classified against RNA virus reference points (Influenza, "
            + "SARS-CoV-2, VSV) regardless of what's actually being analyzed -- misleading for a DNA virus (e.g. a poxvirus) or a "
            + "non-coding marker, where the resulting rate is many orders of magnitude below (or otherwise incomparable to) any "
            + "real RNA virus. Omit only when you genuinely don't know or the marker isn't from a virus genome at all.")
    String genomeType;

    @Option(names = "--organism-class", description = "virus|bacterium|parasite -- what kind of organism this is. Gates the "
            + "indices whose biological premise depends on it. CAI against a HOST codon table only means anything for an "
            + "organism that translates on host ribosomes (a virus); bacteria and eukaryotic parasites have their own "
            + "ribosomes and tRNA pools, so for those CAI is excluded from the composite rather than scored against an "
            + "unrelated organism. Also gates whether GC deviation is described as possible horizontal gene transfer "
            + "(documented in bacteria) or as compositional/mutational-pressure drift (viruses). Omit to disable gating.")
    String organismClass;

    @Option(names = "--pathogen-id", description = "Key into the bundled per-pathogen generation-time table "
            + "(generation_times.yaml), e.g. 'fmd', 'sars_cov_2'. Supplies the serial-interval mean the phylodynamic Re "
            + "fallback needs. Use this INSTEAD of --generation-time-days: a single default applied to every organism "
            + "compresses genuinely different growth rates toward Re=1. If the id has no usable entry the run fails with "
            + "an explanation rather than silently substituting a default.")
    String pathogenId;

    @Option(names = "--generation-times", description = "Override the bundled generation-time table with your own YAML "
            + "file using the same schema.")
    Path generationTimesPath;

    @Option(names = "--trim-to-covered", description = "Trim the alignment to the columns every sequence actually covers "
            + "before computing anything. Use this when the input pools PARTIAL sequences with different coverage windows "
            + "(the same locus sequenced to differing lengths): those align correctly but leave the short ones gap-padded, "
            + "and every position-wise index then compares real bases against absent data. Note this is not a re-alignment "
            + "-- re-running MAFFT will not help, because the aligner already placed the sequences correctly.")
    boolean trimToCovered;

    @Option(names = "--min-trimmed-columns", description = "Refuse to trim below this many columns (default: 100). Guards "
            + "against silently turning a poorly-overlapping set into a tiny alignment that computes cleanly but means "
            + "nothing.")
    int minTrimmedColumns = 100;

    @Option(names = "--figures", description = "Write publication-quality figures (vector SVG, single-column width, "
            + "colour-vision-safe, greyscale-legible) to this directory: composition of the score, an all-index panel, "
            + "and a weight-sensitivity tornado. Excluded indices are drawn hatched rather than omitted, because a GVI "
            + "built from six indices is not the same quantity as one built from nine.")
    Path figuresDir;

    @Option(names = "--auto", description = "Aligned-FASTA-only mode: derive everything derivable from the alignment "
            + "itself, then compute. Collection dates, location and host are read from the FASTA headers; gene "
            + "boundaries are predicted by the native ORF scan; the reference sequence is chosen automatically; and the "
            + "alignment is trimmed to its shared coverage window IF (and only if) it would otherwise fail the "
            + "coverage quality gate. Everything derived is written out beside --json so it can be inspected and reused. "
            + "Two things cannot be read from sequence and are still worth supplying: --organism-class (gates whether "
            + "host-relative CAI is meaningful) and --genome-type (sets which scale mu is judged against).")
    boolean auto;

    @Option(names = "--segment", description = "label=path to an independently-aligned genome SEGMENT, for reassortment "
            + "detection between segments of a segmented-genome pathogen (e.g. Bluetongue's 10 segments, influenza's 8). "
            + "Repeatable -- supply 2 or more (e.g. --segment seg2=vp2.fasta --segment seg10=ns3.fasta) to run a pairwise "
            + "Mantel-test reassortment check across all segment pairs; sequence ids must match exactly across segment "
            + "files for the same isolate to be compared. Entirely separate from --fasta (the single-locus alignment "
            + "every other index runs against) -- pass --fasta as one of the --segment entries too if you want it "
            + "included. Not folded into the composite GVI; reported as its own section.")
    Map<String, Path> segments;

    @Option(names = "--high-accuracy-mu", description = "Estimate mu via GTR(+Gamma) maximum-likelihood branch-length "
            + "optimization instead of plain JC69+Neighbor-Joining -- a real accuracy upgrade, but much slower "
            + "(capped at " + org.gvi.algorithms.mu.EvolutionaryRateCalculator.MAX_TAXA_FOR_ML + " taxa / "
            + org.gvi.algorithms.mu.EvolutionaryRateCalculator.MAX_SITES_FOR_ML + " sites; falls back automatically beyond that).")
    boolean highAccuracyMu;

    @Option(names = "--gamma-alpha", description = "Starting value for Gamma among-site rate heterogeneity's alpha shape "
            + "parameter, for --high-accuracy-mu -- genuinely ML-refined from this seed (not fixed); omit to disable Gamma entirely.")
    Double gammaAlpha;

    @Option(names = "--substitution-model", description = "Substitution model for --high-accuracy-mu: jc69, f81, k80, "
            + "hky85, tn93, gtr (default), or auto to fit all 6 and pick the best by AIC (slower -- up to 6x the cost).")
    String substitutionModel = "gtr";

    @Option(names = "--bootstrap-support", description = "Run a Felsenstein (1985) nonparametric bootstrap ("
            + org.gvi.algorithms.phylo.BootstrapSupportCalculator.DEFAULT_REPLICATES + " replicates) on mu's Neighbor-Joining "
            + "tree and report per-clade support -- the same 'bootstrap support %%' RAxML/IQ-TREE report. Slower (multiplies "
            + "tree-building cost by the replicate count); capped at "
            + org.gvi.algorithms.phylo.BootstrapSupportCalculator.MAX_TAXA_FOR_BOOTSTRAP + " taxa.")
    boolean bootstrapSupport;

    @Option(names = "--lsd-mu", description = "Estimate mu via native least-squares divergence-time dating -- the same objective "
            + "LSD2 (To et al. 2016) solves behind IQ-TREE's --date option: jointly fits the rate AND every internal node's date "
            + "against every edge directly, instead of one root-to-tip regression line. Needs every sequence dated; falls back "
            + "automatically to the standard estimator otherwise.")
    boolean lsdMu;

    @Option(names = "--bdsky-re", description = "Estimate Re via a native birth-death-sampling maximum-likelihood fit -- the same "
            + "generative model BEAST2's BDSKY package uses, fit by optimization instead of full Bayesian MCMC. Only tried when "
            + "no --incidence data is supplied (Cori et al. remains preferred whenever real case counts exist). Needs every "
            + "sequence dated; capped at " + org.gvi.algorithms.re.bdsky.BdskyReEstimator.MAX_TAXA_FOR_BDSKY
            + " taxa; falls back automatically to the phylodynamic fallback otherwise.")
    boolean bdskyRe;

    @Option(names = "--ml-dnds", description = "Also compute dN/dS via a native maximum-likelihood codon-substitution model "
            + "(GY94 -- the same core model PAML's codeml builds its M0 'one-ratio' analysis on): kappa, omega, and every branch "
            + "length jointly ML-fit via Felsenstein pruning over the 61 sense codons, instead of Nei-Gojobori counting. This is "
            + "ADDITIONAL cross-checking detail per gene -- it does not replace or feed into the composite GVI's dN/dS component, "
            + "which remains the already-validated Nei-Gojobori method. Capped at " + org.gvi.algorithms.dnds.ml.MlCodonDnDsEstimator.MAX_TAXA_FOR_ML_DNDS
            + " taxa / " + org.gvi.algorithms.dnds.ml.MlCodonDnDsEstimator.MAX_CODONS_FOR_ML_DNDS
            + " codons per gene (61-state codon likelihood is expensive); genes exceeding the cap are skipped with a clear reason.")
    boolean mlDnds;

    @Option(names = "--generation-time-days", description = "Generation time in days, used by the phylodynamic Re fallback (default: 5.0).")
    double generationTimeDays = 5.0;

    @Option(names = "--beta0", description = "Baseline transmission rate for the beta_genomic output (default: 1.0).")
    double beta0 = 1.0;

    @Option(names = "--scale-factor", description = "GVI-to-beta scale factor (default: 2.0).")
    double scaleFactor = 2.0;

    @Option(names = "--json", description = "Write the full structured report to this JSON file.")
    Path jsonOut;

    @Option(names = "--csv", description = "Write a summary CSV to this file -- one row (the whole-file result) by "
            + "default, or one row per sequence if --per-sequence is also set.")
    Path csvOut;

    @Option(names = "--per-sequence", description = "Also include per-sequence detail in the report -- by default "
            + "only the whole-file result is shown/written (text report, --json, --csv), since the whole-file result "
            + "is the headline answer and per-sequence detail is supporting evidence, not a second result. Pass this "
            + "to see/export the per-sequence breakdown too.")
    boolean perSequence;

    @Option(names = "--calibrate", description = "Fit composite GVI weights against historical data instead of running the pipeline "
            + "(CSV: index-label columns + 'target'; Section 8.5's optional calibration module). Exits after printing the fitted weights.")
    Path calibrateAgainst;

    @Option(names = "--holdout-fraction", description = "Fraction of --calibrate data held out for validation (default: 0.2).")
    double holdoutFraction = 0.2;

    @Option(names = "--weights-out", description = "With --calibrate, write the fitted weights to this JSON file for reuse with --weights.")
    Path weightsOut;

    @Option(names = "--weights", description = "Use custom composite GVI weights instead of the spec's default midpoints, for a real "
            + "run (not --calibrate). JSON object mapping index label to weight, e.g. {\"Re\": 0.30, \"MB\": 0.22}. Accepts exactly "
            + "the format --weights-out writes, so --calibrate --weights-out fitted.json then --weights fitted.json round-trips a "
            + "fitted set into real use. The 9 valid labels: mu, Re, pi, MB, dN/dS, GD, CAI, GC_Deviation, RI. Any label you don't "
            + "mention keeps its spec-default weight; the full set (yours plus whatever defaults filled in) is renormalized to sum "
            + "to 1.0, so you only need to specify what you actually want to change.")
    Path weightsPath;

    public static void main(String[] args) {
        GlobalExceptionHandler.install((thread, throwable) ->
                System.err.println("FATAL: an unexpected internal error occurred. See the log file for details. (" + throwable + ")"));
        int exitCode = new CommandLine(new GviCli()).execute(args);
        System.exit(exitCode);
    }

    @Override
    public Integer call() {
        try {
            if (selfTest) {
                return runSelfTest();
            }
            if (calibrateAgainst != null) {
                return runCalibration();
            }
            if (fasta == null) {
                System.err.println("Error: --fasta is required (or pass --self-test). Use --help for usage.");
                return 1;
            }
            return runPipeline();
        } catch (GviException e) {
            System.err.println("Error: " + e.getMessage());
            log.warn("Handled application error", e);
            return e.getClass().getSimpleName().equals("GviInputException") ? 1 : 2;
        } catch (Exception e) {
            System.err.println("Unexpected internal error: " + e.getMessage() + " (see log file for the full stack trace; this is likely a bug)");
            log.error("Unhandled exception in CLI", e);
            return 9;
        }
    }

    private int runSelfTest() {
        SelfTestReport report = new SelfTestRunner().runAll();
        System.out.println("=== GVI Calculator Self-Test ===");
        System.out.println("Java: " + System.getProperty("java.version") + "  OS: " + System.getProperty("os.name") + " " + System.getProperty("os.version"));
        System.out.println();
        for (SelfTestCase c : report.cases()) {
            System.out.printf("[%s] %s%n", c.passed() ? "PASS" : "FAIL", c.name());
            if (!c.passed()) {
                System.out.println("      " + c.detail());
            }
        }
        System.out.println();
        System.out.println(report.allPassed()
                ? "All self-tests passed. This installation is computing all 8 indices correctly."
                : "SELF-TEST FAILURES DETECTED. Do not trust results from this installation until this is resolved.");
        return report.allPassed() ? 0 : 3;
    }

    private int runCalibration() {
        CalibrationCsvReader.Result data = CalibrationCsvReader.read(calibrateAgainst);
        System.out.println("=== GVI Weight Calibration ===");
        System.out.println(data.observations().size() + " observation(s), fitting weights for: "
                + data.keysInPlay().stream().map(IndexKey::label).reduce((a, b) -> a + ", " + b).orElse(""));

        WeightCalibrator.CalibrationResult result;
        if (data.observations().size() >= 5) {
            result = new WeightCalibrator().calibrateWithHoldout(data.observations(), data.keysInPlay(), holdoutFraction, 42L);
        } else {
            System.out.println("(Fewer than 5 observations -- skipping hold-out split, fitting on all data.)");
            result = new WeightCalibrator().calibrate(data.observations(), data.keysInPlay());
        }

        System.out.println();
        System.out.println("Fitted weights:");
        Map<String, Double> weightsMap = new LinkedHashMap<>();
        for (IndexKey key : data.keysInPlay()) {
            double w = result.weights().get(key);
            weightsMap.put(key.label(), w);
            System.out.printf("  %-14s %.4f%n", key.label(), w);
        }
        System.out.printf("%nTrain RMSE: %.4f%n", result.trainRmse());
        if (result.holdoutRmse() != null) {
            System.out.printf("Hold-out RMSE: %.4f%n", result.holdoutRmse());
        }

        if (weightsOut != null) {
            try {
                new ObjectMapper().enable(SerializationFeature.INDENT_OUTPUT).writeValue(weightsOut.toFile(), weightsMap);
                System.out.println("\nFitted weights written to " + weightsOut);
            } catch (IOException e) {
                throw new org.gvi.core.exception.GviInputException("Could not write weights to " + weightsOut + ": " + e.getMessage(), e);
            }
        }
        return 0;
    }

    private int runPipeline() {
        Set<String> indexSet = new LinkedHashSet<>(Arrays.asList(indices));
        boolean generationTimeExplicit = generationTimeGivenExplicitly();
        double resolvedGenerationTime = resolveGenerationTimeDays(generationTimeExplicit);
        PipelineConfig config = new PipelineConfig(
                fasta, metadata, gff, codonUsage, codonUsageSpecies, incidence, referenceGc, referenceId,
                lower(indexSet), org.gvi.core.model.OrganismClass.parse(organismClass), pathogenId,
                generationTimeExplicit, resolvedGenerationTime, parseGdMethod(gdMethod),
                beta0, scaleFactor, highAccuracyMu, gammaAlpha, substitutionModel, bootstrapSupport, lsdMu, bdskyRe, mlDnds,
                weightsPath, parseGenomeType(genomeType), segments == null ? Map.of() : segments, gffOut,
                trimToCovered, minTrimmedColumns, auto, derivedOutputDir()
        );

        PipelineResult result = new GviPipeline().run(config);

        ReportWriter writer = new ReportWriter();
        writer.writeText(result, System.out, perSequence);
        if (jsonOut != null) {
            writer.writeJson(result, jsonOut, perSequence);
            System.out.println("\nJSON report written to " + jsonOut);
        }
        if (csvOut != null) {
            if (perSequence) {
                writer.writePerSequenceCsv(result, csvOut);
            } else {
                writer.writeCsv(result, csvOut);
            }
            System.out.println("CSV report written to " + csvOut);
        }
        if (figuresDir != null) {
            try {
                String label = figureLabel();
                var files = PublicationFigures.writeAll(result, figuresDir, label);
                System.out.println("\nFigures written to " + figuresDir + ":");
                for (var f : files) System.out.println("  " + f.getFileName());
            } catch (IOException e) {
                System.err.println("Could not write figures to " + figuresDir + ": " + e.getMessage());
            }
        }
        return 0;
    }

    /** True when the user actually typed --generation-time-days, rather than inheriting its default. */
    private boolean generationTimeGivenExplicitly() {
        picocli.CommandLine.Model.CommandSpec spec = picocli.CommandLine.Model.CommandSpec.forAnnotatedObject(this);
        try {
            return spec.commandLine() != null
                    && spec.commandLine().getParseResult().hasMatchedOption("--generation-time-days");
        } catch (RuntimeException e) {
            // No parse result available (direct construction in a test) -- treat as explicit, preserving old behaviour.
            return true;
        }
    }

    /**
     * Resolves the serial-interval mean. Precedence: an explicit --generation-time-days always
     * wins; otherwise --pathogen-id looks it up in the generation-time table, which fails loudly
     * rather than substituting a default. With neither, the historical default is carried forward
     * and {@link GviPipeline} warns about it.
     */
    private double resolveGenerationTimeDays(boolean explicit) {
        var table = generationTimesPath != null
                ? org.gvi.algorithms.re.GenerationTimeTable.fromFile(generationTimesPath)
                : org.gvi.algorithms.re.GenerationTimeTable.bundled();
        return org.gvi.algorithms.re.GenerationTimeTable.resolve(
                explicit, generationTimeDays, pathogenId, generationTimeDays, table);
    }

    /**
     * A filename stem for the figures that actually identifies the dataset.
     * <p>
     * The obvious choice -- the FASTA's own name -- collides constantly in practice, because a
     * corpus organised one-directory-per-pathogen tends to call every file {@code aligned.fasta};
     * two runs then silently overwrite each other's figures. So: prefer the JSON report's name,
     * which users give meaningfully, and otherwise fall back to the FASTA stem qualified by its
     * parent directory when that stem is a generic pipeline artefact name.
     */
    private String figureLabel() {
        if (jsonOut != null) {
            return jsonOut.getFileName().toString().replaceAll("\\.[^.]*$", "");
        }
        String stem = fasta.getFileName().toString().replaceAll("\\.[^.]*$", "");
        boolean generic = java.util.Set.of("aligned", "alignment", "aligned_clean", "input", "sequences", "raw")
                .contains(stem.toLowerCase(java.util.Locale.ROOT));
        java.nio.file.Path parent = fasta.toAbsolutePath().getParent();
        if (generic && parent != null) {
            return parent.getFileName().toString() + "_" + stem;
        }
        return stem;
    }

    /** Where --auto writes derived inputs: beside the JSON report, else beside the input FASTA. */
    private Path derivedOutputDir() {
        if (!auto) return null;
        if (jsonOut != null && jsonOut.getParent() != null) return jsonOut.getParent();
        if (fasta != null && fasta.getParent() != null) return fasta.getParent();
        return Path.of(".");
    }

    private Set<String> lower(Set<String> in) {
        Set<String> out = new LinkedHashSet<>();
        for (String s : in) out.add(s.toLowerCase());
        return out;
    }

    private GdMethod parseGdMethod(String s) {
        try {
            return GdMethod.valueOf(s.toUpperCase());
        } catch (IllegalArgumentException e) {
            throw new org.gvi.core.exception.GviInputException("Unknown --gd-method '" + s + "'; expected hamming, jukes_cantor, or kimura_2_parameter");
        }
    }

    private org.gvi.algorithms.mu.GenomeType parseGenomeType(String s) {
        if (s == null || s.isBlank()) {
            return org.gvi.algorithms.mu.GenomeType.UNSPECIFIED;
        }
        return switch (s.toLowerCase()) {
            case "rna" -> org.gvi.algorithms.mu.GenomeType.RNA;
            case "dna" -> org.gvi.algorithms.mu.GenomeType.DNA;
            default -> throw new org.gvi.core.exception.GviInputException("Unknown --genome-type '" + s + "'; expected rna or dna");
        };
    }
}
