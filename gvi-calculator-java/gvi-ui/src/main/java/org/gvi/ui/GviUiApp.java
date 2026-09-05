package org.gvi.ui;

import javafx.application.Application;
import javafx.application.Platform;
import javafx.concurrent.Task;
import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.Node;
import javafx.scene.Scene;
import javafx.scene.control.Control;
import javafx.scene.control.*;
import javafx.scene.control.cell.PropertyValueFactory;
import javafx.scene.layout.FlowPane;
import javafx.scene.layout.GridPane;
import javafx.scene.layout.HBox;
import javafx.scene.layout.Priority;
import javafx.scene.layout.Region;
import javafx.scene.layout.VBox;
import javafx.scene.input.TransferMode;
import javafx.stage.DirectoryChooser;
import javafx.stage.FileChooser;
import javafx.stage.Stage;
import org.gvi.algorithms.gd.GdMethod;
import org.gvi.cli.DatasetSummary;
import org.gvi.cli.GviPipeline;
import org.gvi.cli.PipelineConfig;
import org.gvi.cli.PipelineResult;
import org.gvi.cli.ReportWriter;
import org.gvi.composite.GviComponent;
import org.gvi.composite.IndexKey;
import org.gvi.core.exception.GviException;
import org.gvi.core.io.FastaReader;
import org.gvi.core.model.SequenceAlignment;
import org.gvi.core.util.GlobalExceptionHandler;
import org.gvi.selftest.SelfTestCase;
import org.gvi.selftest.SelfTestReport;
import org.gvi.selftest.SelfTestRunner;
import org.gvi.ui.cards.MetricCard;
import org.gvi.ui.chart.ContributionBars;
import org.gvi.ui.chart.GaugeView;
import org.gvi.ui.chart.TornadoBars;
import org.gvi.ui.nav.AppShell;
import org.gvi.ui.nav.NavGroup;
import org.gvi.ui.nav.NavRail;
import org.gvi.ui.nav.WorkflowStep;
import org.gvi.ui.pipeline.PipelineStagesView;
import org.gvi.ui.theme.Fonts;
import org.gvi.ui.theme.MetricDescriptions;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.EnumMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;

/**
 * Desktop GUI: the exact same pipeline gvi-cli runs (GviPipeline /
 * PipelineConfig / PipelineResult / ReportWriter, reused directly rather
 * than reimplemented) behind file pickers and a results table instead of
 * command-line flags. Computation always runs on a background
 * {@link Task} so a large alignment never freezes the UI thread, and every
 * error path shows an {@link Alert} instead of the application dying --
 * the same crash-containment philosophy as the CLI, adapted to a GUI.
 *
 * <p>Laid out as a workflow shell ({@link AppShell}): a left nav ({@link NavRail}) groups every real
 * input/option/result into the same 5-group taxonomy a bioinformatics platform brief asked for
 * (Project/Input/Analysis/GVI/Results), and every leaf page is built once at startup and cached, so
 * navigating away and back never loses in-progress field state. The underlying pipeline is still one
 * atomic {@code GviPipeline.run(config)} call -- it isn't internally staged -- so the Analysis group's
 * sub-items (Alignment/Phylogeny/Diversity/Evolution/Recombination/Temporal Signal) are honest
 * explainer pages pointing at where each concern is actually configured, not separate controllable steps.
 */
public final class GviUiApp extends Application {

    private static final Logger log = LoggerFactory.getLogger(GviUiApp.class);

    private static final String HELP_TEXT =
            "GVI Calculator v0.1.0 -- a standalone, offline Genomic Virulence Index tool.\n\n"
                    + "Computes 9 component indices from a pre-aligned viral FASTA -- evolutionary rate (mu), "
                    + "effective reproduction number (Re), nucleotide diversity (pi), mutation burden (MB), "
                    + "selection pressure (dN/dS), genetic distance (GD), codon adaptation (CAI), GC content "
                    + "deviation, and recombination (RI) -- then combines them into one weighted composite GVI(t) "
                    + "score, entirely locally: no data ever leaves this machine.\n\n"
                    + "Basic workflow:\n"
                    + "  1. Input > Sequence Data -- provide a pre-aligned FASTA (required); everything else is optional\n"
                    + "  2. GVI > GVI Metrics -- choose which indices to compute\n"
                    + "  3. GVI > GVI Calculation -- set options and click Run\n"
                    + "  4. Results -- Summary, Metric Dashboard, Visualizations, and Export Report\n\n"
                    + "New to this tool? Click \"Load Example Data\" on the Sequence Data page to load a small "
                    + "real working dataset (FASTA + metadata + incidence + custom weights) and see exactly how "
                    + "each input file should be formatted.\n\n"
                    + "Data-quality checks (e.g. an implausible evolutionary rate, a gap-heavy alignment) are "
                    + "reported automatically in Results > Summary's warnings and the Full Report.";

    private Stage stage;
    private AppShell appShell;
    private NavRail navRail;
    private final Map<String, Node> pageCache = new LinkedHashMap<>();

    private TextField fastaField, metadataField, gffField, codonUsageField, incidenceField, weightsField;
    private TextField referenceGcField, referenceIdField, generationTimeField, scaleFactorField, beta0Field;
    private TextField gammaAlphaField;
    private CheckBox highAccuracyMuCheckbox;
    private CheckBox bootstrapSupportCheckbox;
    private CheckBox lsdMuCheckbox;
    private CheckBox bdskyReCheckbox;
    private CheckBox mlDndsCheckbox;
    private CheckBox perSequenceCheckbox;
    private ComboBox<String> gdMethodCombo;
    private ComboBox<String> substitutionModelCombo;
    private ComboBox<String> codonUsageSpeciesCombo;
    private ComboBox<String> organismClassCombo;
    private ComboBox<String> genomeTypeCombo;
    private ComboBox<String> pathogenIdCombo;
    private CheckBox autoModeCheckbox;
    private CheckBox trimToCoveredCheckbox;
    private final Map<String, CheckBox> indexCheckboxes = new LinkedHashMap<>();
    private VBox datasetGviBox;
    private org.gvi.ui.chart.CompositionStrip compositionStrip;
    private QualityGatingView qualityGatingView;
    private TableView<SequenceRow> table;
    private TableView<SensitivityRow> sensitivityTable;
    private TextArea fullReportArea;
    private TextArea skippedWarningsArea;
    private Button runButton, stopButton, exportJsonButton, exportCsvButton, exportFiguresButton;
    private ProgressIndicator progress;
    private Label statusLabel;
    private PipelineResult lastResult;

    private GaugeView gaugeView;
    private Label seqCountValue, alignLengthValue, timeSpanValue, referenceValue, analysisStatusValue;
    private VBox topDriversBox;
    private FlowPane metricCardsFlow;
    private VBox fastaValidationBox;
    private PipelineStagesView pipelineStages;
    private Set<String> lastRequestedIndices = Set.of();
    private VBox contributionBarsBox, tornadoBarsBox;
    private Task<PipelineResult> runningTask;
    private Thread runningThread;
    private volatile boolean stopRequested;

    public static void main(String[] args) {
        launch(args);
    }

    @Override
    public void start(Stage primaryStage) {
        this.stage = primaryStage;

        Fonts.init();

        GlobalExceptionHandler.install((thread, throwable) -> Platform.runLater(() ->
                showAlert(Alert.AlertType.ERROR, "Unexpected internal error",
                        "Something went wrong that shouldn't have. See the log file for details.\n\n" + throwable)));

        initFields();

        List<WorkflowStep> steps = workflowSteps();
        for (WorkflowStep step : steps) {
            if (!"new_analysis".equals(step.id())) {
                pageCache.put(step.id(), buildPage(step.id()));
            }
        }

        navRail = new NavRail(steps, this::onNavSelect);
        ScrollPane navScroll = new ScrollPane(navRail);
        navScroll.setFitToWidth(true);
        navScroll.getStyleClass().add("gvi-nav-scroll");
        navScroll.setPrefWidth(230);
        navScroll.setMaxWidth(230);

        appShell = new AppShell(buildHeader(), navScroll, buildStatusBar());
        navRail.select("sequence_data");
        appShell.showPage(pageCache.get("sequence_data"));

        Scene scene = new Scene(appShell, 1300, 840);
        scene.getStylesheets().add(getClass().getResource("/theme/app.css").toExternalForm());
        primaryStage.setTitle("GVI Calculator");
        primaryStage.setScene(scene);
        primaryStage.show();
    }

    // ---- Workflow nav ----

    /**
     * The workflow, reduced to steps that actually do something.
     * <p>
     * Eight of the previous seventeen entries were signposts: pages whose only content was a
     * sentence explaining that the thing they named happens somewhere else ("Diversity ... is
     * computed automatically during GVI Calculation", "Alignment ... this tool does not perform
     * sequence alignment"). A navigation entry that exists to say "not here" costs the reader a
     * click to learn nothing, and makes the real steps harder to find. Three more were pages
     * holding one or two fields each, which belong beside the inputs they qualify rather than
     * behind their own click.
     */
    private List<WorkflowStep> workflowSteps() {
        return List.of(
                new WorkflowStep("sequence_data", "1. Input", NavGroup.WORKFLOW),
                new WorkflowStep("gvi_calculation", "2. Configure & Run", NavGroup.WORKFLOW),
                new WorkflowStep("summary", "Summary", NavGroup.RESULTS),
                new WorkflowStep("metric_dashboard", "Metrics", NavGroup.RESULTS),
                new WorkflowStep("visualizations", "Charts", NavGroup.RESULTS),
                new WorkflowStep("export_report", "Report & Export", NavGroup.RESULTS));
    }

    private void onNavSelect(WorkflowStep step) {
        if ("new_analysis".equals(step.id())) {
            clearAllFields();
            navRail.select("sequence_data");
            appShell.showPage(pageCache.get("sequence_data"));
            statusLabel.setText("New analysis started -- all fields cleared.");
            return;
        }
        appShell.showPage(pageCache.get(step.id()));
    }

    private Node buildPage(String id) {
        return switch (id) {
            case "sequence_data" -> buildSequenceDataPage();
            case "gvi_calculation" -> buildGviCalculationPage();
            case "summary" -> buildOverviewPage();
            case "metric_dashboard" -> buildMetricDashboardPage();
            case "visualizations" -> buildVisualizationsPage();
            case "export_report" -> buildResultsArea();
            default -> buildPlaceholderPage(id, "Not yet available in this build.");
        };
    }

    private Node buildHeader() {
        Label title = new Label("GVI CALCULATOR");
        title.getStyleClass().add("gvi-app-title");
        Label subtitle = new Label("Genomic Virulence Index Analysis Platform");
        subtitle.getStyleClass().add("gvi-app-subtitle");
        VBox titleBox = new VBox(2, title, subtitle);

        Button selfTestButton = new Button("Self-Test");
        selfTestButton.setOnAction(e -> runSelfTest());
        Button helpButton = new Button("Help");
        helpButton.setOnAction(e -> showAlert(Alert.AlertType.INFORMATION, "About GVI Calculator", HELP_TEXT));

        Region spacer = new Region();
        HBox.setHgrow(spacer, Priority.ALWAYS);

        HBox actions = new HBox(8, selfTestButton, helpButton);
        actions.setAlignment(Pos.CENTER_RIGHT);

        HBox bar = new HBox(titleBox, spacer, actions);
        bar.setAlignment(Pos.CENTER_LEFT);
        bar.setPadding(new Insets(12, 20, 12, 20));
        bar.getStyleClass().add("gvi-toolbar");
        return bar;
    }

    // ---- Field initialization (once, shared across pages) ----

    private void initFields() {
        fastaField = new TextField();
        metadataField = new TextField();
        gffField = new TextField();
        codonUsageField = new TextField();
        incidenceField = new TextField();
        weightsField = new TextField();
        referenceGcField = new TextField();
        referenceIdField = new TextField();
        generationTimeField = new TextField("5.0");
        scaleFactorField = new TextField("2.0");
        beta0Field = new TextField("1.0");
        gammaAlphaField = new TextField();
        highAccuracyMuCheckbox = new CheckBox();
        bootstrapSupportCheckbox = new CheckBox();
        lsdMuCheckbox = new CheckBox();
        bdskyReCheckbox = new CheckBox();
        mlDndsCheckbox = new CheckBox();
        perSequenceCheckbox = new CheckBox();
        perSequenceCheckbox.setOnAction(e -> {
            if (lastResult != null) refreshFullReportText();
        });
        gdMethodCombo = new ComboBox<>();
        gdMethodCombo.getItems().addAll("jukes_cantor", "kimura_2_parameter", "hamming");
        gdMethodCombo.setValue("jukes_cantor");
        substitutionModelCombo = new ComboBox<>();
        substitutionModelCombo.getItems().addAll("gtr", "jc69", "f81", "k80", "hky85", "tn93", "auto");
        substitutionModelCombo.setValue("gtr");
        codonUsageSpeciesCombo = new ComboBox<>();
        codonUsageSpeciesCombo.getItems().add("(none)");
        codonUsageSpeciesCombo.getItems().addAll(org.gvi.algorithms.cai.BundledCodonUsageTables.availableSpecies());
        codonUsageSpeciesCombo.setValue("(none)");

        // Organism class gates whether host-relative CAI is meaningful at all, and whether a GC
        // shift may be described as horizontal gene transfer. Nothing in an alignment reveals it,
        // and leaving it unset silently disables both checks -- so it belongs in the UI, not a TODO.
        organismClassCombo = new ComboBox<>();
        organismClassCombo.getItems().addAll("(unspecified)", "virus", "bacterium", "parasite");
        organismClassCombo.setValue("(unspecified)");

        // Genome type selects the scale mu is judged against: a DNA genome measured on the RNA
        // scale normalizes to nearly zero however fast it is for a DNA genome.
        genomeTypeCombo = new ComboBox<>();
        genomeTypeCombo.getItems().addAll("(unspecified)", "rna", "dna");
        genomeTypeCombo.setValue("(unspecified)");

        // Generation times come from the bundled table. Entries that cannot yield a usable value --
        // stubs, and organisms with no host-to-host transmission chain -- are listed with the reason
        // rather than hidden, so the absence is legible instead of looking like an oversight.
        pathogenIdCombo = new ComboBox<>();
        pathogenIdCombo.getItems().add("(use generation time below)");
        try {
            var table = org.gvi.algorithms.re.GenerationTimeTable.bundled();
            for (String id : table.usableIds()) {
                var e = table.find(id).orElseThrow();
                pathogenIdCombo.getItems().add(String.format("%s  (T=%.1f d, %s confidence)", id, e.tDays(), e.confidence()));
            }
        } catch (RuntimeException ignored) {
            // A missing bundled table must not stop the window opening; the manual field still works.
        }
        pathogenIdCombo.setValue("(use generation time below)");

        autoModeCheckbox = new CheckBox("Auto mode (aligned FASTA is the only required input)");
        autoModeCheckbox.setSelected(true);
        trimToCoveredCheckbox = new CheckBox("Always trim to shared coverage window");
    }

    private void clearAllFields() {
        fastaField.clear();
        metadataField.clear();
        gffField.clear();
        codonUsageField.clear();
        incidenceField.clear();
        weightsField.clear();
        referenceGcField.clear();
        referenceIdField.clear();
        generationTimeField.setText("5.0");
        scaleFactorField.setText("2.0");
        beta0Field.setText("1.0");
        gammaAlphaField.clear();
        highAccuracyMuCheckbox.setSelected(false);
        bootstrapSupportCheckbox.setSelected(false);
        lsdMuCheckbox.setSelected(false);
        bdskyReCheckbox.setSelected(false);
        mlDndsCheckbox.setSelected(false);
        perSequenceCheckbox.setSelected(false);
        gdMethodCombo.setValue("jukes_cantor");
        substitutionModelCombo.setValue("gtr");
        codonUsageSpeciesCombo.setValue("(none)");
        for (CheckBox cb : indexCheckboxes.values()) cb.setSelected(true);

        lastResult = null;
        if (compositionStrip != null) compositionStrip.clear();
        if (qualityGatingView != null) qualityGatingView.clear();
        table.getItems().clear();
        sensitivityTable.getItems().clear();
        fullReportArea.clear();
        skippedWarningsArea.clear();
        exportJsonButton.setDisable(true);
        exportCsvButton.setDisable(true);
        if (exportFiguresButton != null) exportFiguresButton.setDisable(true);

        gaugeView.setValue(null);
        seqCountValue.setText("--");
        alignLengthValue.setText("--");
        timeSpanValue.setText("--");
        referenceValue.setText("--");
        analysisStatusValue.setText("Not run yet");
        topDriversBox.getChildren().setAll(new Label("Run the pipeline to see which indices drive the GVI score."));
        refreshMetricCards();
        if (contributionBarsBox != null) contributionBarsBox.getChildren().setAll(ContributionBars.build(List.of()));
        if (tornadoBarsBox != null) tornadoBarsBox.getChildren().setAll(TornadoBars.build(List.of()));
        pipelineStages.reset();
    }

    // ---- Pages: Input group ----

    /** Muted explanatory line -- wraps, so long guidance stays readable in a narrow pane. */
    private Label hint(String text) {
        Label l = new Label(text);
        l.setWrapText(true);
        l.getStyleClass().add("hint-text");
        l.setStyle("-fx-text-fill: #6b7280; -fx-font-size: 11px;");
        return l;
    }

    private HBox labelledControl(String caption, Control control) {
        Label l = new Label(caption + ":");
        l.setMinWidth(230);
        HBox row = new HBox(8, l, control);
        row.setAlignment(Pos.CENTER_LEFT);
        return row;
    }

    private Node buildSequenceDataPage() {
        VBox dropZone = buildFastaDropZone();

        fastaValidationBox = new VBox(4);
        validateFastaFile(fastaField.getText());
        fastaField.textProperty().addListener((obs, oldVal, newVal) -> validateFastaFile(newVal));

        // The aligned FASTA is the only required input. Everything else this pane used to demand is
        // derived automatically in auto mode, so the labels now say so rather than presenting four
        // file pickers as if all four were needed.
        VBox required = new VBox(10);
        required.getChildren().add(dropZone);
        required.getChildren().add(fileRow("Aligned FASTA (the only required input)", fastaField,
                "FASTA files", "*.fasta", "*.fa", "*.fna"));
        required.getChildren().add(fastaValidationBox);
        required.getChildren().add(hint("In auto mode, collection dates are read from the FASTA headers, gene "
                + "boundaries are predicted, the reference is chosen, and the alignment is trimmed to its shared "
                + "coverage window only if it would otherwise fail the coverage check. Everything derived is written "
                + "beside your results so you can check it."));

        // Organism identity decides which indices are meaningful at all, and is the one thing that
        // cannot be read off an alignment -- so it belongs beside the input, not several panes away.
        VBox identity = new VBox(8);
        identity.getChildren().add(hint("These two cannot be inferred from sequence, and they change which indices "
                + "are valid: organism class decides whether host-relative CAI means anything, genome type sets the "
                + "scale mu is judged against. Leaving them unset disables both checks."));
        identity.getChildren().add(labelledControl("Organism class", organismClassCombo));
        identity.getChildren().add(labelledControl("Genome type", genomeTypeCombo));
        identity.getChildren().add(labelledControl("Pathogen (supplies generation time for Re)", pathogenIdCombo));
        identity.getChildren().add(autoModeCheckbox);
        identity.getChildren().add(trimToCoveredCheckbox);

        VBox optional = new VBox(10);
        optional.getChildren().add(hint("Supply any of these to override what auto mode derives."));
        optional.getChildren().add(fileRow("GFF3 annotation  (optional - ORFs predicted if omitted)", gffField,
                "GFF3 files", "*.gff", "*.gff3"));
        optional.getChildren().add(fileRow("Codon usage CSV  (optional - viruses only)", codonUsageField,
                "CSV files", "*.csv"));
        optional.getChildren().add(new Label("Or a bundled host codon usage table (Kazusa DB; used only if the CSV above is blank):"));
        optional.getChildren().add(codonUsageSpeciesCombo);
        optional.getChildren().add(fileRow("Incidence CSV  (optional - enables a true epidemiological Re)",
                incidenceField, "CSV files", "*.csv"));

        // Reference and metadata used to be separate nav pages holding two fields and one field
        // respectively. They qualify the input, so they sit with it.
        VBox refMeta = new VBox(8);
        refMeta.getChildren().add(hint("Collection dates are what make mu and Re computable at all. Without them "
                + "those two indices -- the most heavily weighted in the scheme -- are excluded."));
        refMeta.getChildren().add(fileRow("Metadata CSV  (sequence_id, collection_date, location, host)",
                metadataField, "CSV files", "*.csv"));
        refMeta.getChildren().add(labelledControl("Reference sequence id (optional)", referenceIdField));
        refMeta.getChildren().add(labelledControl("Reference GC % (optional)", referenceGcField));

        VBox content = new VBox(14,
                sectionCard("1. Aligned sequences", required),
                sectionCard("2. What this organism is", identity),
                sectionCard("3. Dates & reference", refMeta),
                sectionCard("4. Optional overrides", optional),
                buildExampleDataCard());
        return pageWrap(content);
    }

    private Node buildExampleDataCard() {
        Label desc = new Label("Not sure how your files should be formatted? Load a small real working dataset "
                + "(FASTA + metadata + incidence + custom weights) to see exactly what each input expects.");
        desc.setWrapText(true);
        desc.getStyleClass().add("gvi-placeholder-text");
        Button loadButton = new Button("Load Example Data");
        loadButton.setOnAction(e -> loadExampleData());
        VBox card = new VBox(8, desc, loadButton);
        card.getStyleClass().add("gvi-card");
        return card;
    }

    private void loadExampleData() {
        try {
            Path dir = Files.createTempDirectory("gvi-example-data");
            Path fasta = extractResource("/samples/samples.fasta", dir);
            Path metadata = extractResource("/samples/example_metadata.csv", dir);
            Path incidence = extractResource("/samples/example_incidence.csv", dir);
            Path weights = extractResource("/samples/example_weights.json", dir);

            fastaField.setText(fasta.toString());
            metadataField.setText(metadata.toString());
            incidenceField.setText(incidence.toString());
            weightsField.setText(weights.toString());

            navRail.select("sequence_data");
            appShell.showPage(pageCache.get("sequence_data"));
            statusLabel.setText("Example data loaded from " + dir + " -- 4 toy sequences "
                    + "(reference/query1/query2/query3) with real metadata/incidence/weights formats.");
            showExampleDataPreview(fasta, metadata, incidence, weights);
        } catch (IOException e) {
            showAlert(Alert.AlertType.ERROR, "Could not load example data", e.getMessage());
        }
    }

    private void showExampleDataPreview(Path fasta, Path metadata, Path incidence, Path weights) throws IOException {
        StringBuilder sb = new StringBuilder();
        appendFilePreview(sb, "FASTA -- sequence data (each record: >id, then the aligned sequence)", fasta);
        appendFilePreview(sb, "Metadata CSV -- sequence_id,collection_date,location,host", metadata);
        appendFilePreview(sb, "Incidence CSV -- date,new_cases (for the Re estimator)", incidence);
        appendFilePreview(sb, "Custom weights JSON -- index label -> weight override", weights);

        TextArea area = new TextArea(sb.toString());
        area.setEditable(false);
        area.setWrapText(false);
        area.setPrefSize(640, 460);
        area.getStyleClass().add("gvi-mono-area");

        Alert alert = new Alert(Alert.AlertType.INFORMATION);
        alert.setTitle("Example data -- how your files should look");
        alert.setHeaderText("The 4 files just loaded into this page's fields -- exact real content, so you can "
                + "match your own files to this format.");
        alert.getDialogPane().setContent(area);
        alert.getDialogPane().setPrefWidth(680);
        styleDialog(alert);
        alert.showAndWait();
    }

    private void appendFilePreview(StringBuilder sb, String title, Path file) throws IOException {
        sb.append("== ").append(title).append(" ==\n");
        sb.append(Files.readString(file));
        sb.append("\n\n");
    }

    private Path extractResource(String resourcePath, Path targetDir) throws IOException {
        String fileName = resourcePath.substring(resourcePath.lastIndexOf('/') + 1);
        Path target = targetDir.resolve(fileName);
        try (InputStream in = getClass().getResourceAsStream(resourcePath)) {
            if (in == null) throw new IOException("Bundled example resource not found: " + resourcePath);
            Files.copy(in, target, StandardCopyOption.REPLACE_EXISTING);
        }
        return target;
    }

    private VBox buildFastaDropZone() {
        Label icon = new Label("⇩");
        icon.getStyleClass().add("gvi-drop-icon");
        Label text = new Label("Drop aligned FASTA here, or use Browse below");
        text.getStyleClass().add("gvi-drop-text");
        VBox zone = new VBox(6, icon, text);
        zone.setAlignment(Pos.CENTER);
        zone.getStyleClass().add("gvi-drop-zone");

        zone.setOnDragOver(event -> {
            if (event.getDragboard().hasFiles()) {
                event.acceptTransferModes(TransferMode.COPY);
                zone.getStyleClass().add("gvi-drop-zone-active");
            }
            event.consume();
        });
        zone.setOnDragExited(event -> zone.getStyleClass().remove("gvi-drop-zone-active"));
        zone.setOnDragDropped(event -> {
            var db = event.getDragboard();
            boolean success = false;
            if (db.hasFiles() && !db.getFiles().isEmpty()) {
                fastaField.setText(db.getFiles().get(0).getAbsolutePath());
                success = true;
            }
            event.setDropCompleted(success);
            event.consume();
        });
        return zone;
    }

    private void validateFastaFile(String pathText) {
        if (fastaValidationBox == null) return;
        if (pathText == null || pathText.isBlank()) {
            fastaValidationBox.getChildren().setAll(placeholderLine("Drop or browse to a FASTA file to validate it here."));
            return;
        }
        Task<List<String>> task = new Task<>() {
            @Override
            protected List<String> call() {
                List<String> lines = new ArrayList<>();
                try {
                    FastaReader.Result result = FastaReader.read(Path.of(pathText));
                    lines.add((result.sequences().isEmpty() ? "✗ " : "✓ ") + result.sequences().size() + " sequence(s) detected");
                    for (var w : result.report().getWarnings()) {
                        lines.add("⚠ " + w);
                    }
                    if (!result.sequences().isEmpty()) {
                        SequenceAlignment alignment = SequenceAlignment.of(result.sequences());
                        lines.add(String.format(Locale.ROOT, "✓ Alignment length: %,d bp", alignment.length()));
                        lines.add("✓ Reference sequence identified: " + alignment.getReference().getId());
                    }
                } catch (GviException e) {
                    lines.add("✗ " + e.getMessage());
                } catch (Exception e) {
                    lines.add("✗ Could not read file: " + e.getMessage());
                }
                return lines;
            }
        };
        task.setOnSucceeded(e -> showFastaValidation(task.getValue()));
        task.setOnFailed(e -> showFastaValidation(List.of("✗ Unexpected error validating this file.")));
        new Thread(task, "gvi-fasta-validate").start();
    }

    private void showFastaValidation(List<String> lines) {
        List<Label> labels = new ArrayList<>();
        for (String line : lines) {
            Label l = new Label(line);
            l.getStyleClass().add(line.startsWith("✗") ? "gvi-validation-error"
                    : line.startsWith("⚠") ? "gvi-validation-warning" : "gvi-validation-ok");
            labels.add(l);
        }
        fastaValidationBox.getChildren().setAll(labels);
    }

    private Label placeholderLine(String text) {
        Label l = new Label(text);
        l.getStyleClass().add("gvi-placeholder-text");
        return l;
    }

    // ---- Pages: GVI group ----

    /** The nine index checkboxes, laid out three to a row. */
    private Node indexCheckboxPane() {
        HBox row1 = new HBox(10), row2 = new HBox(10), row3 = new HBox(10);
        for (String key : new String[]{"mu", "re", "pi"}) row1.getChildren().add(checkbox(key));
        for (String key : new String[]{"mb", "dnds", "gd"}) row2.getChildren().add(checkbox(key));
        for (String key : new String[]{"cai", "gc", "ri"}) row3.getChildren().add(checkbox(key));
        return new VBox(8, row1, row2, row3);
    }

    private Node buildGviCalculationPage() {
        GridPane grid = new GridPane();
        grid.setHgap(8);
        grid.setVgap(8);
        int r = 0;
        // Organism identity and input handling now live on the Sequence Data pane, beside the input
        // they describe, rather than here among the numerical-method settings.
        grid.addRow(r++, new Label("GD method:"), gdMethodCombo);
        grid.addRow(r++, new Label("Generation time (days, if no pathogen chosen):"), generationTimeField);
        grid.addRow(r++, new Label("High-accuracy mu (ML, slow):"), highAccuracyMuCheckbox);
        grid.addRow(r++, new Label("Bootstrap tree support (slow):"), bootstrapSupportCheckbox);
        grid.addRow(r++, new Label("Least-squares dating (LSD2-style):"), lsdMuCheckbox);
        grid.addRow(r++, new Label("BDSKY-style Re (ML, no incidence data):"), bdskyReCheckbox);
        grid.addRow(r++, new Label("ML dN/dS (codeml M0-equivalent, cross-check):"), mlDndsCheckbox);
        grid.addRow(r++, new Label("Include per-sequence detail in Full Report/exports:"), perSequenceCheckbox);
        grid.addRow(r++, new Label("Substitution model:"), substitutionModelCombo);
        grid.addRow(r++, new Label("Gamma alpha (optional):"), gammaAlphaField);

        runButton = new Button("Run");
        runButton.setDefaultButton(true);
        runButton.setMaxWidth(Double.MAX_VALUE);
        runButton.setOnAction(e -> runPipeline());
        HBox.setHgrow(runButton, Priority.ALWAYS);

        stopButton = new Button("Stop");
        stopButton.getStyleClass().add("gvi-stop-button");
        stopButton.setDisable(true);
        stopButton.setOnAction(e -> stopPipeline());

        HBox runRow = new HBox(8, runButton, stopButton);

        pipelineStages = new PipelineStagesView();

        // Index selection and custom weights each had their own nav page for a checkbox row and a
        // single file field. They are part of configuring this run, so they are part of this page.
        VBox indices = new VBox(8);
        indices.getChildren().add(hint("Leave all ticked unless you deliberately want a subset. An index with no "
                + "input, or one that fails quality gating, is dropped and reported either way."));
        indices.getChildren().add(indexCheckboxPane());

        VBox weights = new VBox(8);
        weights.getChildren().add(hint("Default weights are the specification's midpoints -- asserted, not fitted. "
                + "Supply your own here, or fit them against observed outcomes with the CLI's --calibrate."));
        weights.getChildren().add(fileRow("Custom composite weights JSON (optional)", weightsField,
                "JSON files", "*.json"));

        VBox content = new VBox(14,
                sectionCard("Indices to compute", indices),
                sectionCard("Method options", grid),
                sectionCard("Weights", weights),
                runRow,
                sectionCard("Pipeline", pipelineStages));
        return pageWrap(content);
    }

    // ---- Pages: Results group ----

    /** "Summary" -- the landing dashboard: GVI gauge, dataset summary, and top-contributing indices. */
    private Node buildOverviewPage() {
        gaugeView = new GaugeView();
        gaugeView.setValue(null);

        Button breakdownButton = new Button("How was this score calculated?");
        breakdownButton.getStyleClass().add("gvi-link-button");
        breakdownButton.setOnAction(e -> showScoreBreakdown());

        VBox gaugeBox = new VBox(10, gaugeView, breakdownButton);
        gaugeBox.setAlignment(Pos.CENTER);

        seqCountValue = new Label("--");
        alignLengthValue = new Label("--");
        timeSpanValue = new Label("--");
        referenceValue = new Label("--");
        analysisStatusValue = new Label("Not run yet");
        for (Label l : List.of(seqCountValue, alignLengthValue, timeSpanValue, referenceValue, analysisStatusValue)) {
            l.getStyleClass().add("gvi-summary-value");
        }

        GridPane summaryGrid = new GridPane();
        summaryGrid.setHgap(16);
        summaryGrid.setVgap(10);
        int r = 0;
        summaryGrid.addRow(r++, summaryLabel("Sequences"), seqCountValue);
        summaryGrid.addRow(r++, summaryLabel("Alignment"), alignLengthValue);
        summaryGrid.addRow(r++, summaryLabel("Time span"), timeSpanValue);
        summaryGrid.addRow(r++, summaryLabel("Reference"), referenceValue);
        summaryGrid.addRow(r++, summaryLabel("Analysis"), analysisStatusValue);

        HBox topRow = new HBox(30, gaugeBox, sectionCard("Dataset summary", summaryGrid));
        topRow.setAlignment(Pos.TOP_LEFT);

        topDriversBox = new VBox(6);
        topDriversBox.getChildren().add(new Label("Run the pipeline to see which indices drive the GVI score."));

        VBox content = new VBox(16, topRow, sectionCard("Top contributing indices", topDriversBox));
        return pageWrap(content, 900);
    }

    private Label summaryLabel(String text) {
        Label l = new Label(text);
        l.getStyleClass().add("gvi-summary-label");
        return l;
    }

    /** "Metric Dashboard" -- one card per index, real value or an honest not-available/not-run state. */
    private Node buildMetricDashboardPage() {
        metricCardsFlow = new FlowPane(14, 14);
        refreshMetricCards();
        VBox content = new VBox(10, sectionLabel("GVI Metrics -- calculated values"), metricCardsFlow);
        return pageWrap(content, 1120);
    }

    private void refreshMetricCards() {
        if (metricCardsFlow == null) return;
        metricCardsFlow.getChildren().clear();
        Map<IndexKey, GviComponent> componentsByKey = new EnumMap<>(IndexKey.class);
        Set<IndexKey> excludedKeys = new HashSet<>();
        if (lastResult != null && lastResult.datasetGvi() != null) {
            for (GviComponent c : lastResult.datasetGvi().components()) componentsByKey.put(c.key(), c);
            excludedKeys.addAll(lastResult.datasetGvi().excludedIndices());
        }
        for (IndexKey key : IndexKey.values()) {
            GviComponent comp = componentsByKey.get(key);
            boolean excluded = excludedKeys.contains(key);
            Runnable onClick = comp != null ? () -> showMetricDetail(key, comp) : null;
            VBox card = MetricCard.build(key, comp, excluded, onClick);
            card.setPrefWidth(330);
            card.setMaxWidth(330);
            metricCardsFlow.getChildren().add(card);
        }
    }

    private void refreshOverview(PipelineResult result) {
        if (gaugeView == null) return;
        Double gvi = result.datasetGvi() != null ? result.datasetGvi().gvi() : null;
        gaugeView.setValue(gvi);

        DatasetSummary summary = result.datasetSummary();
        if (summary != null) {
            seqCountValue.setText(String.valueOf(summary.sequenceCount()));
            alignLengthValue.setText(String.format(Locale.ROOT, "%,d bp", summary.alignmentLengthBp()));
            if (summary.earliestCollectionDate() != null && summary.latestCollectionDate() != null) {
                timeSpanValue.setText(summary.earliestCollectionDate().getYear() + "–" + summary.latestCollectionDate().getYear());
            } else {
                timeSpanValue.setText("Not available (no collection dates)");
            }
            referenceValue.setText(summary.referenceId());
        }
        analysisStatusValue.setText(result.datasetGvi() != null ? "Completed" : "Completed with errors -- see Skipped / Warnings");

        topDriversBox.getChildren().clear();
        if (result.datasetGvi() != null) {
            List<GviComponent> sorted = new ArrayList<>(result.datasetGvi().components());
            sorted.sort((a, b) -> Double.compare(b.contribution(), a.contribution()));
            for (GviComponent c : sorted) {
                Label row = new Label(String.format(Locale.ROOT, "%-28s %.4f",
                        MetricCard.displayName(c.key()), c.contribution()));
                row.getStyleClass().add("gvi-driver-row");
                topDriversBox.getChildren().add(row);
            }
        } else {
            topDriversBox.getChildren().add(new Label("Whole-file GVI could not be computed -- see Skipped / Warnings."));
        }
    }

    private void showScoreBreakdown() {
        if (lastResult == null || lastResult.datasetGvi() == null) {
            showAlert(Alert.AlertType.INFORMATION, "How was this score calculated?",
                    "Run the pipeline first to see a breakdown of the whole-file GVI.");
            return;
        }
        List<GviComponent> components = new ArrayList<>(lastResult.datasetGvi().components());
        components.sort((a, b) -> Double.compare(b.contribution(), a.contribution()));

        StringBuilder sb = new StringBuilder("GVI SCORE\n");
        sb.append("--------------------------------\n");
        for (GviComponent c : components) {
            sb.append(String.format(Locale.ROOT, "%-28s %8.4f  (weight %.2f)%n",
                    MetricCard.displayName(c.key()), c.contribution(), c.effectiveWeight()));
        }
        sb.append("--------------------------------\n");
        sb.append(String.format(Locale.ROOT, "%-28s %8.4f%n", "Final GVI", lastResult.datasetGvi().gvi()));
        if (!lastResult.datasetGvi().excludedIndices().isEmpty()) {
            sb.append("\nExcluded (no data): ").append(lastResult.datasetGvi().excludedIndices().stream()
                    .map(IndexKey::label).reduce((a, b) -> a + ", " + b).orElse(""));
        }

        TextArea area = new TextArea(sb.toString());
        area.setEditable(false);
        area.setPrefSize(500, 340);
        area.getStyleClass().add("gvi-mono-area");

        Alert alert = new Alert(Alert.AlertType.INFORMATION);
        alert.setTitle("How was this score calculated?");
        alert.setHeaderText("Composite GVI breakdown -- real per-index contributions from this run");
        alert.getDialogPane().setContent(area);
        styleDialog(alert);
        alert.showAndWait();
    }

    private void showMetricDetail(IndexKey key, GviComponent component) {
        String body = MetricDescriptions.describe(key)
                + String.format(Locale.ROOT, "%n%nCalculated value: %s = %.6f%nComputed from the selected dataset.",
                        key.label(), component.rawValue())
                + "\n\nInterpretation: depends on the evolutionary model, timescale, sampling design, and viral "
                + "system -- provided for expert review, not as a standalone conclusion.";
        showAlert(Alert.AlertType.INFORMATION, MetricCard.displayName(key), body);
    }

    // ---- Shared page helpers ----

    private Node pageWrap(Node content) {
        return pageWrap(content, 760);
    }

    private Node pageWrap(Node content, double maxWidth) {
        VBox box = new VBox(content);
        box.setPadding(new Insets(20));
        box.setMaxWidth(maxWidth);
        ScrollPane scroll = new ScrollPane(box);
        scroll.setFitToWidth(true);
        scroll.getStyleClass().add("gvi-form-scroll");
        return scroll;
    }

    private Node buildPlaceholderPage(String title, String message) {
        Label msg = new Label(message);
        msg.setWrapText(true);
        msg.getStyleClass().add("gvi-placeholder-text");
        VBox card = new VBox(8, msg);
        return pageWrap(sectionCard(title, card));
    }

    private CheckBox checkbox(String key) {
        CheckBox cb = new CheckBox(key);
        cb.setSelected(true);
        indexCheckboxes.put(key, cb);
        return cb;
    }

    private Label sectionLabel(String text) {
        Label l = new Label(text);
        l.getStyleClass().add("gvi-section-label");
        return l;
    }

    private VBox sectionCard(String title, Node content) {
        VBox card = new VBox(8, sectionLabel(title), content);
        card.getStyleClass().add("gvi-card");
        return card;
    }

    private HBox fileRow(String label, TextField field, String filterName, String... extensions) {
        field.setEditable(false);
        field.setPrefWidth(300);
        Button browse = new Button("Browse…");
        browse.setOnAction(e -> {
            FileChooser fc = new FileChooser();
            fc.setTitle("Select " + label);
            fc.getExtensionFilters().add(new FileChooser.ExtensionFilter(filterName, extensions));
            fc.getExtensionFilters().add(new FileChooser.ExtensionFilter("All files", "*.*"));
            File f = fc.showOpenDialog(stage);
            if (f != null) field.setText(f.getAbsolutePath());
        });
        Button clear = new Button("Clear");
        clear.setOnAction(e -> field.clear());
        VBox labelAndField = new VBox(2, new Label(label), new HBox(6, field, browse, clear));
        return new HBox(labelAndField);
    }

    /** "Visualizations" -- contribution-decomposition bars and the sensitivity tornado chart. */
    private Node buildVisualizationsPage() {
        contributionBarsBox = new VBox(ContributionBars.build(List.of()));
        tornadoBarsBox = new VBox(TornadoBars.build(List.of()));
        VBox content = new VBox(20,
                sectionCard("Contribution decomposition -- each index's share of the whole-file GVI", contributionBarsBox),
                sectionCard("Sensitivity analysis -- how far the GVI moves at ±20% weight", tornadoBarsBox));
        return pageWrap(content, 900);
    }

    private void refreshVisualizations(PipelineResult result) {
        if (contributionBarsBox == null) return;
        List<GviComponent> components = result.datasetGvi() != null ? result.datasetGvi().components() : List.of();
        contributionBarsBox.getChildren().setAll(ContributionBars.build(components));
        tornadoBarsBox.getChildren().setAll(TornadoBars.build(result.sensitivity()));
    }

    // ---- Results (shared by Summary / Metric Dashboard / Visualizations / Export Report) ----

    private VBox buildResultsArea() {
        // The composition strip replaces a bare number: it shows how much of the weighting scheme
        // the score actually rests on, which a GVI value alone cannot convey.
        compositionStrip = new org.gvi.ui.chart.CompositionStrip();
        datasetGviBox = new VBox(compositionStrip);
        datasetGviBox.getStyleClass().add("gvi-banner");

        table = new TableView<>();
        table.getColumns().addAll(
                column("Sequence", "sequenceId", 130),
                column("GVI (this seq. only)", "gvi", 130),
                column("mu", "mu", 80), column("Re", "re", 80), column("pi", "pi", 80),
                column("MB", "mb", 70), column("dN/dS", "dnds", 80), column("GD", "gd", 80),
                column("CAI", "cai", 80), column("GC_Dev", "gc", 80), column("RI", "ri", 80));
        table.setPlaceholder(new Label("Run the pipeline to see results here."));

        sensitivityTable = new TableView<>();
        sensitivityTable.getColumns().addAll(
                sensitivityColumn("Index", "index", 100),
                sensitivityColumn("Base GVI", "baseGvi", 90),
                sensitivityColumn("GVI at -20% weight", "lowWeightGvi", 130),
                sensitivityColumn("GVI at +20% weight", "highWeightGvi", 130),
                sensitivityColumn("Spread", "spread", 90));
        sensitivityTable.setPlaceholder(new Label(
                "Run the pipeline to see how much each index's weight moves the whole-file GVI (most sensitive first)."));

        fullReportArea = new TextArea();
        fullReportArea.setEditable(false);
        fullReportArea.getStyleClass().add("gvi-mono-area");

        skippedWarningsArea = new TextArea();
        skippedWarningsArea.setEditable(false);
        skippedWarningsArea.getStyleClass().add("gvi-mono-area");

        // "Full Report" first so the whole-file result (the headline answer, same as the CLI's default
        // output) is what's shown by default -- TabPane auto-selects whichever tab is added first, and
        // showing "Per-Sequence Detail" by default here would contradict the CLI defaulting to whole-file-only.
        TabPane tabs = new TabPane();
        tabs.getTabs().add(new Tab("Full Report", fullReportArea));
        tabs.getTabs().add(new Tab("Sensitivity Analysis", sensitivityTable));
        tabs.getTabs().add(new Tab("Per-Sequence Detail", table));
        qualityGatingView = new QualityGatingView();
        tabs.getTabs().add(new Tab("Quality Gating", qualityGatingView));
        tabs.getTabs().forEach(t -> t.setClosable(false));

        exportJsonButton = new Button("Export JSON…");
        exportJsonButton.setDisable(true);
        exportJsonButton.setOnAction(e -> exportJson());
        exportCsvButton = new Button("Export CSV (whole-file)…");
        exportCsvButton.setDisable(true);
        exportCsvButton.setOnAction(e -> exportCsv());
        exportFiguresButton = new Button("Export Figures (SVG)…");
        exportFiguresButton.setDisable(true);
        exportFiguresButton.setOnAction(e -> exportFigures());
        HBox exportBar = new HBox(10, exportJsonButton, exportCsvButton, exportFiguresButton);
        exportBar.setPadding(new Insets(8));

        VBox box = new VBox(datasetGviBox, tabs, exportBar);
        VBox.setVgrow(tabs, Priority.ALWAYS);
        return box;
    }

    /**
     * Writes the vector figures for the current result. Vector rather than a screenshot of the
     * on-screen charts: a figure destined for print has to survive being scaled to a journal
     * column, and a raster capture of a JavaFX node does not.
     */
    private void exportFigures() {
        if (lastResult == null) return;
        DirectoryChooser chooser = new DirectoryChooser();
        chooser.setTitle("Choose a folder for the figures");
        java.io.File dir = chooser.showDialog(stage);
        if (dir == null) return;
        try {
            String label = fastaField.getText() == null || fastaField.getText().isBlank()
                    ? "gvi" : Path.of(fastaField.getText()).getFileName().toString().replaceAll("\\.[^.]*$", "");
            var files = org.gvi.cli.PublicationFigures.writeAll(lastResult, dir.toPath(), label);
            StringBuilder sb = new StringBuilder("Wrote " + files.size() + " figure(s):\n\n");
            for (var f : files) sb.append("  ").append(f.getFileName()).append('\n');
            sb.append("\nVector SVG at single-column width, colour-vision-safe, and legible in greyscale. "
                    + "Excluded indices are drawn hatched rather than omitted.");
            showAlert(Alert.AlertType.INFORMATION, "Figures exported", sb.toString());
        } catch (Exception ex) {
            showAlert(Alert.AlertType.ERROR, "Could not write figures", String.valueOf(ex.getMessage()));
        }
    }

    private TableColumn<SequenceRow, String> column(String title, String property, int width) {
        TableColumn<SequenceRow, String> col = new TableColumn<>(title);
        col.setCellValueFactory(new PropertyValueFactory<>(property));
        col.setPrefWidth(width);
        return col;
    }

    private TableColumn<SensitivityRow, String> sensitivityColumn(String title, String property, int width) {
        TableColumn<SensitivityRow, String> col = new TableColumn<>(title);
        col.setCellValueFactory(new PropertyValueFactory<>(property));
        col.setPrefWidth(width);
        return col;
    }

    private HBox buildStatusBar() {
        progress = new ProgressIndicator();
        progress.setVisible(false);
        progress.setMaxSize(20, 20);
        statusLabel = new Label("Ready.");
        HBox bar = new HBox(10, progress, statusLabel);
        bar.setPadding(new Insets(8));
        bar.setAlignment(Pos.CENTER_LEFT);
        bar.getStyleClass().add("gvi-status-bar");
        return bar;
    }

    // ---- Pipeline execution ----

    private void runPipeline() {
        if (fastaField.getText().isBlank()) {
            showAlert(Alert.AlertType.WARNING, "FASTA required", "Please select a FASTA alignment file before running.");
            return;
        }
        Set<String> indices = new LinkedHashSet<>();
        for (var e : indexCheckboxes.entrySet()) if (e.getValue().isSelected()) indices.add(e.getKey());
        if (indices.isEmpty()) {
            showAlert(Alert.AlertType.WARNING, "No indices selected", "Select at least one index to compute.");
            return;
        }

        final PipelineConfig config;
        try {
            config = new PipelineConfig(
                    Path.of(fastaField.getText()),
                    pathOrNull(metadataField), pathOrNull(gffField), pathOrNull(codonUsageField),
                    "(none)".equals(codonUsageSpeciesCombo.getValue()) ? null : codonUsageSpeciesCombo.getValue(),
                    pathOrNull(incidenceField),
                    doubleOrNull(referenceGcField, "Reference GC %"),
                    referenceIdField.getText().isBlank() ? null : referenceIdField.getText().trim(),
                    indices,
                    parseOrganismClassFromUi(),
                    selectedPathogenId(),
                    true,
                    resolveGenerationTimeFromUi(),
                    parseGdMethod(gdMethodCombo.getValue()),
                    parseDouble(beta0Field, "Beta0", 1.0),
                    parseDouble(scaleFactorField, "Beta scale factor", 2.0),
                    highAccuracyMuCheckbox.isSelected(),
                    doubleOrNull(gammaAlphaField, "Gamma alpha"),
                    substitutionModelCombo.getValue(),
                    bootstrapSupportCheckbox.isSelected(),
                    lsdMuCheckbox.isSelected(),
                    bdskyReCheckbox.isSelected(),
                    mlDndsCheckbox.isSelected(),
                    pathOrNull(weightsField),
                    parseGenomeTypeFromUi(),
                    // TODO: no segment/reassortment controls in the UI yet (see GviCli's --segment).
                    java.util.Map.of(),
                    // TODO: no --gff-out control in the UI yet.
                    null,
                    trimToCoveredCheckbox.isSelected(),
                    100,
                    autoModeCheckbox.isSelected(),
                    autoModeCheckbox.isSelected() ? Path.of(fastaField.getText()).toAbsolutePath().getParent() : null
            );
        } catch (IllegalArgumentException e) {
            showAlert(Alert.AlertType.WARNING, "Invalid input", e.getMessage());
            return;
        }

        stopRequested = false;
        lastRequestedIndices = indices;
        setRunning(true);
        pipelineStages.setRunning(indices);
        Task<PipelineResult> task = new Task<>() {
            @Override
            protected PipelineResult call() {
                return new GviPipeline().run(config);
            }
        };
        task.setOnSucceeded(e -> {
            setRunning(false);
            if (stopRequested) {
                statusLabel.setText("Stopped run finished in the background -- result discarded.");
                return;
            }
            lastResult = task.getValue();
            populateResults(lastResult);
            pipelineStages.applyResult(lastResult, lastRequestedIndices);
            statusLabel.setText("Done: " + lastResult.perSequenceIndices().size() + " sequence(s) processed.");
            exportJsonButton.setDisable(false);
            exportCsvButton.setDisable(false);
        if (exportFiguresButton != null) exportFiguresButton.setDisable(false);
        });
        task.setOnFailed(e -> {
            setRunning(false);
            if (stopRequested) {
                statusLabel.setText("Stopped.");
                pipelineStages.reset();
                return;
            }
            pipelineStages.applyFailure();
            Throwable ex = task.getException();
            if (ex instanceof GviException) {
                statusLabel.setText("Failed: " + ex.getMessage());
                showAlert(Alert.AlertType.ERROR, "Could not compute", ex.getMessage());
            } else {
                statusLabel.setText("Unexpected error -- see log.");
                log.error("Pipeline run failed unexpectedly", ex);
                showAlert(Alert.AlertType.ERROR, "Unexpected internal error", String.valueOf(ex));
            }
        });
        task.setOnCancelled(e -> {
            setRunning(false);
            statusLabel.setText("Stopped.");
            pipelineStages.reset();
        });
        runningTask = task;
        runningThread = new Thread(task, "gvi-pipeline-run");
        runningThread.setDaemon(true);
        runningThread.start();
    }

    private void stopPipeline() {
        if (runningTask == null || !runningTask.isRunning()) return;
        stopRequested = true;
        runningTask.cancel(true);
        if (runningThread != null) {
            runningThread.interrupt();
        }
        setRunning(false);
        statusLabel.setText("Stop requested -- the computation has no internal cancellation points, so it may "
                + "keep using CPU briefly in the background, but its result will be discarded and the UI is free again.");
    }

    private void setRunning(boolean running) {
        runButton.setDisable(running);
        stopButton.setDisable(!running);
        progress.setVisible(running);
        statusLabel.setText(running ? "Running…" : statusLabel.getText());
    }

    private void populateResults(PipelineResult result) {
        compositionStrip.show(result.datasetGvi(), result.perSequenceIndices().size());
        qualityGatingView.show(result);

        table.getItems().setAll(SequenceRow.fromResult(result));
        sensitivityTable.getItems().setAll(SensitivityRow.fromResult(result));
        refreshOverview(result);
        refreshMetricCards();
        refreshVisualizations(result);

        refreshFullReportText();

        StringBuilder sw = new StringBuilder();
        sw.append("-- Skipped --\n");
        result.skipped().forEach(s -> sw.append("  - ").append(s).append('\n'));
        sw.append("\n-- Warnings --\n");
        result.warnings().forEach(s -> sw.append("  - ").append(s).append('\n'));
        skippedWarningsArea.setText(sw.toString());
    }

    private void refreshFullReportText() {
        ByteArrayOutputStream buffer = new ByteArrayOutputStream();
        new ReportWriter().writeText(lastResult, new PrintStream(buffer, true, StandardCharsets.UTF_8), perSequenceCheckbox.isSelected());
        fullReportArea.setText(buffer.toString(StandardCharsets.UTF_8));
    }

    private void exportJson() {
        FileChooser fc = new FileChooser();
        fc.setTitle("Export JSON report");
        fc.setInitialFileName("gvi_report.json");
        fc.getExtensionFilters().add(new FileChooser.ExtensionFilter("JSON", "*.json"));
        File f = fc.showSaveDialog(stage);
        if (f == null) return;
        try {
            new ReportWriter().writeJson(lastResult, f.toPath(), perSequenceCheckbox.isSelected());
            statusLabel.setText("JSON report written to " + f.getAbsolutePath());
        } catch (GviException e) {
            showAlert(Alert.AlertType.ERROR, "Export failed", e.getMessage());
        }
    }

    private void exportCsv() {
        FileChooser fc = new FileChooser();
        fc.setTitle("Export CSV summary");
        fc.setInitialFileName("gvi_report.csv");
        fc.getExtensionFilters().add(new FileChooser.ExtensionFilter("CSV", "*.csv"));
        File f = fc.showSaveDialog(stage);
        if (f == null) return;
        try {
            if (perSequenceCheckbox.isSelected()) {
                new ReportWriter().writePerSequenceCsv(lastResult, f.toPath());
            } else {
                new ReportWriter().writeCsv(lastResult, f.toPath());
            }
            statusLabel.setText("CSV report written to " + f.getAbsolutePath());
        } catch (GviException e) {
            showAlert(Alert.AlertType.ERROR, "Export failed", e.getMessage());
        }
    }

    private void runSelfTest() {
        Task<SelfTestReport> task = new Task<>() {
            @Override
            protected SelfTestReport call() {
                return new SelfTestRunner().runAll();
            }
        };
        task.setOnSucceeded(e -> showSelfTestResults(task.getValue()));
        task.setOnFailed(e -> showAlert(Alert.AlertType.ERROR, "Self-test error", String.valueOf(task.getException())));
        new Thread(task, "gvi-self-test").start();
    }

    private void showSelfTestResults(SelfTestReport report) {
        StringBuilder sb = new StringBuilder();
        for (SelfTestCase c : report.cases()) {
            sb.append(c.passed() ? "[PASS] " : "[FAIL] ").append(c.name()).append('\n');
            if (!c.passed()) sb.append("        ").append(c.detail()).append('\n');
        }
        TextArea area = new TextArea(sb.toString());
        area.setEditable(false);
        area.setPrefSize(560, 400);
        area.getStyleClass().add("gvi-mono-area");

        Alert alert = new Alert(report.allPassed() ? Alert.AlertType.INFORMATION : Alert.AlertType.WARNING);
        alert.setTitle("Self-Test Results");
        alert.setHeaderText(report.allPassed()
                ? "All self-tests passed -- this installation is computing correctly."
                : "SELF-TEST FAILURES DETECTED -- do not trust results until resolved.");
        alert.getDialogPane().setContent(area);
        styleDialog(alert);
        alert.showAndWait();
    }

    private void showAlert(Alert.AlertType type, String title, String message) {
        Alert alert = new Alert(type);
        alert.setTitle(title);
        alert.setHeaderText(null);
        alert.setContentText(message);
        alert.getDialogPane().setPrefWidth(480);
        styleDialog(alert);
        alert.showAndWait();
    }

    private void styleDialog(Alert alert) {
        alert.getDialogPane().getStylesheets().add(getClass().getResource("/theme/app.css").toExternalForm());
    }

    private Path pathOrNull(TextField field) {
        return field.getText().isBlank() ? null : Path.of(field.getText());
    }

    private Double doubleOrNull(TextField field, String fieldName) {
        if (field.getText().isBlank()) return null;
        try {
            return Double.parseDouble(field.getText().trim());
        } catch (NumberFormatException e) {
            throw new IllegalArgumentException(fieldName + " must be a number, got '" + field.getText() + "'");
        }
    }

    private double parseDouble(TextField field, String fieldName, double fallback) {
        if (field.getText().isBlank()) return fallback;
        try {
            return Double.parseDouble(field.getText().trim());
        } catch (NumberFormatException e) {
            throw new IllegalArgumentException(fieldName + " must be a number, got '" + field.getText() + "'");
        }
    }

    private GdMethod parseGdMethod(String s) {
        return GdMethod.valueOf(s.toUpperCase());
    }

    private org.gvi.core.model.OrganismClass parseOrganismClassFromUi() {
        String v = organismClassCombo.getValue();
        return (v == null || v.startsWith("(")) ? org.gvi.core.model.OrganismClass.UNSPECIFIED
                : org.gvi.core.model.OrganismClass.parse(v);
    }

    private org.gvi.algorithms.mu.GenomeType parseGenomeTypeFromUi() {
        String v = genomeTypeCombo.getValue();
        if (v == null || v.startsWith("(")) return org.gvi.algorithms.mu.GenomeType.UNSPECIFIED;
        return "dna".equals(v) ? org.gvi.algorithms.mu.GenomeType.DNA : org.gvi.algorithms.mu.GenomeType.RNA;
    }

    /** The bare pathogen id from the combo's descriptive label, or null when none was chosen. */
    private String selectedPathogenId() {
        String v = pathogenIdCombo.getValue();
        if (v == null || v.startsWith("(")) return null;
        int space = v.indexOf(' ');
        return space < 0 ? v : v.substring(0, space);
    }

    /**
     * A chosen pathogen supplies its own generation time from the bundled table; the manual field is
     * only consulted when none was chosen. This keeps the UI from quietly applying a 5-day default to
     * an organism whose real interval is nothing like 5 days -- the failure that made Re read ~1.0
     * for every pathogen in the corpus.
     */
    private double resolveGenerationTimeFromUi() {
        String id = selectedPathogenId();
        if (id != null) {
            try {
                return org.gvi.algorithms.re.GenerationTimeTable.bundled().requireGenerationTimeDays(id);
            } catch (RuntimeException e) {
                throw new IllegalArgumentException("Generation time for '" + id + "': " + e.getMessage());
            }
        }
        return parseDouble(generationTimeField, "Generation time", 5.0);
    }
}
