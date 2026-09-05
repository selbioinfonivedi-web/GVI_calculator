package org.gvi.ui.pipeline;

import javafx.geometry.Pos;
import javafx.scene.control.Label;
import javafx.scene.layout.HBox;
import javafx.scene.layout.Priority;
import javafx.scene.layout.Region;
import javafx.scene.layout.VBox;
import org.gvi.cli.PipelineResult;
import org.gvi.composite.GviComponent;
import org.gvi.composite.IndexKey;

import java.util.EnumMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Shows the analysis as a stage list, per the product brief's pipeline-visualization request.
 *
 * <p>Honest caveat baked into the design: {@code GviPipeline.run(config)} is one atomic call, not
 * internally staged, so this view cannot show true live per-stage progress. Before a run, every
 * applicable stage shows Pending. During a run, every applicable stage shows a single shared Running
 * indicator (there is no way to know which stage is actually executing at any instant -- this is a
 * deliberate, documented limitation, not an oversight). After a run, each stage's status is derived
 * from the REAL {@link PipelineResult}: did its index end up in {@code datasetGvi().components()}
 * (Completed), in {@code datasetGvi().excludedIndices()} (Warning -- attempted, no data), or was it
 * never selected this run at all (Not selected) -- so the after-the-fact readout is accurate, even
 * though the during-the-run animation necessarily isn't truly per-stage.
 */
public final class PipelineStagesView extends VBox {

    public enum Status { PENDING, RUNNING, COMPLETED, WARNING, FAILED, NOT_SELECTED }

    private static final List<PipelineStage> STAGES = List.of(
            new PipelineStage("sequence_input", "Sequence Input", List.of()),
            new PipelineStage("alignment", "Alignment", List.of()),
            new PipelineStage("phylogeny", "Phylogeny", List.of(IndexKey.MU, IndexKey.RE, IndexKey.RI)),
            new PipelineStage("diversity", "Diversity", List.of(IndexKey.PI)),
            new PipelineStage("evolution", "Evolution", List.of(IndexKey.MU, IndexKey.DNDS)),
            new PipelineStage("recombination", "Recombination", List.of(IndexKey.RI)),
            new PipelineStage("temporal", "Temporal Analysis", List.of(IndexKey.RE)),
            new PipelineStage("scoring", "GVI Scoring", List.of()),
            new PipelineStage("report", "Final Report", List.of())
    );

    private final Map<String, Label> iconLabels = new LinkedHashMap<>();
    private final Map<String, Label> textLabels = new LinkedHashMap<>();

    public PipelineStagesView() {
        setSpacing(4);
        for (PipelineStage stage : STAGES) {
            Label icon = new Label();
            icon.getStyleClass().add("gvi-stage-icon");
            Label text = new Label(stage.label());
            text.getStyleClass().add("gvi-stage-label");
            iconLabels.put(stage.id(), icon);
            textLabels.put(stage.id(), text);

            Region spacer = new Region();
            HBox.setHgrow(spacer, Priority.ALWAYS);
            HBox row = new HBox(8, icon, text);
            row.setAlignment(Pos.CENTER_LEFT);
            getChildren().add(row);
        }

        Label footnote = new Label("Mutation Burden / Genetic Distance / Codon Adaptation / GC Deviation are "
                + "computed directly per sequence and feed straight into GVI Scoring -- they aren't shown as "
                + "separate stages here.");
        footnote.setWrapText(true);
        footnote.getStyleClass().add("gvi-placeholder-text");
        getChildren().add(footnote);

        reset();
    }

    public void reset() {
        for (PipelineStage stage : STAGES) apply(stage.id(), Status.PENDING, null);
    }

    /** @param requestedIndices the lowercase index keys the user actually selected for this run. */
    public void setRunning(Set<String> requestedIndices) {
        for (PipelineStage stage : STAGES) {
            if (stage.indexKeys().isEmpty() || stageRequested(stage, requestedIndices)) {
                apply(stage.id(), Status.RUNNING, null);
            } else {
                apply(stage.id(), Status.NOT_SELECTED, "not selected");
            }
        }
    }

    public void applyFailure() {
        for (PipelineStage stage : STAGES) apply(stage.id(), Status.FAILED, "run failed -- see the error dialog");
    }

    public void applyResult(PipelineResult result, Set<String> requestedIndices) {
        Map<IndexKey, GviComponent> componentsByKey = new EnumMap<>(IndexKey.class);
        List<IndexKey> excluded = result.datasetGvi() != null ? result.datasetGvi().excludedIndices() : List.of();
        if (result.datasetGvi() != null) {
            for (GviComponent c : result.datasetGvi().components()) componentsByKey.put(c.key(), c);
        }

        for (PipelineStage stage : STAGES) {
            if (!stage.indexKeys().isEmpty() && !stageRequested(stage, requestedIndices)) {
                apply(stage.id(), Status.NOT_SELECTED, "not selected this run");
                continue;
            }
            switch (stage.id()) {
                case "sequence_input", "alignment" -> apply(stage.id(), Status.COMPLETED, null);
                case "scoring" -> {
                    if (result.datasetGvi() != null) {
                        apply(stage.id(), Status.COMPLETED, null);
                    } else {
                        apply(stage.id(), Status.FAILED, "whole-file GVI could not be computed -- see Skipped / Warnings");
                    }
                }
                case "report" -> apply(stage.id(), Status.COMPLETED, null);
                default -> applyIndexBackedStage(stage, componentsByKey, excluded);
            }
        }
    }

    private void applyIndexBackedStage(PipelineStage stage, Map<IndexKey, GviComponent> componentsByKey, List<IndexKey> excluded) {
        boolean anyComputed = stage.indexKeys().stream().anyMatch(componentsByKey::containsKey);
        boolean anyExcluded = stage.indexKeys().stream().anyMatch(excluded::contains);
        if (anyComputed) {
            apply(stage.id(), Status.COMPLETED, null);
        } else if (anyExcluded) {
            apply(stage.id(), Status.WARNING, "excluded this run -- see Skipped / Warnings");
        } else {
            apply(stage.id(), Status.FAILED, "could not be computed -- see Skipped / Warnings");
        }
    }

    private boolean stageRequested(PipelineStage stage, Set<String> requestedIndices) {
        return stage.indexKeys().stream().anyMatch(k -> requestedIndices.contains(requestKey(k)));
    }

    private String requestKey(IndexKey key) {
        return switch (key) {
            case MU -> "mu";
            case RE -> "re";
            case PI -> "pi";
            case MB -> "mb";
            case DNDS -> "dnds";
            case GD -> "gd";
            case CAI -> "cai";
            case GC -> "gc";
            case RI -> "ri";
        };
    }

    private void apply(String stageId, Status status, String note) {
        Label icon = iconLabels.get(stageId);
        Label text = textLabels.get(stageId);
        if (icon == null || text == null) return;
        icon.getStyleClass().removeAll("stage-pending", "stage-running", "stage-completed", "stage-warning",
                "stage-failed", "stage-skipped");
        String glyph;
        String styleClass;
        switch (status) {
            case PENDING -> { glyph = "○"; styleClass = "stage-pending"; }
            case RUNNING -> { glyph = "●"; styleClass = "stage-running"; }
            case COMPLETED -> { glyph = "✓"; styleClass = "stage-completed"; }
            case WARNING -> { glyph = "⚠"; styleClass = "stage-warning"; }
            case FAILED -> { glyph = "✕"; styleClass = "stage-failed"; }
            default -> { glyph = "—"; styleClass = "stage-skipped"; }
        }
        icon.setText(glyph);
        icon.getStyleClass().add(styleClass);
        PipelineStage stage = STAGES.stream().filter(s -> s.id().equals(stageId)).findFirst().orElseThrow();
        text.setText(stage.label() + (note != null ? "  (" + note + ")" : ""));
    }
}
