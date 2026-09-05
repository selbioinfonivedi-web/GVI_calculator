package org.gvi.ui;

import javafx.geometry.Insets;
import javafx.scene.control.Label;
import javafx.scene.control.ScrollPane;
import javafx.scene.layout.VBox;
import javafx.scene.paint.Color;
import javafx.scene.text.Font;
import javafx.scene.text.FontWeight;
import org.gvi.cli.PipelineResult;

import java.util.List;

/**
 * Why indices were dropped from the score, grouped by cause and paired with the remedy.
 * <p>
 * This information previously existed only as an undifferentiated text dump in a tab named
 * "Skipped / Warnings", where an exclusion that changes the meaning of the headline number sat
 * beside routine notices about defaulted options. Yet the exclusions are the most actionable
 * output the tool produces: each one names a specific, fixable problem with the input, and several
 * are fixable in minutes once the cause is legible.
 * <p>
 * Grouping is by cause rather than by index because the causes are what the user acts on -- five
 * indices failing for one bad alignment is one problem to fix, not five.
 */
public final class QualityGatingView extends ScrollPane {

    private static final Color INK = Color.web("#1a1a1a");
    private static final Color MUTED = Color.web("#6b7280");
    private static final Color WARN = Color.web("#b2560d");
    private static final Color OK = Color.web("#2f6b4f");

    private final VBox body = new VBox(12);

    public QualityGatingView() {
        body.setPadding(new Insets(14));
        setContent(body);
        setFitToWidth(true);
        clear();
    }

    public void clear() {
        body.getChildren().setAll(muted("No result yet."));
    }

    public void show(PipelineResult result) {
        body.getChildren().clear();

        List<String> exclusions = result.skipped().stream()
                .filter(s -> s.contains("excluded from composite")).toList();

        if (exclusions.isEmpty()) {
            body.getChildren().add(heading("Every computed index contributed to the score.", OK));
            body.getChildren().add(muted("Nothing was dropped by quality gating. Indices absent from the report were "
                    + "not computed at all -- usually because an optional input was not supplied."));
            addRemainingNotices(result);
            return;
        }

        body.getChildren().add(heading(exclusions.size() + " index/indices excluded from the composite GVI", WARN));
        body.getChildren().add(muted("These were computed and are reported in full -- they are only removed from the "
                + "weighted score, and the remaining weights were renormalized. Each entry names what to fix."));

        // Grouping and wording come from the shared explainer, so the desktop app, the CLI report
        // and the web front end all describe the same failure identically.
        for (var group : org.gvi.cli.ExclusionExplainer.explain(result.skipped())) {
            VBox card = new VBox(4);
            card.setPadding(new Insets(10, 12, 10, 12));
            card.setStyle("-fx-background-color: #fff; -fx-border-color: #e3e5e8; -fx-border-radius: 3; -fx-background-radius: 3;");
            Label cause = new Label(group.cause());
            cause.setFont(Font.font("System", FontWeight.BOLD, 12));
            cause.setTextFill(INK);
            cause.setWrapText(true);
            Label affects = new Label("Affects: " + String.join(", ", group.affectedIndices()));
            affects.setFont(Font.font(11));
            affects.setTextFill(MUTED);
            Label fix = new Label(group.remedy());
            fix.setFont(Font.font(11));
            fix.setTextFill(WARN);
            fix.setWrapText(true);
            card.getChildren().addAll(cause, affects, fix);
            body.getChildren().add(card);
        }
        addRemainingNotices(result);
    }

    /** Everything that is not an exclusion still matters, but belongs below it. */
    private void addRemainingNotices(PipelineResult result) {
        List<String> notices = result.warnings();
        if (notices.isEmpty()) return;
        body.getChildren().add(heading("Other notices (" + notices.size() + ")", INK));
        for (String w : notices) {
            Label l = muted("• " + w);
            l.setWrapText(true);
            body.getChildren().add(l);
        }
    }

    private static Label heading(String text, Color color) {
        Label l = new Label(text);
        l.setFont(Font.font("System", FontWeight.BOLD, 13));
        l.setTextFill(color);
        l.setWrapText(true);
        return l;
    }

    private static Label muted(String text) {
        Label l = new Label(text);
        l.setFont(Font.font(11));
        l.setTextFill(MUTED);
        l.setWrapText(true);
        return l;
    }
}
