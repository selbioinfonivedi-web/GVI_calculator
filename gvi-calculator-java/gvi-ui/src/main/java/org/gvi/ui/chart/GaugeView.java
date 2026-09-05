package org.gvi.ui.chart;

import javafx.geometry.Pos;
import javafx.scene.control.Label;
import javafx.scene.layout.StackPane;
import javafx.scene.layout.VBox;
import javafx.scene.shape.Arc;
import javafx.scene.shape.ArcType;
import javafx.scene.shape.StrokeLineCap;

/**
 * A semicircular progress gauge for the whole-file composite GVI score. Deliberately a single neutral
 * color (the app's primary blue) with no risk-tier color bands: the real backend defines no severity
 * scale for the composite score anywhere ({@code PipelineResult}/{@code GviResult}/
 * {@code CompositeGviEngine} carry no "low/moderate/high" concept), so inventing a color-coded scale
 * here would be presenting an interpretation the tool itself doesn't make -- exactly what this
 * redesign's own "scientific transparency" requirement warns against.
 */
public final class GaugeView extends StackPane {

    private static final double RADIUS = 90;
    private static final double STROKE_WIDTH = 16;

    private final Arc valueArc;
    private final Label valueLabel;
    private final Label captionLabel;

    public GaugeView() {
        double size = RADIUS * 2 + STROKE_WIDTH;
        setPrefSize(size, RADIUS + STROKE_WIDTH);
        setMinSize(size, RADIUS + STROKE_WIDTH);
        setMaxSize(size, RADIUS + STROKE_WIDTH);

        double cx = size / 2;
        double cy = RADIUS + STROKE_WIDTH / 2;

        Arc track = new Arc(cx, cy, RADIUS, RADIUS, 180, -180);
        track.setType(ArcType.OPEN);
        track.setFill(null);
        track.getStyleClass().add("gvi-gauge-track");
        track.setStrokeWidth(STROKE_WIDTH);
        track.setStrokeLineCap(StrokeLineCap.ROUND);

        valueArc = new Arc(cx, cy, RADIUS, RADIUS, 180, 0);
        valueArc.setType(ArcType.OPEN);
        valueArc.setFill(null);
        valueArc.getStyleClass().add("gvi-gauge-value");
        valueArc.setStrokeWidth(STROKE_WIDTH);
        valueArc.setStrokeLineCap(StrokeLineCap.ROUND);

        valueLabel = new Label("--");
        valueLabel.getStyleClass().add("gvi-gauge-value-label");
        captionLabel = new Label("No analysis run yet");
        captionLabel.getStyleClass().add("gvi-gauge-caption");
        VBox labels = new VBox(2, valueLabel, captionLabel);
        labels.setAlignment(Pos.BOTTOM_CENTER);
        labels.setMaxHeight(RADIUS);
        StackPane.setAlignment(labels, Pos.BOTTOM_CENTER);

        getChildren().addAll(track, valueArc, labels);
    }

    /** @param fraction the GVI score on its native 0..1 scale, or {@code null} for the not-yet-run state. */
    public void setValue(Double fraction) {
        if (fraction == null) {
            valueArc.setLength(0);
            valueLabel.setText("--");
            captionLabel.setText("No analysis run yet");
            return;
        }
        double clamped = Math.max(0, Math.min(1, fraction));
        valueArc.setLength(-180 * clamped);
        valueLabel.setText(String.format("%.1f%%", clamped * 100));
        captionLabel.setText("Genomic Vulnerability");
    }
}
