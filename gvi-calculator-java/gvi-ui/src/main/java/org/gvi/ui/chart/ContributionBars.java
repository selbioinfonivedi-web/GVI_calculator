package org.gvi.ui.chart;

import javafx.geometry.Pos;
import javafx.scene.Node;
import javafx.scene.control.Label;
import javafx.scene.layout.HBox;
import javafx.scene.layout.Region;
import javafx.scene.layout.StackPane;
import javafx.scene.layout.VBox;
import org.gvi.composite.GviComponent;
import org.gvi.ui.cards.MetricCard;

import java.util.Comparator;
import java.util.List;
import java.util.Locale;

/**
 * Hand-rolled CSS-styled {@link Region} bars for the composite GVI's per-index contribution
 * decomposition -- no {@code javafx.scene.chart.BarChart} precedent exists anywhere in this codebase,
 * and a BarChart's default axis/legend/animation stack would need more stripping-down than building
 * these simple single-series bars from scratch. Bar width is unit-tested indirectly via
 * {@link #barWidthPx} being a pure function, isolated from the JavaFX Application thread.
 */
public final class ContributionBars {

    private static final double TRACK_WIDTH = 260;

    private ContributionBars() {
    }

    public static Node build(List<GviComponent> components) {
        if (components.isEmpty()) {
            Label empty = new Label("Run the pipeline to see the per-index contribution breakdown.");
            empty.getStyleClass().add("gvi-placeholder-text");
            return empty;
        }
        List<GviComponent> sorted = components.stream()
                .sorted(Comparator.comparingDouble((GviComponent c) -> -Math.abs(c.contribution())))
                .toList();
        double max = sorted.stream().mapToDouble(c -> Math.abs(c.contribution())).max().orElse(0);
        if (max <= 0) max = 1;

        VBox rows = new VBox(10);
        for (GviComponent c : sorted) {
            rows.getChildren().add(buildRow(c, max));
        }
        return rows;
    }

    private static Node buildRow(GviComponent c, double max) {
        Label name = new Label(MetricCard.displayName(c.key()));
        name.getStyleClass().add("gvi-bar-label");
        name.setPrefWidth(190);

        Region track = new Region();
        track.getStyleClass().add("gvi-bar-track");
        track.setPrefSize(TRACK_WIDTH, 14);
        track.setMinSize(TRACK_WIDTH, 14);
        track.setMaxSize(TRACK_WIDTH, 14);

        double w = barWidthPx(c.contribution(), max, TRACK_WIDTH);
        Region fill = new Region();
        fill.getStyleClass().add("gvi-bar-fill");
        fill.setPrefSize(w, 14);
        fill.setMinSize(w, 14);
        fill.setMaxSize(w, 14);

        StackPane barStack = new StackPane(track, fill);
        StackPane.setAlignment(fill, Pos.CENTER_LEFT);

        Label value = new Label(String.format(Locale.ROOT, "%+.4f", c.contribution()));
        value.getStyleClass().add("gvi-bar-value");

        HBox row = new HBox(10, name, barStack, value);
        row.setAlignment(Pos.CENTER_LEFT);
        return row;
    }

    /** Pure math, isolated from JavaFX nodes so it's directly unit-testable. */
    public static double barWidthPx(double contribution, double maxAbsContribution, double trackWidth) {
        if (maxAbsContribution <= 0) return 2;
        return Math.max(2, (Math.abs(contribution) / maxAbsContribution) * trackWidth);
    }
}
