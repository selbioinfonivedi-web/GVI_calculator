package org.gvi.ui.chart;

import javafx.geometry.Pos;
import javafx.scene.Node;
import javafx.scene.control.Label;
import javafx.scene.layout.HBox;
import javafx.scene.layout.Pane;
import javafx.scene.layout.Region;
import org.gvi.composite.SensitivityResult;
import org.gvi.ui.cards.MetricCard;

import java.util.Comparator;
import java.util.List;
import java.util.Locale;
import java.util.stream.Stream;

/**
 * Hand-rolled tornado chart for the sensitivity analysis: one row per index, a shared horizontal
 * value scale, and two bar segments running from the (shared, real) base GVI out to the GVI at -20%
 * and +20% weight -- rows sorted most-sensitive first, the standard tornado-chart shape. Same
 * hand-rolled-{@link Region} approach as {@link ContributionBars}, for the same reasons (no
 * BarChart/Canvas precedent in this codebase).
 */
public final class TornadoBars {

    private static final double TRACK_WIDTH = 300;
    private static final double BAR_HEIGHT = 12;

    private TornadoBars() {
    }

    public static Node build(List<SensitivityResult> sensitivity) {
        if (sensitivity.isEmpty()) {
            Label empty = new Label("Run the pipeline to see how much each index's weight moves the whole-file GVI.");
            empty.getStyleClass().add("gvi-placeholder-text");
            return empty;
        }
        List<SensitivityResult> sorted = sensitivity.stream()
                .sorted(Comparator.comparingDouble((SensitivityResult s) -> -Math.abs(s.spread())))
                .toList();

        double base = sorted.get(0).baseGvi();
        double lo = Stream.concat(sorted.stream().map(SensitivityResult::gviAtLowWeight),
                        sorted.stream().map(SensitivityResult::gviAtHighWeight))
                .mapToDouble(Double::doubleValue).min().orElse(base);
        double hi = Stream.concat(sorted.stream().map(SensitivityResult::gviAtLowWeight),
                        sorted.stream().map(SensitivityResult::gviAtHighWeight))
                .mapToDouble(Double::doubleValue).max().orElse(base);
        lo = Math.min(lo, base);
        hi = Math.max(hi, base);
        double pad = (hi - lo) > 0 ? (hi - lo) * 0.1 : 0.01;
        lo -= pad;
        hi += pad;

        javafx.scene.layout.VBox rows = new javafx.scene.layout.VBox(10);
        for (SensitivityResult s : sorted) {
            rows.getChildren().add(buildRow(s, lo, hi));
        }
        return rows;
    }

    private static Node buildRow(SensitivityResult s, double lo, double hi) {
        double baseX = xFor(s.baseGvi(), lo, hi);
        double lowX = xFor(s.gviAtLowWeight(), lo, hi);
        double highX = xFor(s.gviAtHighWeight(), lo, hi);

        Pane track = new Pane();
        track.getStyleClass().add("gvi-tornado-track");
        track.setPrefSize(TRACK_WIDTH, BAR_HEIGHT + 4);
        track.setMinSize(TRACK_WIDTH, BAR_HEIGHT + 4);
        track.setMaxSize(TRACK_WIDTH, BAR_HEIGHT + 4);

        Region lowBar = new Region();
        lowBar.getStyleClass().add("gvi-tornado-bar-low");
        lowBar.setLayoutX(Math.min(baseX, lowX));
        lowBar.setLayoutY(2);
        lowBar.setPrefSize(Math.max(0, Math.abs(baseX - lowX)), BAR_HEIGHT);

        Region highBar = new Region();
        highBar.getStyleClass().add("gvi-tornado-bar-high");
        highBar.setLayoutX(Math.min(baseX, highX));
        highBar.setLayoutY(2);
        highBar.setPrefSize(Math.max(0, Math.abs(highX - baseX)), BAR_HEIGHT);

        Region baseMarker = new Region();
        baseMarker.getStyleClass().add("gvi-tornado-base-marker");
        baseMarker.setLayoutX(baseX - 1);
        baseMarker.setLayoutY(0);
        baseMarker.setPrefSize(2, BAR_HEIGHT + 4);

        track.getChildren().addAll(lowBar, highBar, baseMarker);

        Label name = new Label(MetricCard.displayName(s.key()));
        name.getStyleClass().add("gvi-bar-label");
        name.setPrefWidth(190);

        Label spreadLabel = new Label(String.format(Locale.ROOT, "spread %.4f", Math.abs(s.spread())));
        spreadLabel.getStyleClass().add("gvi-bar-value");

        HBox row = new HBox(10, name, track, spreadLabel);
        row.setAlignment(Pos.CENTER_LEFT);
        return row;
    }

    private static double xFor(double value, double lo, double hi) {
        if (hi <= lo) return TRACK_WIDTH / 2;
        return ((value - lo) / (hi - lo)) * TRACK_WIDTH;
    }
}
