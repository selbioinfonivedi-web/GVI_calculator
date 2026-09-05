package org.gvi.ui.chart;

import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.canvas.Canvas;
import javafx.scene.canvas.GraphicsContext;
import javafx.scene.control.Label;
import javafx.scene.control.Tooltip;
import javafx.scene.layout.HBox;
import javafx.scene.layout.Priority;
import javafx.scene.layout.VBox;
import javafx.scene.paint.Color;
import javafx.scene.text.Font;
import javafx.scene.text.FontWeight;
import org.gvi.composite.GviComponent;
import org.gvi.composite.GviResult;
import org.gvi.composite.IndexKey;

import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

/**
 * The headline result: the score, what it was built from, and whether it can be compared to
 * anything else.
 * <p>
 * A bare GVI number is not interpretable on its own, because the composite renormalizes whatever
 * survived quality gating to sum to 1 -- so a score built from two indices and a score built from
 * nine both land in [0,1] and look identical. Enterotoxaemia's 0.053 from 6% of the weighting
 * scheme sitting beside haemorrhagic septicaemia's 0.680 reads as "much lower risk" when it
 * actually means the data was too broken to score. This strip makes that difference the first thing
 * on screen rather than something buried in a warnings tab.
 * <p>
 * Contributing indices are drawn as filled segments proportional to their actual contribution;
 * excluded ones follow as hatched blocks sized by the weight they would have carried, so the gap in
 * the picture is visible rather than merely absent.
 */
public final class CompositionStrip extends VBox {

    private static final Color SCORED = Color.web("#2166ac");
    private static final Color TRACK = Color.web("#eef0f3");
    private static final Color RULE = Color.web("#c9ccd1");
    private static final Color MUTED = Color.web("#6b7280");
    private static final Color WARN = Color.web("#b2560d");

    private final Canvas canvas = new Canvas(760, 34);
    private final Label headline = new Label();
    private final Label coverage = new Label();
    private final Label caveat = new Label();

    public CompositionStrip() {
        setSpacing(6);
        setPadding(new Insets(10, 12, 10, 12));

        headline.setFont(Font.font("System", FontWeight.BOLD, 15));
        coverage.setFont(Font.font(11));
        coverage.setTextFill(MUTED);
        caveat.setFont(Font.font("System", FontWeight.BOLD, 11));
        caveat.setTextFill(WARN);
        caveat.setWrapText(true);
        caveat.setManaged(false);
        caveat.setVisible(false);

        HBox top = new HBox(12, headline, coverage);
        top.setAlignment(Pos.BASELINE_LEFT);
        HBox.setHgrow(coverage, Priority.ALWAYS);

        getChildren().addAll(top, caveat, canvas);
        clear();
    }

    public void clear() {
        headline.setText("No result yet");
        coverage.setText("");
        caveat.setVisible(false);
        caveat.setManaged(false);
        GraphicsContext g = canvas.getGraphicsContext2D();
        g.clearRect(0, 0, canvas.getWidth(), canvas.getHeight());
    }

    /** @param sequenceCount only for the headline; the strip itself is about index composition. */
    public void show(GviResult gvi, int sequenceCount) {
        if (gvi == null) {
            headline.setText("Whole-file GVI could not be computed");
            coverage.setText("See the Quality Gating tab for why.");
            caveat.setVisible(false);
            caveat.setManaged(false);
            canvas.getGraphicsContext2D().clearRect(0, 0, canvas.getWidth(), canvas.getHeight());
            return;
        }

        headline.setText(String.format(Locale.ROOT, "GVI = %.4f", gvi.gvi()));
        coverage.setText(String.format(Locale.ROOT, "%s  ·  %d sequence(s)", gvi.coverageSummary(), sequenceCount));

        boolean comparable = gvi.comparable();
        caveat.setVisible(!comparable);
        caveat.setManaged(!comparable);
        if (!comparable) {
            caveat.setText("Not comparable with other datasets: too little of the weighting scheme had usable data. "
                    + "This score summarises only the indices that survived gating, and must not be ranked against a "
                    + "fully-populated GVI.");
        }
        draw(gvi);
    }

    private void draw(GviResult gvi) {
        GraphicsContext g = canvas.getGraphicsContext2D();
        double w = canvas.getWidth(), h = canvas.getHeight();
        g.clearRect(0, 0, w, h);

        List<GviComponent> comps = new ArrayList<>(gvi.components());
        comps.sort((a, b) -> Double.compare(b.contribution(), a.contribution()));
        List<IndexKey> excluded = gvi.excludedIndices() == null ? List.of() : gvi.excludedIndices();

        // The strip spans the whole weighting scheme, not just the part that scored: contributing
        // segments are sized by effective weight, excluded ones by the weight they would have had.
        double totalWeight = 0;
        for (GviComponent c : comps) totalWeight += c.effectiveWeight();
        double excludedShare = Math.max(0, 1.0 - gvi.effectiveWeightSum());
        double perExcluded = excluded.isEmpty() ? 0 : excludedShare / excluded.size();
        double denom = totalWeight * gvi.effectiveWeightSum() + excludedShare;
        if (denom <= 0) denom = 1;

        g.setFill(TRACK);
        g.fillRect(0, 0, w, h);

        double x = 0;
        int i = 0;
        for (GviComponent c : comps) {
            double share = (c.effectiveWeight() * gvi.effectiveWeightSum()) / denom;
            double segW = share * w;
            // Opacity ramp keeps adjacent segments separable without needing nine distinct hues.
            g.setFill(SCORED.deriveColor(0, 1, 1, 0.45 + 0.55 * (1.0 - (double) i / Math.max(comps.size(), 1))));
            g.fillRect(x, 0, segW, h);
            g.setStroke(Color.WHITE);
            g.setLineWidth(0.8);
            g.strokeRect(x, 0, segW, h);
            if (segW > 34) {
                g.setFill(Color.WHITE);
                g.setFont(Font.font("System", FontWeight.BOLD, 9));
                g.fillText(c.key().label(), x + 4, h / 2 + 3);
            }
            x += segW;
            i++;
        }

        for (IndexKey k : excluded) {
            double segW = (perExcluded / denom) * w;
            g.setFill(Color.WHITE);
            g.fillRect(x, 0, segW, h);
            g.setStroke(MUTED);
            g.setLineWidth(0.7);
            for (double d = -h; d < segW; d += 5) {
                double x0 = Math.max(x, x + d), y0 = x + d < x ? h - (x - (x + d)) : h;
                double x1 = Math.min(x + segW, x + d + h), y1 = h - (x1 - (x + d));
                if (x1 > x0) g.strokeLine(x0, y0, x1, Math.max(0, y1));
            }
            g.setStroke(RULE);
            g.setLineWidth(0.8);
            g.strokeRect(x, 0, segW, h);
            if (segW > 30) {
                g.setFill(MUTED);
                g.setFont(Font.font("System", FontWeight.NORMAL, 9));
                g.fillText(k.label(), x + 4, h / 2 + 3);
            }
            x += segW;
        }

        StringBuilder tip = new StringBuilder("Contributing:\n");
        for (GviComponent c : comps) {
            tip.append(String.format(Locale.ROOT, "  %-14s weight %.3f -> %.4f%n",
                    c.key().label(), c.effectiveWeight(), c.contribution()));
        }
        if (!excluded.isEmpty()) {
            tip.append("\nExcluded (reported, not scored):\n");
            for (IndexKey k : excluded) tip.append("  ").append(k.label()).append('\n');
        }
        Tooltip.install(canvas, new Tooltip(tip.toString().trim()));
    }
}
