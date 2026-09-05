package org.gvi.ui.cards;

import javafx.geometry.Pos;
import javafx.scene.layout.HBox;
import javafx.scene.layout.Priority;
import javafx.scene.layout.Region;
import javafx.scene.layout.VBox;
import javafx.scene.control.Label;
import org.gvi.composite.GviComponent;
import org.gvi.composite.IndexKey;
import org.gvi.ui.theme.MetricDescriptions;

import java.util.EnumMap;
import java.util.Locale;
import java.util.Map;

/**
 * One metric card: display name, the real computed value from this run (or an honest not-available /
 * not-run-yet state -- never fabricated), a short static description, and a status badge. Clicking a
 * card invokes a caller-supplied callback (used to open a detail dialog) rather than owning any
 * dialog/Alert logic itself, keeping this class a plain, easily-testable node builder.
 */
public final class MetricCard {

    private static final Map<IndexKey, String> DISPLAY_NAMES = new EnumMap<>(IndexKey.class);

    static {
        DISPLAY_NAMES.put(IndexKey.MU, "Evolutionary Rate");
        DISPLAY_NAMES.put(IndexKey.RE, "Effective Reproduction Number");
        DISPLAY_NAMES.put(IndexKey.PI, "Genetic Diversity");
        DISPLAY_NAMES.put(IndexKey.MB, "Mutation Burden");
        DISPLAY_NAMES.put(IndexKey.DNDS, "Selection Pressure (dN/dS)");
        DISPLAY_NAMES.put(IndexKey.GD, "Genetic Distance");
        DISPLAY_NAMES.put(IndexKey.CAI, "Codon Adaptation");
        DISPLAY_NAMES.put(IndexKey.GC, "GC Content Deviation");
        DISPLAY_NAMES.put(IndexKey.RI, "Recombination");
    }

    private MetricCard() {
    }

    public static String displayName(IndexKey key) {
        return DISPLAY_NAMES.getOrDefault(key, key.label());
    }

    /**
     * @param component real computed component for this run, or {@code null} if unavailable
     * @param excluded  true if this index was attempted this run but excluded for missing data
     * @param onClick   invoked when the card is clicked (e.g. to open a detail dialog); may be null
     */
    public static VBox build(IndexKey key, GviComponent component, boolean excluded, Runnable onClick) {
        Label title = new Label(displayName(key));
        title.getStyleClass().add("gvi-card-metric-title");

        Label status = new Label();
        status.getStyleClass().add("gvi-badge");

        Label value = new Label();
        value.getStyleClass().add("gvi-card-metric-value");

        if (component != null) {
            value.setText(key.label() + " = " + formatValue(component.rawValue()));
            status.setText("Calculated");
            status.getStyleClass().add("badge-pass");
        } else if (excluded) {
            value.setText("--");
            status.setText("Excluded (no data)");
            status.getStyleClass().add("badge-excluded");
        } else {
            value.setText("--");
            status.setText("Not run yet");
            status.getStyleClass().add("badge-default");
        }

        Label description = new Label(MetricDescriptions.describe(key));
        description.setWrapText(true);
        description.getStyleClass().add("gvi-card-metric-desc");

        Region spacer = new Region();
        HBox.setHgrow(spacer, Priority.ALWAYS);
        HBox header = new HBox(8, title, spacer, status);
        header.setAlignment(Pos.CENTER_LEFT);

        VBox card = new VBox(6, header, value, description);
        card.getStyleClass().addAll("gvi-card", "gvi-metric-card");
        if (onClick != null) {
            card.setOnMouseClicked(e -> onClick.run());
        }
        return card;
    }

    private static String formatValue(double v) {
        if (v != 0 && (Math.abs(v) < 0.001 || Math.abs(v) >= 100000)) {
            return String.format(Locale.ROOT, "%.4g", v);
        }
        return String.format(Locale.ROOT, "%.4f", v);
    }
}
