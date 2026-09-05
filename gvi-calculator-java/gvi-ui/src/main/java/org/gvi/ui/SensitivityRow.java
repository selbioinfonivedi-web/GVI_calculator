package org.gvi.ui;

import org.gvi.cli.PipelineResult;
import org.gvi.composite.SensitivityResult;

import java.util.Comparator;
import java.util.List;

/**
 * One row of the sensitivity-analysis table: how much the whole-file GVI
 * would move if a single index's weight were perturbed +/-20% (Section
 * 5.9/8), sorted most-sensitive-first so the "tornado chart" ordering is
 * visible without the user having to sort it themselves. Plain JavaBean
 * getters so JavaFX's {@code PropertyValueFactory} can bind by name.
 */
public final class SensitivityRow {

    private final SensitivityResult result;

    SensitivityRow(SensitivityResult result) {
        this.result = result;
    }

    public static List<SensitivityRow> fromResult(PipelineResult result) {
        return result.sensitivity().stream()
                .sorted(Comparator.comparingDouble(SensitivityResult::spread).reversed())
                .map(SensitivityRow::new)
                .toList();
    }

    public String getIndex() {
        return result.key().label();
    }

    public String getBaseGvi() {
        return String.format("%.4f", result.baseGvi());
    }

    public String getLowWeightGvi() {
        return String.format("%.4f", result.gviAtLowWeight());
    }

    public String getHighWeightGvi() {
        return String.format("%.4f", result.gviAtHighWeight());
    }

    public String getSpread() {
        return String.format("%.4f", result.spread());
    }
}
