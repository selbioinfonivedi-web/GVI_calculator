package org.gvi.ui.nav;

import javafx.geometry.Pos;
import javafx.scene.control.Button;
import javafx.scene.control.Label;
import javafx.scene.layout.VBox;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Consumer;

/**
 * Left workflow navigation: grouped, clickable steps. Uses plain {@link Button}s with a manually
 * tracked "selected" style class rather than a {@link javafx.scene.control.ToggleGroup} -- a
 * ToggleGroup deselects its toggle when the already-selected one is clicked again, which would blank
 * the content pane; tracking selection by hand avoids that pitfall entirely.
 */
public final class NavRail extends VBox {

    private final Map<String, Button> buttons = new LinkedHashMap<>();
    private Button selectedButton;

    public NavRail(List<WorkflowStep> steps, Consumer<WorkflowStep> onSelect) {
        getStyleClass().add("gvi-nav-rail");
        setSpacing(2);

        NavGroup currentGroup = null;
        for (WorkflowStep step : steps) {
            if (currentGroup != step.group()) {
                currentGroup = step.group();
                Label groupLabel = new Label(currentGroup.label().toUpperCase());
                groupLabel.getStyleClass().add("gvi-nav-group-label");
                getChildren().add(groupLabel);
            }

            Button btn = new Button(step.label());
            btn.getStyleClass().add("gvi-nav-item");
            btn.setMaxWidth(Double.MAX_VALUE);
            btn.setAlignment(Pos.CENTER_LEFT);
            btn.setOnAction(e -> {
                select(step.id());
                onSelect.accept(step);
            });
            buttons.put(step.id(), btn);
            getChildren().add(btn);
        }
    }

    public void select(String stepId) {
        Button btn = buttons.get(stepId);
        if (btn == null) return;
        if (selectedButton != null) {
            selectedButton.getStyleClass().remove("gvi-nav-item-selected");
        }
        btn.getStyleClass().add("gvi-nav-item-selected");
        selectedButton = btn;
    }
}
