package org.gvi.ui.nav;

import javafx.scene.Node;
import javafx.scene.layout.BorderPane;
import javafx.scene.layout.StackPane;

/**
 * Application shell: header on top, nav rail on the left, a single content area in the center that
 * pages get swapped into, status bar on the bottom. Pages are built once by the caller and swapped in
 * via {@link #showPage}, never rebuilt on nav clicks -- rebuilding would lose in-progress form state
 * (e.g. file paths already typed into the Sequence Data page).
 */
public final class AppShell extends BorderPane {

    private final StackPane content = new StackPane();

    public AppShell(Node header, Node navRail, Node statusBar) {
        setTop(header);
        setLeft(navRail);
        content.getStyleClass().add("gvi-content-area");
        setCenter(content);
        setBottom(statusBar);
    }

    public void showPage(Node page) {
        content.getChildren().setAll(page);
    }
}
