package org.gvi.ui;

/**
 * Plain (non-Application) main class. Launching {@link javafx.application.Application}
 * subclasses directly as a jar's Main-Class can fail with "JavaFX runtime
 * components are missing" on some classpath (non-modular) setups; going
 * through an indirection class that just calls {@code Application.launch}
 * avoids that class-loading pitfall.
 */
public final class GviUiLauncher {
    public static void main(String[] args) {
        GviUiApp.main(args);
    }
}
