package org.gvi.ui.theme;

import javafx.scene.text.Font;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.InputStream;

/**
 * Loads the bundled IBM Plex Sans/Mono weights so the CSS in {@code theme/app.css} can reference
 * them by family name. Must run before {@code scene.getStylesheets().add(...)} -- JavaFX CSS resolves
 * {@code -fx-font-family} against whatever font families are already registered at stylesheet-apply
 * time, so loading fonts afterward risks a silent fallback to the system default.
 *
 * <p>TrueType only supports Regular/Bold/Italic/BoldItalic as style variants of one family name --
 * IBM Plex's Medium and SemiBold weights are shipped as their own separate family names rather than
 * as weight variants of "IBM Plex Sans" (confirmed by inspecting each file's embedded name table).
 */
public final class Fonts {

    private static final Logger log = LoggerFactory.getLogger(Fonts.class);

    private static final String[] FILES = {
            "/fonts/IBMPlexSans-Regular.ttf",
            "/fonts/IBMPlexSans-Medium.ttf",
            "/fonts/IBMPlexSans-SemiBold.ttf",
            "/fonts/IBMPlexSans-Bold.ttf",
            "/fonts/IBMPlexMono-Regular.ttf",
    };

    private Fonts() {
    }

    public static void init() {
        for (String path : FILES) {
            try (InputStream in = Fonts.class.getResourceAsStream(path)) {
                if (in == null) {
                    log.warn("Font resource not found: {} -- falling back to system default for this weight", path);
                    continue;
                }
                Font font = Font.loadFont(in, 12);
                if (font == null) {
                    log.warn("Font failed to load: {} -- falling back to system default for this weight", path);
                } else {
                    log.info("Loaded font: {} -> family \"{}\"", path, font.getFamily());
                }
            } catch (Exception e) {
                log.warn("Error loading font {}: {}", path, e.getMessage());
            }
        }
    }
}
