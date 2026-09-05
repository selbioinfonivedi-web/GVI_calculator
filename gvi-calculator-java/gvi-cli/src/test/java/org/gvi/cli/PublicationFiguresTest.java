package org.gvi.cli;

import org.gvi.composite.GviComponent;
import org.gvi.composite.GviResult;
import org.gvi.composite.IndexKey;
import org.gvi.composite.SensitivityResult;
import org.gvi.core.spi.SimpleIndexResult;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.EnumMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static org.assertj.core.api.Assertions.assertThat;

/**
 * SVG clips silently: anything drawn outside the viewBox simply disappears, with no error and no
 * visual hint that a figure is incomplete. A real case -- the excluded-index block on FMD was drawn
 * at y=230 on a 211-tall canvas and vanished, while every value in the figure was still correct.
 * These tests assert that nothing is drawn outside the declared canvas.
 */
class PublicationFiguresTest {

    private static final Pattern HEIGHT = Pattern.compile("height='([0-9.]+)'");
    private static final Pattern TEXT_Y = Pattern.compile("<text[^>]*y='([0-9.]+)'");
    private static final Pattern RECT = Pattern.compile("<rect[^>]*y='([0-9.]+)'[^>]*height='([0-9.]+)'");

    private static PipelineResult resultWith(int contributing, int excludedCount, boolean comparable) {
        List<GviComponent> comps = new ArrayList<>();
        IndexKey[] keys = IndexKey.values();
        double w = 1.0 / Math.max(contributing, 1);
        for (int i = 0; i < contributing; i++) {
            comps.add(new GviComponent(keys[i], 0.5, 0.5, w, w, w * 0.5));
        }
        List<IndexKey> excluded = new ArrayList<>();
        for (int i = contributing; i < contributing + excludedCount && i < keys.length; i++) excluded.add(keys[i]);

        GviResult gvi = new GviResult(0.5, comps, excluded, comparable ? 0.9 : 0.2, List.of());
        Map<IndexKey, org.gvi.core.spi.IndexResult> dataset = new EnumMap<>(IndexKey.class);
        for (IndexKey k : excluded) dataset.put(k, new SimpleIndexResult(k.label(), 1.23, "cat", List.of()));

        List<SensitivityResult> sens = new ArrayList<>();
        for (GviComponent c : comps) sens.add(new SensitivityResult(c.key(), 0.5, 0.48, 0.52));

        return new PipelineResult(Map.of(), dataset, gvi, Map.of(), Map.of(), Map.of(), Map.of(), Map.of(),
                Map.of(), sens, List.of(), List.of(), null, List.of());
    }

    /** Largest y coordinate any element actually paints to. */
    private static double lowestDrawnY(String svg) {
        double max = 0;
        Matcher t = TEXT_Y.matcher(svg);
        while (t.find()) max = Math.max(max, Double.parseDouble(t.group(1)));
        Matcher r = RECT.matcher(svg);
        while (r.find()) max = Math.max(max, Double.parseDouble(r.group(1)) + Double.parseDouble(r.group(2)));
        return max;
    }

    private static double canvasHeight(String svg) {
        Matcher m = HEIGHT.matcher(svg);
        assertThat(m.find()).isTrue();
        return Double.parseDouble(m.group(1));
    }

    private static void assertNothingClipped(String svg, String what) {
        assertThat(lowestDrawnY(svg))
                .as("%s draws below its own canvas, so content is silently clipped", what)
                .isLessThanOrEqualTo(canvasHeight(svg));
    }

    @Test
    void compositionFigureFitsItsCanvasWithExcludedIndices() {
        // The exact shape that broke: 8 contributing, 1 excluded.
        assertNothingClipped(PublicationFigures.compositionFigure(resultWith(8, 1, true), "fmd"), "composition");
    }

    @Test
    void compositionFigureFitsAcrossTheFullRangeOfShapes() {
        for (int contributing = 1; contributing <= 9; contributing++) {
            for (int excluded = 0; excluded + contributing <= 9; excluded++) {
                for (boolean comparable : new boolean[]{true, false}) {
                    String svg = PublicationFigures.compositionFigure(
                            resultWith(contributing, excluded, comparable), "d");
                    assertNothingClipped(svg, contributing + " scored / " + excluded + " excluded");
                }
            }
        }
    }

    @Test
    void indexPanelAndSensitivityFitTheirCanvases() {
        PipelineResult r = resultWith(6, 3, true);
        assertNothingClipped(PublicationFigures.indexPanelFigure(r, "d"), "index panel");
        assertNothingClipped(PublicationFigures.sensitivityFigure(r, "d"), "sensitivity");
    }

    @Test
    void aNonComparableScoreCarriesTheCaveatOnTheFigureItself() {
        String svg = PublicationFigures.compositionFigure(resultWith(2, 7, false), "entero");
        assertThat(svg).contains("Not comparable across datasets");
    }

    @Test
    void aComparableScoreDoesNotCarryTheCaveat() {
        assertThat(PublicationFigures.compositionFigure(resultWith(8, 1, true), "fmd"))
                .doesNotContain("Not comparable");
    }

    @Test
    void excludedIndicesAreHatchedNotOmitted() {
        String svg = PublicationFigures.compositionFigure(resultWith(6, 3, true), "d");
        assertThat(svg).contains("Excluded from the score");
        // Each excluded index gets its own hatched swatch, drawn as diagonal strokes.
        assertThat(svg.split("stroke-width='0.70'", -1).length - 1)
                .as("diagonal hatch strokes across three excluded swatches")
                .isGreaterThanOrEqualTo(3);
    }

    /**
     * Hatching must be plain geometry, not an SVG pattern or a clipPath. Both are ignored by
     * common rasterizers (ImageMagick/rsvg): a pattern fill renders as flat black, turning
     * "excluded" into a solid blob, and an ignored clipPath lets the diagonals spill across
     * neighbouring labels. Both were observed on real output before this was changed, so the
     * absence of those constructs is the actual guarantee worth pinning.
     */
    @Test
    void hatchingUsesNoPatternOrClipPath() {
        PipelineResult r = resultWith(6, 3, true);
        for (String svg : List.of(PublicationFigures.compositionFigure(r, "d"),
                PublicationFigures.indexPanelFigure(r, "d"))) {
            assertThat(svg).doesNotContain("<pattern").doesNotContain("url(#hatch)").doesNotContain("clip-path");
        }
    }

    @Test
    void outputIsWellFormedXml() throws Exception {
        var f = javax.xml.parsers.DocumentBuilderFactory.newInstance();
        f.setNamespaceAware(true);
        PipelineResult r = resultWith(5, 4, true);
        for (String svg : List.of(PublicationFigures.compositionFigure(r, "d & <test>"),
                PublicationFigures.indexPanelFigure(r, "d"), PublicationFigures.sensitivityFigure(r, "d"))) {
            f.newDocumentBuilder().parse(new org.xml.sax.InputSource(new java.io.StringReader(svg)));
        }
    }
}
