package org.gvi.cli;

import org.gvi.composite.GviComponent;
import org.gvi.composite.GviResult;
import org.gvi.composite.IndexKey;
import org.gvi.composite.SensitivityResult;
import org.gvi.core.spi.IndexResult;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

/**
 * Writes publication-quality figures as standalone SVG.
 * <p>
 * SVG rather than PNG because a figure that will be printed must be resolution-independent -- a
 * raster chart sized for a screen reproduces badly at 300 dpi in a journal column, and reviewers
 * routinely ask for vector. SVG is also hand-editable, so a coauthor can adjust a label in
 * Illustrator or Inkscape without regenerating anything. Written directly rather than through a
 * charting library: the output is a few hundred elements, and a dependency that emits its own
 * styling would be harder to make print-correct than the markup itself.
 * <p>
 * Design constraints applied throughout:
 * <ul>
 *   <li><b>Single-column width.</b> 89 mm at 3.78 px/mm, the standard narrow figure width for most
 *       journals, so the figure needs no rescaling when placed.</li>
 *   <li><b>Colour is never the only encoding.</b> Excluded indices are hatched as well as greyed,
 *       so the distinction survives greyscale printing and colour-vision deficiency.</li>
 *   <li><b>Type at print sizes.</b> Labels at 7-9 pt equivalent, which stay legible at final size
 *       rather than being sized for a screen and shrinking to nothing.</li>
 *   <li><b>The exclusion is the finding.</b> These figures show which indices were scored and which
 *       were not, because a GVI built from six indices is not the same quantity as one built from
 *       nine, and a figure that hides that invites a false comparison.</li>
 * </ul>
 */
public final class PublicationFigures {

    /** 89 mm single-column width at 3.78 px/mm. */
    private static final double COLUMN_W = 336.0;
    private static final String FONT = "Helvetica, Arial, 'Liberation Sans', sans-serif";

    // Colour-vision-safe: blue / orange / teal / grey, distinguishable under deuteranopia and in greyscale.
    private static final String INK = "#1a1a1a";
    private static final String MUTED = "#6b7280";
    private static final String RULE = "#d4d4d8";
    private static final String SCORED = "#2166ac";
    private static final String ACCENT = "#b2560d";
    private static final String EXCLUDED_FILL = "#c9ccd1";

    private PublicationFigures() {
    }

    /** Writes every figure for one result. Returns the files created. */
    public static List<Path> writeAll(PipelineResult result, Path dir, String datasetLabel) throws IOException {
        Files.createDirectories(dir);
        List<Path> written = new ArrayList<>();
        if (result.datasetGvi() != null) {
            written.add(write(dir.resolve(datasetLabel + "_fig1_composition.svg"),
                    compositionFigure(result, datasetLabel)));
            written.add(write(dir.resolve(datasetLabel + "_fig2_indices.svg"),
                    indexPanelFigure(result, datasetLabel)));
        }
        if (result.sensitivity() != null && !result.sensitivity().isEmpty()) {
            written.add(write(dir.resolve(datasetLabel + "_fig3_sensitivity.svg"),
                    sensitivityFigure(result, datasetLabel)));
        }
        return written;
    }

    private static Path write(Path p, String svg) throws IOException {
        Files.writeString(p, svg);
        return p;
    }

    /**
     * Figure 1 -- what the score is made of. A horizontal stacked bar of each index's actual
     * contribution, with excluded indices shown at the end as hatched blocks sized by the weight
     * they would have carried. The reader sees both the score and how much of the intended picture
     * it rests on.
     */
    static String compositionFigure(PipelineResult result, String label) {
        GviResult gvi = result.datasetGvi();
        List<GviComponent> comps = new ArrayList<>(gvi.components());
        comps.sort((a, b) -> Double.compare(b.contribution(), a.contribution()));
        List<IndexKey> excluded = gvi.excludedIndices() == null ? List.of() : gvi.excludedIndices();

        // Derive the canvas height from the same layout arithmetic the body uses, rather than from a
        // constant that has to be remembered whenever a row is added. Getting this wrong does not
        // fail loudly -- SVG simply clips whatever falls outside the viewBox, so the excluded-index
        // block silently vanished from the figure while every value in it was still correct.
        double barTop = gvi.comparable() ? 52 : 66;
        double legendTop = barTop + 26 + 34;
        double h = legendTop + 12 + comps.size() * 11
                + (excluded.isEmpty() ? 0 : 6 + 12 + excluded.size() * 11) + 14;
        StringBuilder s = header(COLUMN_W, h);

        s.append(text(14, 22, 11, "600", INK, esc(label)));
        s.append(text(14, 36, 8, "400", MUTED,
                String.format(Locale.ROOT, "Composite GVI = %.4f  |  %s", gvi.gvi(), gvi.coverageSummary())));

        // A score built from a sliver of the weighting scheme is not comparable to a fully-populated
        // one, and by value alone the two are indistinguishable. Say so on the figure itself, so the
        // caveat travels with the image into a manuscript rather than living only in the report text.
        if (!gvi.comparable()) {
            s.append(rect(14, 42, COLUMN_W - 28, 15, "#fdf1e7", ACCENT, 0.7));
            s.append(text(19, 52.5, 6.8, "600", ACCENT,
                    "Not comparable across datasets: too little of the weighting scheme had usable data."));
        }

        // Stacked contribution bar, scaled so the drawn width is the GVI on a 0-1 axis.
        double barY = barTop, barH = 26, barX = 14, barW = COLUMN_W - 28;
        s.append(rect(barX, barY, barW, barH, "none", RULE, 0.8));
        double x = barX;
        for (GviComponent c : comps) {
            double w = c.contribution() * barW;
            if (w <= 0) continue;
            s.append(rect(x, barY, w, barH, SCORED, "#ffffff", 0.6, opacityFor(comps.indexOf(c), comps.size())));
            if (w > 20) {
                s.append(text(x + w / 2, barY + barH / 2 + 3, 6.5, "600", "#ffffff", c.key().label(), "middle"));
            }
            x += w;
        }
        s.append(line(barX + gvi.gvi() * barW, barY - 4, barX + gvi.gvi() * barW, barY + barH + 4, ACCENT, 1.2));
        s.append(text(barX, barY + barH + 16, 7, "400", MUTED, "0.0"));
        s.append(text(barX + barW, barY + barH + 16, 7, "400", MUTED, "1.0", "end"));

        // Legend of contributions, then the excluded block.
        double ly = barY + barH + 34;
        s.append(text(14, ly, 7.5, "600", INK, "Contributing indices"));
        ly += 12;
        for (GviComponent c : comps) {
            s.append(rect(14, ly - 6, 8, 8, SCORED, "none", 0, opacityFor(comps.indexOf(c), comps.size())));
            s.append(text(27, ly, 7, "400", INK, c.key().label()));
            s.append(text(120, ly, 7, "400", MUTED,
                    String.format(Locale.ROOT, "w=%.3f", c.effectiveWeight())));
            s.append(text(COLUMN_W - 14, ly, 7, "400", INK,
                    String.format(Locale.ROOT, "%.4f", c.contribution()), "end"));
            ly += 11;
        }
        if (!excluded.isEmpty()) {
            ly += 6;
            s.append(text(14, ly, 7.5, "600", ACCENT, "Excluded from the score (reported, not scored)"));
            ly += 12;
            for (IndexKey k : excluded) {
                s.append(hatched(14, ly - 6, 8, 8));
                s.append(text(27, ly, 7, "400", MUTED, k.label()));
                ly += 11;
            }
        }
        return s.append("</svg>\n").toString();
    }

    /**
     * Figure 2 -- every index at a glance, scored or not. A small-multiple panel rather than a
     * table: the eye compares bar lengths far faster than it reads numbers, and the hatching makes
     * the excluded ones unmistakable without needing a legend lookup.
     */
    static String indexPanelFigure(PipelineResult result, String label) {
        Map<IndexKey, Double> normalized = new LinkedHashMap<>();
        Map<IndexKey, Double> raw = new LinkedHashMap<>();
        for (GviComponent c : result.datasetGvi().components()) {
            normalized.put(c.key(), c.normalizedValue());
            raw.put(c.key(), c.rawValue());
        }
        List<IndexKey> excluded = result.datasetGvi().excludedIndices() == null
                ? List.of() : result.datasetGvi().excludedIndices();
        for (IndexKey k : excluded) {
            IndexResult r = result.datasetIndices().get(k);
            if (r == null) r = result.populationIndices().get(k);
            if (r != null) raw.put(k, r.primaryValue());
        }

        List<IndexKey> order = new ArrayList<>(normalized.keySet());
        for (IndexKey k : excluded) if (!order.contains(k)) order.add(k);

        double rowH = 17, top = 56;
        double h = top + order.size() * rowH + 26;
        StringBuilder s = header(COLUMN_W, h);
        s.append(text(14, 22, 11, "600", INK, esc(label)));
        s.append(text(14, 36, 8, "400", MUTED, "Normalized index values (0-1). Hatched = excluded from the composite."));

        double barX = 96, barW = COLUMN_W - barX - 58;
        for (int i = 0; i < order.size(); i++) {
            IndexKey k = order.get(i);
            double y = top + i * rowH;
            boolean isExcluded = excluded.contains(k);
            s.append(text(14, y + 9, 7.5, isExcluded ? "400" : "600", isExcluded ? MUTED : INK, k.label()));
            s.append(rect(barX, y, barW, 11, "#f4f4f5", RULE, 0.5));
            Double v = normalized.get(k);
            if (v != null && v > 0) {
                s.append(rect(barX, y, v * barW, 11, SCORED, "none", 0));
            } else if (isExcluded) {
                s.append(hatched(barX, y, barW, 11));
            }
            Double rv = raw.get(k);
            if (rv != null) {
                s.append(text(COLUMN_W - 14, y + 9, 6.8, "400", isExcluded ? MUTED : INK, fmt(rv), "end"));
            }
        }
        s.append(line(barX, top - 4, barX, top + order.size() * rowH, RULE, 0.6));
        return s.append("</svg>\n").toString();
    }

    /**
     * Figure 3 -- a tornado of how far the score moves when each weight is perturbed. This is the
     * honest companion to any GVI number in a paper: the weights are asserted rather than fitted,
     * so the reader is entitled to see how much the conclusion depends on them.
     */
    static String sensitivityFigure(PipelineResult result, String label) {
        List<SensitivityResult> sens = new ArrayList<>(result.sensitivity());
        sens.sort((a, b) -> Double.compare(span(b), span(a)));

        double rowH = 16, top = 58;
        double h = top + sens.size() * rowH + 30;
        StringBuilder s = header(COLUMN_W, h);
        s.append(text(14, 22, 11, "600", INK, esc(label)));
        s.append(text(14, 36, 8, "400", MUTED, "Sensitivity of GVI to a ±20% change in each weight"));

        double base = sens.isEmpty() ? 0 : sens.get(0).baseGvi();
        double maxSpan = sens.stream().mapToDouble(PublicationFigures::span).max().orElse(0.001);
        double cx = COLUMN_W / 2 + 20, half = (COLUMN_W - cx - 20);
        double scale = maxSpan == 0 ? 0 : half / maxSpan;

        s.append(line(cx, top - 6, cx, top + sens.size() * rowH + 2, INK, 0.8));
        s.append(text(cx, top - 10, 6.5, "400", MUTED, String.format(Locale.ROOT, "base %.4f", base), "middle"));

        for (int i = 0; i < sens.size(); i++) {
            SensitivityResult r = sens.get(i);
            double y = top + i * rowH;
            double lo = (r.gviAtLowWeight() - r.baseGvi()) * scale;
            double hi = (r.gviAtHighWeight() - r.baseGvi()) * scale;
            double x0 = cx + Math.min(lo, hi), x1 = cx + Math.max(lo, hi);
            s.append(text(14, y + 9, 7.5, "400", INK, r.key().label()));
            s.append(rect(x0, y + 1, Math.max(x1 - x0, 0.6), 10, ACCENT, "none", 0, 0.85));
            s.append(text(COLUMN_W - 14, y + 9, 6.5, "400", MUTED,
                    String.format(Locale.ROOT, "%+.4f", span(r)), "end"));
        }
        return s.append("</svg>\n").toString();
    }

    private static double span(SensitivityResult r) {
        return Math.abs(r.gviAtHighWeight() - r.gviAtLowWeight());
    }

    // ---- SVG primitives -------------------------------------------------------------------

    private static StringBuilder header(double w, double h) {
        StringBuilder s = new StringBuilder(4096);
        s.append("<?xml version='1.0' encoding='UTF-8'?>\n");
        s.append("<svg xmlns='http://www.w3.org/2000/svg' width='").append(f(w)).append("' height='").append(f(h))
         .append("' viewBox='0 0 ").append(f(w)).append(' ').append(f(h)).append("' font-family=\"").append(FONT).append("\">\n");
        // Hatching is drawn as explicit diagonals rather than an SVG <pattern>. Patterns are the
        // tidier markup, but several common rasterizers (ImageMatick/rsvg among them) render a
        // pattern fill as flat black, which turns "excluded" into a solid blob -- the exact opposite
        // of the distinction the hatching exists to make. Explicit lines render identically
        // everywhere, and these figures are meant to survive whatever toolchain a journal uses.
        s.append(rect(0, 0, w, h, "#ffffff", "none", 0));
        return s;
    }

    /**
     * A hatched block: white ground, thin outline, and 45-degree diagonals whose endpoints are
     * computed to land exactly on the block's edges.
     * <p>
     * The obvious implementations both fail in practice. An SVG {@code <pattern>} fill renders as
     * flat black in several common rasterizers, turning "excluded" into a solid blob -- the opposite
     * of the distinction the hatching exists to draw. Wrapping the diagonals in a {@code clipPath}
     * fails the same way: the clip is ignored and the lines spill across neighbouring labels. So the
     * clipping is done arithmetically here, and the emitted markup is nothing but plain lines, which
     * every renderer treats identically. These figures have to survive whatever toolchain a journal
     * puts them through.
     */
    private static String hatched(double x, double y, double w, double h) {
        StringBuilder b = new StringBuilder();
        b.append(rect(x, y, w, h, "#ffffff", MUTED, 0.5));
        // Each diagonal runs parallel to (1,-1) from (x+d, y+h); solve for the parameter range that
        // keeps it inside [x, x+w] and draw only that segment.
        for (double d = -h; d < w; d += 3.2) {
            double tLo = Math.max(0.0, -d / h);
            double tHi = Math.min(1.0, (w - d) / h);
            if (tHi <= tLo) continue;
            double x0 = x + d + tLo * h, y0 = y + h - tLo * h;
            double x1 = x + d + tHi * h, y1 = y + h - tHi * h;
            b.append(line(x0, y0, x1, y1, MUTED, 0.7));
        }
        // Redraw the outline so the diagonals sit inside a crisp edge.
        b.append(rect(x, y, w, h, "none", MUTED, 0.5));
        return b.toString();
    }

    private static String rect(double x, double y, double w, double h, String fill, String stroke, double sw) {
        return rect(x, y, w, h, fill, stroke, sw, 1.0);
    }

    private static String rect(double x, double y, double w, double h, String fill, String stroke, double sw, double op) {
        return "<rect x='" + f(x) + "' y='" + f(y) + "' width='" + f(w) + "' height='" + f(h)
                + "' fill='" + fill + "'" + (sw > 0 ? " stroke='" + stroke + "' stroke-width='" + f(sw) + "'" : "")
                + (op < 1.0 ? " opacity='" + f(op) + "'" : "") + "/>";
    }

    private static String line(double x1, double y1, double x2, double y2, String c, double w) {
        return "<line x1='" + f(x1) + "' y1='" + f(y1) + "' x2='" + f(x2) + "' y2='" + f(y2)
                + "' stroke='" + c + "' stroke-width='" + f(w) + "'/>";
    }

    private static String text(double x, double y, double size, String weight, String fill, String content) {
        return text(x, y, size, weight, fill, content, "start");
    }

    private static String text(double x, double y, double size, String weight, String fill, String content, String anchor) {
        return "<text x='" + f(x) + "' y='" + f(y) + "' font-size='" + f(size) + "' font-weight='" + weight
                + "' fill='" + fill + "' text-anchor='" + anchor + "'>" + esc(content) + "</text>";
    }

    /** Opacity ramp so adjacent stacked segments stay separable without needing nine hues. */
    private static double opacityFor(int i, int n) {
        return n <= 1 ? 1.0 : 0.45 + 0.55 * (1.0 - (double) i / n);
    }

    private static String fmt(double v) {
        double a = Math.abs(v);
        if (a != 0 && (a < 1e-3 || a >= 1e5)) return String.format(Locale.ROOT, "%.2e", v);
        if (a >= 100) return String.format(Locale.ROOT, "%.1f", v);
        return String.format(Locale.ROOT, "%.4f", v);
    }

    private static String f(double d) {
        return String.format(Locale.ROOT, "%.2f", d);
    }

    private static String esc(String s) {
        return s == null ? "" : s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
                .replace("'", "&apos;").replace("\"", "&quot;");
    }
}
