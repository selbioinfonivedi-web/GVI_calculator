package org.gvi.web;

import java.util.List;
import java.util.Map;

/**
 * The shape the browser renders.
 * <p>
 * This is deliberately NOT the raw {@link org.gvi.cli.PipelineResult}. Two things the pipeline reports
 * as loose text are promoted to first-class, structured fields here, because a browser UI that buries
 * them reproduces the exact failure mode this project has been fixing:
 * <ul>
 *   <li>{@code comparable}/{@code coverageSummary} -- a GVI built from two indices is renormalized to
 *       [0,1] and is numerically indistinguishable from one built from nine. The flag has to travel
 *       with the number so the UI can refuse to present it as a rankable score.</li>
 *   <li>{@code exclusions} -- which indices quality gating dropped, and why. In the text report these
 *       are strings in a "Skipped" list below the result; here they are split out so the UI can show
 *       each excluded index next to the score it did not contribute to.</li>
 * </ul>
 */
public class AnalyzeResponse {
    public Double gvi;
    public boolean comparable;
    public String coverageSummary;
    public double effectiveWeightSum;
    public Double betaGenomic;
    public boolean betaExcludesRe;

    public List<Component> components;
    public List<Exclusion> exclusions;
    public List<String> otherSkipped;
    public List<String> warnings;
    public List<Sensitivity> sensitivity;

    public Map<String, IndexView> indices;
    public Summary summary;
    public String reportText;

    public static class Component {
        public String key;
        public double rawValue;
        public double normalizedValue;
        public double effectiveWeight;
        public double contribution;
    }

    /** One index quality gating removed from the weighted score, with the reason shown verbatim. */
    public static class Exclusion {
        public String key;
        public String reason;
    }

    public static class Sensitivity {
        public String key;
        public double low;
        public double high;
        public double spread;
    }

    /** An index's reported value and interpretation, whether or not it was allowed to score. */
    public static class IndexView {
        public String indexName;
        public double value;
        public String category;
        public List<String> diagnostics;
        public boolean scored;
        /**
         * Present only where the estimator produces one -- currently the birth-death Re profile.
         * Null elsewhere, which the UI must render as "no interval", never as zero.
         */
        public Interval interval;
    }

    /** A two-sided interval with its confidence level, e.g. Re's 95% likelihood-ratio interval. */
    public static class Interval {
        public double lower;
        public double upper;
        public double level;
    }

    public static class Summary {
        public int sequences;
        public int alignmentLength;
        public double gapFraction;
        public String referenceId;
        public String organismClass;
        public String genomeType;
        /** The generation time the run actually used, and where it came from. */
        public Double generationTimeDays;
        public String generationTimeSource;
    }
}
