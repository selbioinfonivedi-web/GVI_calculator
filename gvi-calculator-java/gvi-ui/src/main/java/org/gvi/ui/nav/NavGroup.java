package org.gvi.ui.nav;

/**
 * Top-level groups in the left navigation rail.
 * <p>
 * Reduced from five to two. PROJECT and ANALYSIS held only pages that explained the work happened
 * elsewhere, and INPUT/GVI each ended up with a single entry once the one- and two-field pages were
 * folded into the step they qualify -- a heading above a single item is noise, not structure. What
 * remains matches how the tool is actually used: set the run up, then read the result.
 */
public enum NavGroup {
    WORKFLOW("Workflow"),
    RESULTS("Results");

    private final String label;

    NavGroup(String label) {
        this.label = label;
    }

    public String label() {
        return label;
    }
}
