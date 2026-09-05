package org.gvi.composite;

/**
 * The 9 scalar components the composite GVI combines (Section 5.9): the 8
 * named indices, with CAI split into its two sub-components (CAI and GC
 * Content Deviation) since the spec weights them together but they are
 * computed and normalized separately.
 */
public enum IndexKey {
    MU("mu"),
    RE("Re"),
    PI("pi"),
    MB("MB"),
    DNDS("dN/dS"),
    GD("GD"),
    CAI("CAI"),
    GC("GC_Deviation"),
    RI("RI");

    private final String label;

    IndexKey(String label) {
        this.label = label;
    }

    public String label() {
        return label;
    }

    public static IndexKey fromLabel(String label) {
        for (IndexKey k : values()) {
            if (k.label.equalsIgnoreCase(label)) return k;
        }
        throw new IllegalArgumentException("Unknown index label '" + label + "'; expected one of: mu, Re, pi, MB, dN/dS, GD, CAI, GC_Deviation, RI");
    }
}
