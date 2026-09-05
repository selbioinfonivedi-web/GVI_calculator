package org.gvi.core.model;

/**
 * A single coding-sequence (CDS) region from a GFF3 annotation, 1-based
 * inclusive coordinates as per GFF3 convention. Used to extract in-frame
 * codons for dN/dS and CAI, and to attribute a gene name so gene-specific
 * reference tables (Section 5.5, 5.8) can be applied.
 */
public record GeneAnnotation(String geneName, long start, long end, char strand) {

    public GeneAnnotation {
        if (end < start) {
            throw new IllegalArgumentException("gene '" + geneName + "': end (" + end + ") < start (" + start + ")");
        }
        if (strand != '+' && strand != '-') {
            throw new IllegalArgumentException("gene '" + geneName + "': strand must be '+' or '-', got '" + strand + "'");
        }
    }

    public long length() {
        return end - start + 1;
    }
}
