package org.gvi.core.model;

/**
 * One variant call from a VCF: a SNP (ref/alt both length 1) or an indel
 * (ref/alt differ in length). {@code depth} is the read depth at this site
 * (VCF INFO DP, or the first sample's FORMAT DP) used for the Mutation
 * Burden coverage filter (Section 5.4); {@code null} means depth was not
 * reported in the file.
 */
public record VariantRecord(String chrom, long pos, String ref, String alt, Integer depth, String sampleId) {

    public boolean isIndel() {
        return ref.length() != alt.length();
    }

    public boolean isSnp() {
        return ref.length() == 1 && alt.length() == 1;
    }
}
