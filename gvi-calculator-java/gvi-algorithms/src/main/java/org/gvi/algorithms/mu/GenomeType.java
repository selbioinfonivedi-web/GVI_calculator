package org.gvi.algorithms.mu;

/**
 * A pathogen's genome chemistry, supplied by the user (there is nothing in a
 * FASTA alignment itself that reveals whether the sequence came from an RNA
 * or DNA genome). Determines which {@link MuReferenceTable} band set -- and
 * which category wording -- a computed {@code mu} is classified against.
 * <p>
 * Without this, every mu estimate was classified against RNA-virus reference
 * points (Influenza, SARS-CoV-2, VSV) regardless of the actual pathogen: a
 * dsDNA poxvirus or a non-coding rRNA locus would still come back labelled
 * "RNA (low-fidelity)" with an "immune escape 6-18 months" narrative that
 * doesn't apply to either. {@link #UNSPECIFIED} preserves that historical
 * default (for backward compatibility and when the user genuinely doesn't
 * know/care), but the returned category text says plainly that it's an
 * RNA-calibrated assumption rather than presenting it as fact.
 */
public enum GenomeType {
    RNA,
    DNA,
    UNSPECIFIED
}
