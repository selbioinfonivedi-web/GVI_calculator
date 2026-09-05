package org.gvi.core.model;

/**
 * What kind of organism the alignment came from. Supplied by the user -- as
 * with {@code GenomeType}, nothing in a FASTA alignment itself reveals it.
 * <p>
 * This exists because several indices are only biologically meaningful for
 * some organism classes, and were previously applied uniformly:
 * <ul>
 *   <li><b>CAI</b> against a <em>host</em> codon table measures adaptation to
 *       the host's tRNA pool. That premise holds only for organisms that
 *       translate their proteins on host ribosomes -- i.e. viruses. Bacteria
 *       and eukaryotic parasites carry their own ribosomes and tRNA genes, so
 *       scoring them against a host (or an unrelated <em>E. coli</em>) table
 *       measures similarity to a third organism, not host adaptation.
 *       Note this is a statement about the <em>reference set</em>, not about
 *       CAI itself: Sharp &amp; Li (1987) define CAI against highly expressed
 *       genes of the same organism, and their original application was
 *       <em>E. coli</em> against its own ribosomal proteins. Scoring a
 *       bacterium against its own highly-expressed gene set would be valid;
 *       scoring it against its host is what is not.</li>
 *   <li><b>GC deviation</b> is reported with "HGT suspected" wording. Classical
 *       horizontal gene transfer (conjugation, transformation, transduction) is
 *       a real and central mechanism in bacteria -- <em>B. anthracis</em>
 *       virulence is plasmid-borne -- but is not how viruses typically acquire
 *       compositional shifts, where mutational pressure (APOBEC/ADAR editing)
 *       is the usual explanation.</li>
 * </ul>
 * {@link #UNSPECIFIED} preserves the historical behaviour (no gating) so
 * existing invocations keep working, but indices whose premise depends on the
 * organism class say plainly in their output that the premise was assumed
 * rather than confirmed.
 */
public enum OrganismClass {
    /** Translates on host ribosomes -- host-relative CAI is meaningful. */
    VIRUS,
    /** Own ribosomes and tRNA pool; HGT-competent. */
    BACTERIUM,
    /** Eukaryotic parasite (apicomplexan, kinetoplastid, helminth) -- own translation machinery. */
    EUKARYOTIC_PARASITE,
    /** Not supplied; no organism-class gating is applied. */
    UNSPECIFIED;

    /** True when this organism translates its own proteins rather than using host ribosomes. */
    public boolean hasOwnTranslationMachinery() {
        return this == BACTERIUM || this == EUKARYOTIC_PARASITE;
    }

    /** True when classical horizontal gene transfer is a documented mechanism for this class. */
    public boolean hgtCompetent() {
        return this == BACTERIUM;
    }

    public static OrganismClass parse(String raw) {
        if (raw == null || raw.isBlank()) return UNSPECIFIED;
        String v = raw.trim().toLowerCase(java.util.Locale.ROOT).replace('-', '_');
        return switch (v) {
            case "virus", "viral" -> VIRUS;
            case "bacterium", "bacteria", "bacterial" -> BACTERIUM;
            case "parasite", "eukaryotic_parasite", "eukaryote", "protozoan", "helminth" -> EUKARYOTIC_PARASITE;
            case "unspecified", "unknown" -> UNSPECIFIED;
            default -> throw new org.gvi.core.exception.GviInputException(
                    "Unrecognized --organism-class '" + raw + "'. Expected one of: virus, bacterium, parasite (or unspecified).");
        };
    }
}
