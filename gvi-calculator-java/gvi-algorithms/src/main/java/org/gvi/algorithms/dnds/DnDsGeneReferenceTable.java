package org.gvi.algorithms.dnds;

import java.util.List;
import java.util.Optional;

/**
 * Gene-specific dN/dS reference bands from build spec Section 5.4. Matched
 * against a gene name by case-insensitive keyword, since GFF3 gene naming
 * conventions vary across pathogens. Falls back to the generic selection
 * classification (Section 5.2) when no keyword matches.
 */
public final class DnDsGeneReferenceTable {

    public record GeneBand(String geneCategory, String dnDsRange, String selectionType,
                            String biologicalRole, String pandemicRisk, List<String> keywords) {
        boolean matches(String geneName) {
            String lower = geneName.toLowerCase();
            return keywords.stream().anyMatch(lower::contains);
        }
    }

    public static final List<GeneBand> BANDS = List.of(
            new GeneBand("RNA Polymerase", "<0.3", "Strong purifying", "Replication fidelity; essential",
                    "Low (resistance rare)", List.of("pol", "rdrp", "polymerase", "l gene")),
            new GeneBand("Protease", "0.3-0.6", "Moderate purifying", "Substrate specificity; essential",
                    "Low-moderate", List.of("protease", "3clpro", "mpro", "nsp5")),
            new GeneBand("Nucleocapsid", "0.6-1.0", "Neutral-weak purifying", "Packaging; less constrained",
                    "Moderate", List.of("nucleocapsid", "capsid", "n gene", "nucleoprotein", "np")),
            new GeneBand("Spike / Surface", "1.0-3.0", "Positive (often) / relaxed", "Host cell entry; immune target",
                    "High (escape hotspot)", List.of("spike", "surface", "s gene", "hemagglutinin", "ha", "gp120", "env glycoprotein")),
            new GeneBand("Envelope", "1.5-2.5", "Positive", "Membrane incorporation; variable",
                    "High (immune pressure)", List.of("envelope", "e gene", "matrix"))
    );

    private DnDsGeneReferenceTable() {
    }

    public static Optional<GeneBand> match(String geneName) {
        if (geneName == null || geneName.isBlank()) return Optional.empty();
        return BANDS.stream().filter(b -> b.matches(geneName)).findFirst();
    }

    /** Generic fallback classification (Section 5.2) when the gene name doesn't match a known category. */
    public static String genericClassification(double omega) {
        if (omega < 0.95) return "Purifying selection (dN/dS < 1): functional constraint";
        if (omega > 1.05) return "Positive selection (dN/dS > 1): adaptive advantage";
        return "Neutral evolution (dN/dS ~= 1): relaxed selection";
    }
}
