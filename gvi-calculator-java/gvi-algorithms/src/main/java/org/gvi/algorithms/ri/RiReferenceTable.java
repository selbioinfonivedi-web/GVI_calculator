package org.gvi.algorithms.ri;

import java.util.List;

/**
 * Reference thresholds from build spec Section 7.4.
 * <p>
 * <b>These bands describe RECOMBINATION only -- never reassortment.</b> Both
 * of {@link RecombinationIndexCalculator}'s methods operate on a single
 * locus/alignment: the PHI test measures intra-locus homoplasy structure, and
 * the isolate-classification variant counts isolates an upstream tool already
 * flagged as recombinant. Neither compares segments against each other, which
 * is what detecting reassortment requires -- that is
 * {@link ReassortmentIndexCalculator}'s job (driven by {@code --segment}), and
 * it carries its own separate category vocabulary.
 * <p>
 * Deviation from the source spec, deliberate: Section 5.7's own wording labels
 * the top band "reassortment-prone / segmented" with "influenza pandemic
 * strains" as its example. That label was reachable from a single-alignment
 * PHI test that cannot observe reassortment at all, and it fired exactly that
 * way on real data (a Bluetongue run was labelled "Reassortment-prone /
 * segmented (Very high mixing)" from one segment's alignment). The numeric
 * thresholds are preserved as specified; only the wording is corrected to
 * describe what was actually measured.
 */
public final class RiReferenceTable {

    public record Band(double upperBound, String pathogenContext, String interpretation, String controlImplication, String example) {
    }

    public static final List<Band> BANDS = List.of(
            new Band(0.02, "Single-lineage outbreak", "No/minimal recombination", "Simple phylogeny; standard analysis", "Early Ebola outbreaks"),
            new Band(0.10, "Endemic multi-type circulation", "Occasional recombinants", "Exclude from phylodynamic analysis", "Dengue endemic countries"),
            new Band(0.30, "High co-circulation", "Common recombination", "Partition analysis; expect chimeric variants", "HIV-1 circulating recombinant forms (CRFs)"),
            new Band(Double.POSITIVE_INFINITY, "Extensive recombination / mosaic genomes", "Very high mosaicism",
                    "Breakpoint-aware partitioning required; single-tree phylogenies unreliable. If this genome is segmented, "
                            + "reassortment is NOT assessed here -- run the segment-aware analysis (--segment) separately",
                    "Enteroviruses; HIV-1 inter-subtype recombinants")
    );

    private RiReferenceTable() {
    }

    public static Band classify(double ri) {
        for (Band b : BANDS) {
            if (ri <= b.upperBound()) return b;
        }
        return BANDS.get(BANDS.size() - 1);
    }
}
