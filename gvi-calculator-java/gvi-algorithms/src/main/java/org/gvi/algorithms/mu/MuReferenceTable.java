package org.gvi.algorithms.mu;

import java.util.List;

/** Reference thresholds from build spec Section 1.4. */
public final class MuReferenceTable {

    public record Band(double upperBound, String pathogenExample, String polymeraseType,
                        String immuneEscapeRate, String clinicalImplication) {
    }

    /** RNA virus reference points -- the original (and, for RNA genomes, still appropriate) band set. */
    public static final List<Band> BANDS = List.of(
            new Band(1e-4, "Measles, Polio (DNA-like)", "High-fidelity", "Years-decades", "Lifelong immunity; vaccine durable"),
            new Band(1e-3, "Dengue, WNV, Zika", "RNA (moderate)", "1-5 years", "Serotype circulation; periodic outbreaks"),
            new Band(1e-2, "Influenza, SARS-CoV-2, HIV", "RNA (low-fidelity)", "6-18 months", "Annual/seasonal variants; escape mutations"),
            new Band(Double.POSITIVE_INFINITY, "Vesicular Stomatitis Virus (VSV)", "Ultra high error", "Days-weeks", "Quasispecies; rapid resistance; hard to target")
    );

    /**
     * DNA genome reference points. Real DNA viruses -- lacking RNA-dependent-RNA-polymerase's error rate, and
     * usually (except small ssDNA genomes replicated without proofreading) using a high-fidelity host or viral
     * DNA polymerase -- sit several orders of magnitude below the RNA bands above. A DNA-genome mu landing above
     * ~1e-5 subs/site/yr essentially never reflects a real substitution rate at that scale; it is a strong signal
     * of lineage-mixing, misalignment, or a marker under such extreme diversifying selection that the clock
     * assumption itself has broken down -- so the top band says so explicitly instead of naming a virus example.
     */
    public static final List<Band> DNA_BANDS = List.of(
            new Band(1e-8, "Herpesvirus, large Poxvirus (e.g. Variola, Vaccinia)", "High-fidelity dsDNA (proofreading)", "Decades", "Antigenically stable; live long-term in host/reservoir"),
            new Band(1e-7, "Adenovirus, Papillomavirus, African Swine Fever Virus", "Moderate-fidelity DNA", "Years", "Slow antigenic drift; genotype-level surveillance adequate"),
            new Band(1e-6, "Capripoxvirus (LSDV, Sheep/Goat Pox), Baculovirus", "Lower-fidelity DNA / larger indel tolerance", "Years", "Strain-level drift detectable; still far slower than RNA viruses"),
            new Band(1e-5, "Circovirus, Parvovirus, Geminivirus (small ssDNA, rolling-circle)", "Small ssDNA (no proofreading)", "Months-years", "Can approach RNA-virus-like turnover despite being DNA"),
            new Band(Double.POSITIVE_INFINITY, "No known DNA virus replicates this fast", "Implausible for a DNA genome", "N/A", "Likely lineage-mixing, misalignment, or a broken clock assumption -- verify the input before trusting this rate")
    );

    /**
     * Upper edge of the fastest band {@link #DNA_BANDS} treats as a real DNA substitution rate (small ssDNA
     * genomes replicated without proofreading -- circoviruses, parvoviruses, geminiviruses). Above this, the
     * table's own top band says "No known DNA virus replicates this fast".
     */
    public static final double MAX_CREDIBLE_DNA_RATE = 1e-5;

    /**
     * True when the reference table classifies this rate as not physically credible for the declared genome
     * chemistry, i.e. the terminal DNA band whose clinical implication is literally "verify the input before
     * trusting this rate".
     * <p>
     * That verdict was being computed and printed while the rate it condemned still scored at full weight.
     * Lumpy Skin Disease is the case that exposed it: mu came back at 1.85e-4 sub/site/yr for a capripoxvirus
     * whose real rate the table puts around 1e-6, so the category read "Implausible for a DNA genome (No known
     * DNA virus replicates this fast)" -- and because the composite normalizes mu against a 1e-4 DNA ceiling,
     * that same implausible rate saturated at 1.0 and became the single largest contributor to the headline
     * score. The more impossible the rate, the harder it pushed the GVI up.
     * <p>
     * Note this is a distinct failure from the ones {@code QualityGate} already caught. The date-randomization
     * test <em>passed</em> here (empirical p=0.0099) and R² was 0.47: there is genuine temporal signal, the
     * dates really do explain the divergence. What is wrong is the magnitude, and no amount of temporal-signal
     * evidence makes a rate two orders of magnitude above the biological ceiling credible. Only a scale check
     * against the organism's own chemistry catches it.
     * <p>
     * RNA genomes have no equivalent ceiling: {@link #BANDS}' top band (VSV, "ultra high error") is a real,
     * observed regime rather than an impossibility, so nothing is gated for RNA.
     */
    public static boolean isImplausible(double muSubPerSiteYear, GenomeType genomeType) {
        return genomeType == GenomeType.DNA && muSubPerSiteYear > MAX_CREDIBLE_DNA_RATE;
    }

    private MuReferenceTable() {
    }

    /** @deprecated use {@link #classify(double, GenomeType)}; kept for callers that haven't specified a genome type. */
    @Deprecated
    public static Band classify(double muSubPerSiteYear) {
        return classify(muSubPerSiteYear, GenomeType.UNSPECIFIED);
    }

    public static Band classify(double muSubPerSiteYear, GenomeType genomeType) {
        List<Band> bands = genomeType == GenomeType.DNA ? DNA_BANDS : BANDS;
        for (Band b : bands) {
            if (muSubPerSiteYear <= b.upperBound()) return b;
        }
        return bands.get(bands.size() - 1);
    }

    /**
     * Category text for a classified band, with an explicit caveat when the genome type wasn't actually
     * specified -- the RNA bands above are a real reference table for RNA genomes, not a safe default assumption
     * for whatever pathogen happens to be running through the pipeline.
     */
    public static String describe(Band band, GenomeType genomeType) {
        String base = band.polymeraseType() + " (" + band.pathogenExample() + "); immune escape " + band.immuneEscapeRate();
        if (genomeType == GenomeType.UNSPECIFIED) {
            return base + " -- assumes an RNA genome (pass --genome-type dna if this pathogen's genome is DNA; "
                    + "these examples/timescales will be wrong otherwise)";
        }
        return base;
    }
}
