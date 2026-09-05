package org.gvi.algorithms.cai;

import org.gvi.core.exception.GviInputException;
import org.gvi.core.io.CodonUsageTableReader;
import org.gvi.core.model.CodonUsageTable;

import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 * Precomputed host codon-usage reference tables, bundled directly into the
 * jar (no network access needed at run time) as the same kind of built-in
 * default CodonW itself ships (its own default reference set, historically
 * E. coli high-expression genes) -- so CAI has a real fallback for common
 * hosts without requiring the user to source and format their own
 * {@code --codon-usage} CSV first.
 * <p>
 * <b>Provenance:</b> every table here is computed by THIS project directly
 * from real, complete CDS FASTA files downloaded from NCBI RefSeq
 * (each species' current reference genome assembly, resolved via the NCBI
 * Datasets API on 2026-08-17) -- every codon in every annotated coding
 * sequence for that assembly is tallied and converted to a per-thousand
 * frequency, the same computation the Kazusa Codon Usage Database itself
 * performs, just against a current, much larger RefSeq annotation rather
 * than an older static snapshot (an earlier pass of this feature used
 * Kazusa's own precomputed tables directly; this pass recomputes from the
 * raw CDS sequences instead, which mattered most for species Kazusa had
 * very little data for -- see buffalo below, a >1000x sample increase):
 * <ul>
 *   <li>{@code human} -- <i>Homo sapiens</i>, RefSeq GCF_000001405.40
 *       (GRCh38.p14), 146,337 CDS's / 99,625,497 codons</li>
 *   <li>{@code mouse} -- <i>Mus musculus</i>, RefSeq GCF_000001635.27
 *       (GRCm39), 98,005 CDS's / 66,614,578 codons</li>
 *   <li>{@code pig} / {@code wild_boar} -- <i>Sus scrofa</i>, RefSeq
 *       GCF_054392235.1, 111,261 CDS's / 81,974,901 codons. Domestic pig
 *       and wild boar are the SAME species taxonomically (one NCBI taxid,
 *       9823, no separate subspecies-level assembly/table exists) -- so
 *       {@code wild_boar} is an alias for exactly the same real data as
 *       {@code pig}, not a distinct fetch.</li>
 *   <li>{@code cattle} -- <i>Bos taurus</i>, RefSeq GCF_002263795.3
 *       (ARS-UCD2.0), 64,900 CDS's / 44,401,208 codons</li>
 *   <li>{@code buffalo} -- <i>Bubalus bubalis</i> (water buffalo), RefSeq
 *       GCF_019923935.1, 64,670 CDS's / 46,756,930 codons</li>
 *   <li>{@code sheep} -- <i>Ovis aries</i>, RefSeq GCF_016772045.2
 *       (ARS-UI_Ramb_v3.0), 76,912 CDS's / 56,263,957 codons</li>
 *   <li>{@code goat} -- <i>Capra hircus</i>, RefSeq GCF_001704415.2
 *       (ARS1.2), 42,815 CDS's / 28,640,502 codons</li>
 *   <li>{@code horse} -- <i>Equus caballus</i>, RefSeq GCF_041296265.1
 *       (TB-T2T), 97,963 CDS's / 71,968,894 codons</li>
 *   <li>{@code ecoli} -- <i>Escherichia coli</i> str. K-12 substr. MG1655,
 *       RefSeq GCF_000005845.2 (ASM584v2), 4,318 CDS's / 1,342,295 codons</li>
 *   <li>{@code aedes_aegypti} -- <i>Aedes aegypti</i>, RefSeq
 *       GCF_002204515.2 (AaegL5.0), 28,317 CDS's / 20,227,358 codons</li>
 * </ul>
 * A bundled table is always a real default, never a substitute for a
 * dataset-specific host reference built from the ACTUAL host organism's own
 * highly-expressed genes -- prefer a genuine {@code --codon-usage} table
 * whenever you have one; use these when you don't and one of these species
 * is a reasonable match for the actual host under study.
 */
public final class BundledCodonUsageTables {

    private static final Map<String, String> RESOURCE_BY_SPECIES = new LinkedHashMap<>();

    static {
        RESOURCE_BY_SPECIES.put("human", "/codon_usage/human.csv");
        RESOURCE_BY_SPECIES.put("mouse", "/codon_usage/mouse.csv");
        RESOURCE_BY_SPECIES.put("pig", "/codon_usage/pig.csv");
        RESOURCE_BY_SPECIES.put("wild_boar", "/codon_usage/pig.csv"); // same species (Sus scrofa) -- see class javadoc
        RESOURCE_BY_SPECIES.put("cattle", "/codon_usage/cattle.csv");
        RESOURCE_BY_SPECIES.put("buffalo", "/codon_usage/buffalo.csv");
        RESOURCE_BY_SPECIES.put("sheep", "/codon_usage/sheep.csv");
        RESOURCE_BY_SPECIES.put("goat", "/codon_usage/goat.csv");
        RESOURCE_BY_SPECIES.put("horse", "/codon_usage/horse.csv");
        RESOURCE_BY_SPECIES.put("ecoli", "/codon_usage/ecoli.csv");
        RESOURCE_BY_SPECIES.put("aedes_aegypti", "/codon_usage/aedes_aegypti.csv");
    }

    private BundledCodonUsageTables() {
    }

    public static java.util.Set<String> availableSpecies() {
        return RESOURCE_BY_SPECIES.keySet();
    }

    public static CodonUsageTable load(String species) {
        String key = species.toLowerCase().trim();
        String resource = RESOURCE_BY_SPECIES.get(key);
        if (resource == null) {
            throw new GviInputException("Unknown bundled codon usage species '" + species + "'; available: "
                    + String.join(", ", RESOURCE_BY_SPECIES.keySet()));
        }
        try (InputStream in = BundledCodonUsageTables.class.getResourceAsStream(resource)) {
            if (in == null) {
                throw new GviInputException("Bundled codon usage resource '" + resource + "' is missing from the jar");
            }
            return CodonUsageTableReader.read(new InputStreamReader(in, StandardCharsets.UTF_8), "bundled:" + key);
        } catch (IOException e) {
            throw new GviInputException("Could not read bundled codon usage table '" + key + "': " + e.getMessage(), e);
        }
    }
}
