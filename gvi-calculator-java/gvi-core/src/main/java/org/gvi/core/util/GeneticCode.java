package org.gvi.core.util;

import org.gvi.core.exception.GviComputationException;

import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * The standard nuclear genetic code (NCBI translation table 1). Bundled as
 * a resource-free static table (not a config file) because it is a fixed
 * biological constant, not something a user should be editing per Section
 * 5.5's Nei-Gojobori method and Section 5.8's CAI.
 */
public final class GeneticCode {

    public static final char STOP = '*';

    /** codon (uppercase, T not U) -> single-letter amino acid, '*' = stop */
    public static final Map<String, Character> CODON_TABLE;

    /** amino acid -> list of codons encoding it (includes '*' -> stop codons) */
    public static final Map<Character, List<String>> SYNONYMOUS_CODONS;

    static {
        Map<String, Character> table = new LinkedHashMap<>();
        put(table, "TTT", 'F'); put(table, "TTC", 'F'); put(table, "TTA", 'L'); put(table, "TTG", 'L');
        put(table, "CTT", 'L'); put(table, "CTC", 'L'); put(table, "CTA", 'L'); put(table, "CTG", 'L');
        put(table, "ATT", 'I'); put(table, "ATC", 'I'); put(table, "ATA", 'I'); put(table, "ATG", 'M');
        put(table, "GTT", 'V'); put(table, "GTC", 'V'); put(table, "GTA", 'V'); put(table, "GTG", 'V');
        put(table, "TCT", 'S'); put(table, "TCC", 'S'); put(table, "TCA", 'S'); put(table, "TCG", 'S');
        put(table, "CCT", 'P'); put(table, "CCC", 'P'); put(table, "CCA", 'P'); put(table, "CCG", 'P');
        put(table, "ACT", 'T'); put(table, "ACC", 'T'); put(table, "ACA", 'T'); put(table, "ACG", 'T');
        put(table, "GCT", 'A'); put(table, "GCC", 'A'); put(table, "GCA", 'A'); put(table, "GCG", 'A');
        put(table, "TAT", 'Y'); put(table, "TAC", 'Y'); put(table, "TAA", STOP); put(table, "TAG", STOP);
        put(table, "CAT", 'H'); put(table, "CAC", 'H'); put(table, "CAA", 'Q'); put(table, "CAG", 'Q');
        put(table, "AAT", 'N'); put(table, "AAC", 'N'); put(table, "AAA", 'K'); put(table, "AAG", 'K');
        put(table, "GAT", 'D'); put(table, "GAC", 'D'); put(table, "GAA", 'E'); put(table, "GAG", 'E');
        put(table, "TGT", 'C'); put(table, "TGC", 'C'); put(table, "TGA", STOP); put(table, "TGG", 'W');
        put(table, "CGT", 'R'); put(table, "CGC", 'R'); put(table, "CGA", 'R'); put(table, "CGG", 'R');
        put(table, "AGT", 'S'); put(table, "AGC", 'S'); put(table, "AGA", 'R'); put(table, "AGG", 'R');
        put(table, "GGT", 'G'); put(table, "GGC", 'G'); put(table, "GGA", 'G'); put(table, "GGG", 'G');
        CODON_TABLE = Collections.unmodifiableMap(table);

        Map<Character, List<String>> syn = new LinkedHashMap<>();
        for (Map.Entry<String, Character> e : table.entrySet()) {
            syn.computeIfAbsent(e.getValue(), k -> new ArrayList<>()).add(e.getKey());
        }
        Map<Character, List<String>> synUnmod = new LinkedHashMap<>();
        for (Map.Entry<Character, List<String>> e : syn.entrySet()) {
            synUnmod.put(e.getKey(), Collections.unmodifiableList(e.getValue()));
        }
        SYNONYMOUS_CODONS = Collections.unmodifiableMap(synUnmod);
    }

    private static void put(Map<String, Character> m, String codon, char aa) {
        m.put(codon, aa);
    }

    private GeneticCode() {
    }

    public static char translate(String codon) {
        String key = codon.toUpperCase().replace('U', 'T');
        Character aa = CODON_TABLE.get(key);
        if (aa == null) {
            throw new GviComputationException("Cannot translate codon '" + codon + "' (contains gap/ambiguity code or is not a valid triplet)");
        }
        return aa;
    }

    public static boolean isStandardCodon(String codon) {
        String key = codon.toUpperCase().replace('U', 'T');
        return CODON_TABLE.containsKey(key);
    }

    public static boolean isStop(String codon) {
        return translate(codon) == STOP;
    }
}
