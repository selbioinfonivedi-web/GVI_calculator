package org.gvi.algorithms.dnds.ml;

import org.gvi.core.util.GeneticCode;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * The 61 sense codons (all 64 standard codons minus the 3 stops), in a
 * fixed canonical order, plus the pairwise classification the GY94 codon
 * substitution model ({@link CodonModel}) needs: how many nucleotide
 * positions two codons differ at, which position if exactly one, whether
 * that single difference is a transition or transversion, and whether the
 * two codons encode the same amino acid (synonymous) or not.
 */
public final class CodonAlphabet {

    public static final List<String> SENSE_CODONS;
    public static final int SIZE;

    private static final Map<String, Integer> INDEX;

    static {
        List<String> codons = new ArrayList<>();
        for (Map.Entry<String, Character> e : GeneticCode.CODON_TABLE.entrySet()) {
            if (e.getValue() != GeneticCode.STOP) codons.add(e.getKey());
        }
        Collections.sort(codons);
        SENSE_CODONS = Collections.unmodifiableList(codons);
        SIZE = SENSE_CODONS.size(); // 61

        Map<String, Integer> idx = new HashMap<>();
        for (int i = 0; i < SENSE_CODONS.size(); i++) idx.put(SENSE_CODONS.get(i), i);
        INDEX = Collections.unmodifiableMap(idx);
    }

    private CodonAlphabet() {
    }

    public static int index(String codon) {
        Integer i = INDEX.get(codon.toUpperCase());
        if (i == null) {
            throw new IllegalArgumentException("'" + codon + "' is not one of the 61 standard sense codons");
        }
        return i;
    }

    public static boolean isSenseCodon(String codon) {
        return INDEX.containsKey(codon.toUpperCase());
    }

    /** Number of nucleotide positions (0-3) at which two same-length codons differ. */
    public static int ntDifferences(String a, String b) {
        int diff = 0;
        for (int i = 0; i < 3; i++) {
            if (a.charAt(i) != b.charAt(i)) diff++;
        }
        return diff;
    }

    /** The 0-based position of the single difference between {@code a} and {@code b}; only valid when {@link #ntDifferences} == 1. */
    public static int differingPosition(String a, String b) {
        for (int i = 0; i < 3; i++) {
            if (a.charAt(i) != b.charAt(i)) return i;
        }
        return -1;
    }

    /** A<->G and C<->T are transitions; every other single-nucleotide substitution is a transversion. */
    public static boolean isTransition(char from, char to) {
        return (from == 'A' && to == 'G') || (from == 'G' && to == 'A')
                || (from == 'C' && to == 'T') || (from == 'T' && to == 'C');
    }

    /** Same encoded amino acid (both codons are sense codons by construction, so translate() never sees a stop here). */
    public static boolean isSynonymous(String a, String b) {
        return GeneticCode.translate(a) == GeneticCode.translate(b);
    }
}
