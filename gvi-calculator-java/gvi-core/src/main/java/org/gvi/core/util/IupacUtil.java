package org.gvi.core.util;

import java.util.Set;

/** Helpers for recognizing gaps and IUPAC nucleotide ambiguity codes in aligned sequences. */
public final class IupacUtil {

    private static final Set<Character> UNAMBIGUOUS = Set.of('A', 'C', 'G', 'T');
    private static final Set<Character> AMBIGUITY_CODES = Set.of('N', 'R', 'Y', 'W', 'S', 'K', 'M', 'B', 'D', 'H', 'V');
    private static final Set<Character> GAP_CHARS = Set.of('-', '.', '?');

    private IupacUtil() {
    }

    public static boolean isGap(char c) {
        return GAP_CHARS.contains(c);
    }

    public static boolean isAmbiguous(char c) {
        return AMBIGUITY_CODES.contains(Character.toUpperCase(c));
    }

    public static boolean isUnambiguousBase(char c) {
        return UNAMBIGUOUS.contains(Character.toUpperCase(c));
    }

    /** True if this position is usable for pairwise-deletion comparisons (real base, not gap/N/ambiguous). */
    public static boolean isComparable(char c) {
        return isUnambiguousBase(c);
    }

    public static boolean isTransition(char a, char b) {
        a = Character.toUpperCase(a);
        b = Character.toUpperCase(b);
        return (a == 'A' && b == 'G') || (a == 'G' && b == 'A') || (a == 'C' && b == 'T') || (a == 'T' && b == 'C');
    }

    public static boolean isTransversion(char a, char b) {
        return isUnambiguousBase(a) && isUnambiguousBase(b) && Character.toUpperCase(a) != Character.toUpperCase(b) && !isTransition(a, b);
    }
}
