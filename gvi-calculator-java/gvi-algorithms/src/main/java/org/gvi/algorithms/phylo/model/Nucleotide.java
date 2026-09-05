package org.gvi.algorithms.phylo.model;

/**
 * The 4-state nucleotide alphabet (A=0, C=1, G=2, T=3) used throughout the
 * GTR+Gamma likelihood machinery, plus IUPAC-ambiguity-aware "partial
 * likelihood" vectors (1.0 for every state the ambiguity code is
 * consistent with, 0.0 otherwise) for Felsenstein pruning at tip nodes.
 */
public final class Nucleotide {

    public static final int A = 0, C = 1, G = 2, T = 3;
    public static final int STATES = 4;

    private Nucleotide() {
    }

    public static int index(char base) {
        return switch (Character.toUpperCase(base)) {
            case 'A' -> A;
            case 'C' -> C;
            case 'G' -> G;
            case 'T', 'U' -> T;
            default -> -1; // ambiguous/gap
        };
    }

    /** 1.0 for every state consistent with this (possibly ambiguous) IUPAC code, else 0.0; all-1 for gaps/N. */
    public static double[] partialLikelihood(char base) {
        double[] p = new double[STATES];
        switch (Character.toUpperCase(base)) {
            case 'A' -> p[A] = 1;
            case 'C' -> p[C] = 1;
            case 'G' -> p[G] = 1;
            case 'T', 'U' -> p[T] = 1;
            case 'R' -> { p[A] = 1; p[G] = 1; }             // puRine
            case 'Y' -> { p[C] = 1; p[T] = 1; }             // pYrimidine
            case 'S' -> { p[C] = 1; p[G] = 1; }             // Strong
            case 'W' -> { p[A] = 1; p[T] = 1; }             // Weak
            case 'K' -> { p[G] = 1; p[T] = 1; }             // Keto
            case 'M' -> { p[A] = 1; p[C] = 1; }             // aMino
            case 'B' -> { p[C] = 1; p[G] = 1; p[T] = 1; }   // not A
            case 'D' -> { p[A] = 1; p[G] = 1; p[T] = 1; }   // not C
            case 'H' -> { p[A] = 1; p[C] = 1; p[T] = 1; }   // not G
            case 'V' -> { p[A] = 1; p[C] = 1; p[G] = 1; }   // not T
            default -> { p[A] = 1; p[C] = 1; p[G] = 1; p[T] = 1; } // N, gap, or anything else: fully ambiguous
        }
        return p;
    }
}
