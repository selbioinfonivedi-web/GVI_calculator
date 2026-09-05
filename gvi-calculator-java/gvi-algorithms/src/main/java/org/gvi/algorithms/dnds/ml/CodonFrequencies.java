package org.gvi.algorithms.dnds.ml;

/**
 * F3x4 codon frequency estimation -- codeml's own default (CodonFreq=2):
 * empirical nucleotide frequency is tallied separately at each of the 3
 * codon positions across every observed codon, and a codon's frequency is
 * the product of its 3 position-specific nucleotide frequencies, renormalized
 * over just the 61 sense codons (the raw per-position product also implies
 * some probability mass on the 3 stop-codon combinations, which has to be
 * dropped and the remainder rescaled to sum to 1).
 */
final class CodonFrequencies {

    private CodonFrequencies() {
    }

    static double[] f3x4(Iterable<String> observedSenseCodons) {
        double[][] positionCounts = new double[3][4]; // [position][A,C,G,T]
        double total = 0;
        for (String codon : observedSenseCodons) {
            for (int pos = 0; pos < 3; pos++) {
                int base = baseIndex(codon.charAt(pos));
                if (base < 0) continue;
                positionCounts[pos][base]++;
            }
            total++;
        }
        if (total == 0) {
            throw new IllegalArgumentException("Cannot estimate F3x4 codon frequencies: no observed sense codons");
        }

        double[][] positionFreq = new double[3][4];
        for (int pos = 0; pos < 3; pos++) {
            for (int b = 0; b < 4; b++) {
                positionFreq[pos][b] = positionCounts[pos][b] / total;
            }
        }

        double[] pi = new double[CodonAlphabet.SIZE];
        double sum = 0;
        for (int i = 0; i < CodonAlphabet.SIZE; i++) {
            String codon = CodonAlphabet.SENSE_CODONS.get(i);
            double p = 1.0;
            for (int pos = 0; pos < 3; pos++) {
                p *= positionFreq[pos][baseIndex(codon.charAt(pos))];
            }
            // A codon whose position-wise nucleotides never co-occurred in the data gets exactly 0 mass,
            // which would make the model treat it as literally unreachable (and blow up the symmetric
            // eigendecomposition's sqrt(pi) term) -- float it just above zero instead, same convention
            // as CodonUsageTable's zero-frequency floor for CAI.
            pi[i] = Math.max(p, 1e-8);
            sum += pi[i];
        }
        for (int i = 0; i < pi.length; i++) pi[i] /= sum;
        return pi;
    }

    private static int baseIndex(char c) {
        return switch (Character.toUpperCase(c)) {
            case 'A' -> 0;
            case 'C' -> 1;
            case 'G' -> 2;
            case 'T', 'U' -> 3;
            default -> -1;
        };
    }
}
