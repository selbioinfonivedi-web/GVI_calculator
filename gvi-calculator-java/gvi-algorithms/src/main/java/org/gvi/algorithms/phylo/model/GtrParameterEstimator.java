package org.gvi.algorithms.phylo.model;

import org.gvi.core.exception.GviComputationException;
import org.gvi.core.model.NucleotideSequence;
import org.gvi.core.model.SequenceAlignment;
import org.gvi.core.util.MemoryGuard;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * Empirical (method-of-moments) estimate of GTR parameters from an
 * alignment: observed base composition for the frequencies, and observed
 * pairwise substitution-type counts -- corrected for base-composition bias
 * (dividing by the expected co-occurrence rate pi_i*pi_j under
 * independence) -- for the 6 exchangeability rates. This is a real,
 * standard technique (the same idea behind empirically-derived
 * substitution matrices generally), and a materially better starting
 * point than assuming JC69's "every rate equal" when the data demonstrably
 * isn't.
 * <p>
 * It is <b>not</b> full maximum-likelihood joint estimation of GTR
 * parameters (which would optimize rates and frequencies together against
 * the tree likelihood, iterating with branch lengths/topology) -- that is
 * explicitly out of scope here; see the package documentation for why.
 * This estimate is intended as a solid, honestly-labeled starting point,
 * good enough to materially outperform the single-rate JC69/K80 models
 * used elsewhere in this codebase for distance calculations.
 */
public final class GtrParameterEstimator {

    private static final double PSEUDOCOUNT = 0.5;
    private static final long MAX_PAIRS = 200_000L;

    public record Estimate(double[] baseFrequencies, double[] exchangeabilityRates, List<String> diagnostics) {
    }

    private GtrParameterEstimator() {
    }

    public static Estimate estimate(SequenceAlignment alignment) {
        return estimate(alignment, 42L);
    }

    public static Estimate estimate(SequenceAlignment alignment, long seed) {
        List<NucleotideSequence> sequences = alignment.getSequences();
        List<String> diagnostics = new ArrayList<>();

        long[] baseCounts = new long[4];
        for (NucleotideSequence seq : sequences) {
            String s = seq.getSequence();
            for (int i = 0; i < s.length(); i++) {
                int idx = Nucleotide.index(s.charAt(i));
                if (idx >= 0) baseCounts[idx]++;
            }
        }
        long totalBases = baseCounts[0] + baseCounts[1] + baseCounts[2] + baseCounts[3];
        if (totalBases == 0) {
            throw new GviComputationException("Cannot estimate GTR parameters: no unambiguous bases found in the alignment");
        }
        // Pseudocount so a base that's genuinely never observed (plausible on a small or compositionally
        // skewed alignment, and hit by this class's own test data) still gets a small positive frequency --
        // GtrModel's eigendecomposition requires strictly positive frequencies (division by pi_i, sqrt(pi_j/pi_i)).
        double[] freq = new double[4];
        for (int i = 0; i < 4; i++) freq[i] = (baseCounts[i] + PSEUDOCOUNT) / (totalBases + 4 * PSEUDOCOUNT);

        int n = sequences.size();
        boolean subsample = MemoryGuard.shouldSubsample(n, alignment.length());
        double[][] pairCounts = new double[4][4];
        long pairsUsed = 0;

        if (!subsample) {
            for (int i = 0; i < n; i++) {
                for (int j = i + 1; j < n; j++) {
                    accumulateSubstitutionCounts(sequences.get(i), sequences.get(j), pairCounts);
                    pairsUsed++;
                }
            }
        } else {
            diagnostics.add("Large sequence set (n=" + n + "); subsampled up to " + MAX_PAIRS + " pairs for GTR rate estimation");
            Random rng = new Random(seed);
            long totalPairs = MemoryGuard.pairCount(n);
            long attempts = Math.min(MAX_PAIRS, totalPairs);
            for (long k = 0; k < attempts; k++) {
                int i = rng.nextInt(n), j = rng.nextInt(n);
                if (i == j) continue;
                accumulateSubstitutionCounts(sequences.get(Math.min(i, j)), sequences.get(Math.max(i, j)), pairCounts);
                pairsUsed++;
            }
        }
        diagnostics.add("Base frequencies and exchangeability rates estimated from " + totalBases
                + " bases across " + pairsUsed + " pairwise comparison(s) (empirical/method-of-moments, not full ML)");

        double[] rates = new double[]{
                rate(pairCounts, freq, Nucleotide.A, Nucleotide.C),
                rate(pairCounts, freq, Nucleotide.A, Nucleotide.G),
                rate(pairCounts, freq, Nucleotide.A, Nucleotide.T),
                rate(pairCounts, freq, Nucleotide.C, Nucleotide.G),
                rate(pairCounts, freq, Nucleotide.C, Nucleotide.T),
                rate(pairCounts, freq, Nucleotide.G, Nucleotide.T)
        };

        return new Estimate(freq, rates, diagnostics);
    }

    private static void accumulateSubstitutionCounts(NucleotideSequence a, NucleotideSequence b, double[][] pairCounts) {
        String sa = a.getSequence(), sb = b.getSequence();
        int len = Math.min(sa.length(), sb.length());
        for (int pos = 0; pos < len; pos++) {
            int ia = Nucleotide.index(sa.charAt(pos));
            int ib = Nucleotide.index(sb.charAt(pos));
            if (ia < 0 || ib < 0 || ia == ib) continue;
            pairCounts[ia][ib]++;
            pairCounts[ib][ia]++;
        }
    }

    private static double rate(double[][] pairCounts, double[] freq, int i, int j) {
        double observed = pairCounts[i][j]; // symmetric accumulation, so [i][j] == [j][i]
        double expectedCoOccurrence = freq[i] * freq[j];
        return (observed + PSEUDOCOUNT) / expectedCoOccurrence;
    }
}
