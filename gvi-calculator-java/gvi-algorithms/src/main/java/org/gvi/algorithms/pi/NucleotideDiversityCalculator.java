package org.gvi.algorithms.pi;

import org.gvi.core.exception.GviComputationException;
import org.gvi.core.model.ConfidenceInterval;
import org.gvi.core.model.NucleotideSequence;
import org.gvi.core.model.SequenceAlignment;
import org.gvi.core.util.IupacUtil;
import org.gvi.core.util.MemoryGuard;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * Index 3 - Nucleotide Diversity (pi), Section 5.3.
 * <p>
 * pi = average, over all pairs of sequences, of (pairwise differences /
 * comparable sites) -- equivalent to the spec's
 * "(sum of pairwise differences) / (pairs x alignment length)" when there
 * are no gaps/ambiguous sites, and generalizes correctly (via pairwise
 * deletion) when there are.
 * <p>
 * Guards (Section 2 / 5.3): for large sequence sets, falls back to random
 * pair subsampling instead of full O(n^2) enumeration. 95% CI is obtained
 * by a nonparametric bootstrap over alignment columns (sites), the standard
 * approach for this statistic, using a seeded RNG for reproducibility.
 */
public final class NucleotideDiversityCalculator {

    private static final long DEFAULT_MAX_PAIRS = 100_000L;

    private final long seed;
    private final int bootstrapReplicates;
    private final long maxPairs;

    public NucleotideDiversityCalculator() {
        this(42L, 1000, DEFAULT_MAX_PAIRS);
    }

    public NucleotideDiversityCalculator(long seed, int bootstrapReplicates, long maxPairs) {
        this.seed = seed;
        this.bootstrapReplicates = bootstrapReplicates;
        this.maxPairs = maxPairs;
    }

    public PiResult compute(SequenceAlignment alignment) {
        int n = alignment.size();
        if (n < 2) {
            throw new GviComputationException("Nucleotide diversity requires at least 2 sequences, got " + n);
        }
        List<String> diagnostics = new ArrayList<>();
        List<NucleotideSequence> sequences = alignment.getSequences();
        int length = alignment.length();

        boolean subsample = MemoryGuard.shouldSubsample(n, length);
        Long pairCap = subsample ? Math.min(maxPairs, MemoryGuard.pairCount(n)) : null;
        if (subsample) {
            diagnostics.add("Sequence set is large (n=" + n + "); using random subsampling of up to "
                    + pairCap + " pairs instead of the full " + MemoryGuard.pairCount(n) + " pairwise comparisons. Result is an estimate.");
        }

        Random rng = new Random(seed);
        PairwiseEstimate main = estimate(sequences, pairCap, rng);
        diagnostics.add(main.pairsWithData + " of " + main.pairsAttempted + " sampled pairs had at least one comparable site");

        ConfidenceInterval ci = bootstrapCi(sequences, length, pairCap, diagnostics);

        var band = PiReferenceTable.classify(main.pi);
        String category = band.scenario() + " (" + band.interpretation() + "); " + band.implication();

        return new PiResult(main.pi, n, length, main.pairsAttempted, subsample, ci, category, diagnostics);
    }

    private ConfidenceInterval bootstrapCi(List<NucleotideSequence> sequences, int length, Long pairCap, List<String> diagnostics) {
        double[] replicatePis = new double[bootstrapReplicates];
        Random rng = new Random(seed + 1);
        int validReplicates = 0;
        for (int r = 0; r < bootstrapReplicates; r++) {
            List<NucleotideSequence> resampled = resampleColumns(sequences, length, rng);
            try {
                PairwiseEstimate est = estimate(resampled, pairCap, rng);
                replicatePis[validReplicates++] = est.pi;
            } catch (GviComputationException e) {
                // a pathological resample with zero comparable sites anywhere; skip this replicate rather than aborting the CI
            }
        }
        if (validReplicates < bootstrapReplicates / 2) {
            diagnostics.add("Bootstrap CI unreliable: only " + validReplicates + "/" + bootstrapReplicates + " replicates were usable");
        }
        if (validReplicates < 2) {
            return null;
        }
        double[] sorted = java.util.Arrays.copyOf(replicatePis, validReplicates);
        java.util.Arrays.sort(sorted);
        double lower = percentile(sorted, 0.025);
        double upper = percentile(sorted, 0.975);
        return new ConfidenceInterval(lower, upper, 0.95);
    }

    private static double percentile(double[] sortedValues, double p) {
        if (sortedValues.length == 1) return sortedValues[0];
        double idx = p * (sortedValues.length - 1);
        int lo = (int) Math.floor(idx);
        int hi = (int) Math.ceil(idx);
        if (lo == hi) return sortedValues[lo];
        double frac = idx - lo;
        return sortedValues[lo] * (1 - frac) + sortedValues[hi] * frac;
    }

    private List<NucleotideSequence> resampleColumns(List<NucleotideSequence> sequences, int length, Random rng) {
        int[] columns = new int[length];
        for (int i = 0; i < length; i++) columns[i] = rng.nextInt(length);
        List<NucleotideSequence> resampled = new ArrayList<>(sequences.size());
        for (NucleotideSequence seq : sequences) {
            char[] out = new char[length];
            String s = seq.getSequence();
            for (int i = 0; i < length; i++) out[i] = s.charAt(columns[i]);
            resampled.add(new NucleotideSequence(seq.getId(), new String(out)));
        }
        return resampled;
    }

    private record PairwiseEstimate(double pi, long pairsAttempted, long pairsWithData) {
    }

    private PairwiseEstimate estimate(List<NucleotideSequence> sequences, Long pairCap, Random rng) {
        int n = sequences.size();
        double sumP = 0.0;
        long pairsWithData = 0;
        long pairsAttempted = 0;

        if (pairCap == null) {
            for (int i = 0; i < n; i++) {
                for (int j = i + 1; j < n; j++) {
                    pairsAttempted++;
                    Double p = pairwiseP(sequences.get(i), sequences.get(j));
                    if (p != null) {
                        sumP += p;
                        pairsWithData++;
                    }
                }
            }
        } else {
            long attempts = Math.max(pairCap, 1);
            for (long k = 0; k < attempts; k++) {
                int i = rng.nextInt(n);
                int j = rng.nextInt(n);
                if (i == j) continue;
                pairsAttempted++;
                Double p = pairwiseP(sequences.get(i), sequences.get(j));
                if (p != null) {
                    sumP += p;
                    pairsWithData++;
                }
            }
        }

        if (pairsWithData == 0) {
            throw new GviComputationException("No sequence pair had any comparable (unambiguous) site in common");
        }
        return new PairwiseEstimate(sumP / pairsWithData, pairsAttempted, pairsWithData);
    }

    /** Pairwise-deletion proportion of differing sites, or null if the two sequences share no comparable site. */
    private Double pairwiseP(NucleotideSequence a, NucleotideSequence b) {
        String sa = a.getSequence();
        String sb = b.getSequence();
        int len = Math.min(sa.length(), sb.length());
        long comparable = 0, diffs = 0;
        for (int i = 0; i < len; i++) {
            char x = sa.charAt(i);
            char y = sb.charAt(i);
            if (!IupacUtil.isComparable(x) || !IupacUtil.isComparable(y)) continue;
            comparable++;
            if (x != y) diffs++;
        }
        return comparable == 0 ? null : (double) diffs / comparable;
    }
}
