package org.gvi.algorithms.ri;

import org.gvi.algorithms.gd.GdMethod;
import org.gvi.algorithms.gd.GdResult;
import org.gvi.algorithms.gd.GeneticDistanceCalculator;
import org.gvi.core.exception.GviComputationException;
import org.gvi.core.model.NucleotideSequence;
import org.gvi.core.model.SequenceAlignment;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

/**
 * Detects reassortment between two segments of a segmented-genome pathogen (e.g. Bluetongue's 10
 * segments, influenza's 8, rotavirus's 11) -- a genuinely different question from within-locus
 * recombination ({@link RecombinationIndexCalculator}'s PHI test). PHI only sees one contiguous
 * alignment and cannot represent "this isolate's segment 2 traces one lineage while its segment 10
 * traces a different one," which is what reassortment actually is; testing for it needs two
 * INDEPENDENTLY-aligned segments for the SAME set of isolates (matched by sequence id).
 * <p>
 * Method: a Mantel test (Mantel 1967) -- the standard technique for testing congruence between two
 * distance matrices from different loci, widely used in cophylogenetics to compare gene trees. Under
 * normal joint inheritance (no reassortment), both segments travel down the same transmission chain, so
 * isolates that are close on segment A should generally also be close on segment B -- the two segments'
 * pairwise-distance matrices should be positively correlated. Reassortment breaks that correspondence
 * for whichever isolates picked up a segment from a different lineage, pulling the correlation down.
 * <p>
 * The observed statistic is the Pearson correlation ("Mantel r") between the two matrices' matched-taxon
 * pairwise distances. Its significance is assessed by permutation: repeatedly relabel which taxon goes
 * with which row/column of one matrix (destroying any true taxon-to-taxon correspondence while
 * preserving each matrix's own internal distance distribution) and see how often a fully-randomized
 * pairing produces a correlation at least as high as what was actually observed. A low p-value means the
 * observed concordance is unlikely to be a coincidence -- real shared signal, consistent with joint
 * inheritance. A p-value that isn't significant means the two segments' relationships among these taxa
 * are statistically indistinguishable from randomly-matched ones -- consistent with reassortment having
 * decoupled their histories (though, same caveat as {@link RecombinationIndexCalculator}, it can also
 * simply mean too few/too-similar taxa to detect a real signal either way).
 */
public final class ReassortmentIndexCalculator {

    /** Mantel tests are only meaningful with enough taxa for the permutation null to be informative -- matches the general spirit of {@link RecombinationIndexCalculator}'s own minimum. */
    private static final int MIN_COMMON_TAXA = 4;
    private static final double SIGNIFICANCE_THRESHOLD = 0.05;

    private final int permutations;
    private final long seed;

    public ReassortmentIndexCalculator() {
        this(1000, 42L);
    }

    public ReassortmentIndexCalculator(int permutations, long seed) {
        this.permutations = permutations;
        this.seed = seed;
    }

    public ReassortmentResult compare(String segmentALabel, SequenceAlignment segmentA,
                                       String segmentBLabel, SequenceAlignment segmentB) {
        List<String> commonIds = commonTaxa(segmentA, segmentB);
        if (commonIds.size() < MIN_COMMON_TAXA) {
            throw new GviComputationException("Reassortment test between '" + segmentALabel + "' and '" + segmentBLabel
                    + "' needs at least " + MIN_COMMON_TAXA + " sequence ids shared by both segment alignments (matched "
                    + "by id -- the same isolate must use the same sequence_id in both files), found only "
                    + commonIds.size());
        }

        double[][] distA = distanceMatrix(segmentA, commonIds);
        double[][] distB = distanceMatrix(segmentB, commonIds);
        int n = commonIds.size();

        double observedR = pearsonUpperTriangle(distA, distB, null);

        Random rng = new Random(seed);
        int[] order = new int[n];
        for (int i = 0; i < n; i++) order[i] = i;
        int atOrAboveObserved = 0;
        for (int rep = 0; rep < permutations; rep++) {
            shuffle(order, rng);
            double permutedR = pearsonUpperTriangle(distA, distB, order);
            if (permutedR >= observedR - 1e-12) atOrAboveObserved++;
        }
        double pValue = (double) (atOrAboveObserved + 1) / (permutations + 1); // add-one smoothing, standard for permutation p-values

        double reassortmentIndex = clamp01(1.0 - observedR);
        boolean significant = pValue < SIGNIFICANCE_THRESHOLD;

        List<String> diagnostics = new ArrayList<>();
        diagnostics.add(n + " taxa shared by both segments (matched by sequence id); " + permutations + " permutations run");
        diagnostics.add(String.format(java.util.Locale.ROOT, "Mantel r=%.4f, one-sided permutation p=%.4f", observedR, pValue));

        String category = significant
                ? String.format(java.util.Locale.ROOT,
                        "Concordant segments (Mantel r=%.2f, p=%.3f) -- consistent with joint inheritance; little evidence "
                                + "of reassortment between '%s' and '%s' across these %d taxa",
                        observedR, pValue, segmentALabel, segmentBLabel, n)
                : String.format(java.util.Locale.ROOT,
                        "No significant concordance detected (Mantel r=%.2f, p=%.3f) between '%s' and '%s' across these %d "
                                + "taxa -- consistent with reassortment having decoupled these segments' histories for at "
                                + "least some isolates, though this can also simply reflect too few or too-similar taxa to "
                                + "detect a real signal either way",
                        observedR, pValue, segmentALabel, segmentBLabel, n);

        return new ReassortmentResult(segmentALabel, segmentBLabel, n, observedR, pValue, reassortmentIndex, category, diagnostics);
    }

    private List<String> commonTaxa(SequenceAlignment a, SequenceAlignment b) {
        Set<String> idsB = new LinkedHashSet<>();
        for (NucleotideSequence seq : b.getSequences()) idsB.add(seq.getId());
        List<String> common = new ArrayList<>();
        for (NucleotideSequence seq : a.getSequences()) {
            if (idsB.contains(seq.getId())) common.add(seq.getId());
        }
        return common;
    }

    private double[][] distanceMatrix(SequenceAlignment alignment, List<String> orderedIds) {
        Map<String, NucleotideSequence> byId = new HashMap<>();
        for (NucleotideSequence seq : alignment.getSequences()) byId.put(seq.getId(), seq);
        int n = orderedIds.size();
        double[][] matrix = new double[n][n];
        GeneticDistanceCalculator calc = new GeneticDistanceCalculator();
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                GdResult gd = calc.compute(byId.get(orderedIds.get(i)), byId.get(orderedIds.get(j)), GdMethod.JUKES_CANTOR);
                double d = (gd.saturated() || gd.distance() == null) ? gd.pDistance() : gd.distance();
                matrix[i][j] = d;
                matrix[j][i] = d;
            }
        }
        return matrix;
    }

    /** Pearson correlation between distA's upper triangle and distB's (optionally taxon-permuted) upper triangle. */
    private double pearsonUpperTriangle(double[][] distA, double[][] distB, int[] permutedOrder) {
        int n = distA.length;
        List<Double> xs = new ArrayList<>();
        List<Double> ys = new ArrayList<>();
        for (int i = 0; i < n; i++) {
            int pi = permutedOrder == null ? i : permutedOrder[i];
            for (int j = i + 1; j < n; j++) {
                int pj = permutedOrder == null ? j : permutedOrder[j];
                xs.add(distA[i][j]);
                ys.add(distB[pi][pj]);
            }
        }
        int m = xs.size();
        double meanX = xs.stream().mapToDouble(Double::doubleValue).average().orElse(0);
        double meanY = ys.stream().mapToDouble(Double::doubleValue).average().orElse(0);
        double num = 0, denomX = 0, denomY = 0;
        for (int k = 0; k < m; k++) {
            double dx = xs.get(k) - meanX;
            double dy = ys.get(k) - meanY;
            num += dx * dy;
            denomX += dx * dx;
            denomY += dy * dy;
        }
        if (denomX <= 0 || denomY <= 0) return 0.0;
        return num / Math.sqrt(denomX * denomY);
    }

    private void shuffle(int[] array, Random rng) {
        for (int i = array.length - 1; i > 0; i--) {
            int j = rng.nextInt(i + 1);
            int tmp = array[i];
            array[i] = array[j];
            array[j] = tmp;
        }
    }

    private static double clamp01(double v) {
        return Math.max(0.0, Math.min(1.0, v));
    }
}
