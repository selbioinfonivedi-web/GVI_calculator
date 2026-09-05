package org.gvi.algorithms.ri;

import org.gvi.core.exception.GviComputationException;
import org.gvi.core.exception.GviInputException;
import org.gvi.core.model.NucleotideSequence;
import org.gvi.core.model.SequenceAlignment;
import org.gvi.core.util.IupacUtil;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * Index 7 - Recombination Index (Section 5.7).
 * <p>
 * Primary method: the PHI test (Pairwise Homoplasy Index; Bruen, Bryant
 * &amp; Poss 2006). Rationale: under clonal (mutation-only) evolution,
 * homoplasy/incompatibility between two polymorphic sites is unrelated to
 * how physically close they are on the genome. Recombination breaks that
 * independence -- nearby sites stay co-inherited (rarely incompatible)
 * while distant sites are free to recombine apart (more often
 * incompatible). So: compute the mean pairwise incompatibility (four-gamete
 * test) among site pairs within a physical window, compare it to the
 * global mean incompatibility, and test significance by permuting site
 * order (which destroys the position/incompatibility association while
 * preserving the incompatibility values themselves).
 * <p>
 * Simplification vs. the reference PhiPack implementation: only biallelic
 * informative sites are used (multi-allelic sites are excluded, flagged in
 * diagnostics) so the four-gamete test is unambiguous, and true breakpoint
 * mapping is not attempted -- {@code ri} is a derived 0..1 signal-strength
 * score (relative drop in local vs. global incompatibility), not a literal
 * breakpoint count. The p-value from the permutation test is the
 * statistically rigorous part of the result; {@code ri} is a convenience
 * scalar for composite weighting.
 */
public final class RecombinationIndexCalculator {

    private static final int MAX_INFORMATIVE_SITES = 1500;
    /** Matches {@link RiResult#significant()}'s own threshold -- kept as one definition to avoid the two silently diverging. */
    private static final double SIGNIFICANCE_THRESHOLD = 0.05;

    private final int windowBp;
    private final int permutations;
    private final long seed;

    public RecombinationIndexCalculator() {
        this(100, 1000, 42L);
    }

    public RecombinationIndexCalculator(int windowBp, int permutations, long seed) {
        this.windowBp = windowBp;
        this.permutations = permutations;
        this.seed = seed;
    }

    private record InformativeSite(int position, char[] alleles) {
    }

    public RiResult computeFromAlignment(SequenceAlignment alignment) {
        List<InformativeSite> sites = extractBiallelicInformativeSites(alignment);
        List<String> diagnostics = new ArrayList<>();

        Random rng = new Random(seed);
        if (sites.size() > MAX_INFORMATIVE_SITES) {
            diagnostics.add("Found " + sites.size() + " informative biallelic sites; randomly subsampled to "
                    + MAX_INFORMATIVE_SITES + " to bound the O(k^2) incompatibility scan");
            java.util.Collections.shuffle(sites, rng);
            sites = new ArrayList<>(sites.subList(0, MAX_INFORMATIVE_SITES));
            sites.sort((a, b) -> Integer.compare(a.position(), b.position()));
        }

        int k = sites.size();
        if (k < 4) {
            throw new GviComputationException("PHI test needs at least 4 informative biallelic sites, found " + k
                    + ". This alignment may be too conserved or too small to assess recombination.");
        }

        boolean[][] incompatible = new boolean[k][k];
        for (int i = 0; i < k; i++) {
            for (int j = i + 1; j < k; j++) {
                boolean inc = fourGameteIncompatible(sites.get(i).alleles(), sites.get(j).alleles());
                incompatible[i][j] = inc;
                incompatible[j][i] = inc;
            }
        }

        int[] positions = new int[k];
        for (int i = 0; i < k; i++) positions[i] = sites.get(i).position();

        double globalPhi = meanIncompatibility(incompatible, null, k);
        WindowedPhi observed = windowedPhi(incompatible, positions, windowBp, k);
        if (observed.pairCount() == 0) {
            throw new GviComputationException("No informative site pairs fall within the " + windowBp
                    + "bp window; widen the window or supply more densely sampled informative sites");
        }

        int atOrBelowObserved = 0;
        int[] shuffledPositions = positions.clone();
        for (int r = 0; r < permutations; r++) {
            shuffle(shuffledPositions, rng);
            WindowedPhi perm = windowedPhi(incompatible, shuffledPositions, windowBp, k);
            if (perm.pairCount() > 0 && perm.phi() <= observed.phi() + 1e-12) {
                atOrBelowObserved++;
            }
        }
        double pValue = (double) (atOrBelowObserved + 1) / (permutations + 1); // add-one smoothing, standard for permutation p-values

        // (globalPhi - windowedPhi) / globalPhi is a *relative* drop, so when globalPhi itself is small (few
        // informative sites / little incompatibility overall) the ratio is numerically unstable -- a tiny absolute
        // difference near zero can swing all the way to 0.0 or 1.0 even though the permutation test found nothing
        // statistically distinguishable from noise. Report that raw ratio only when the permutation test actually
        // cleared significance; otherwise ri is 0.0 (no detected signal), matching what RiResult.significant()
        // already says about the same result, instead of contradicting it with a confident-looking band label.
        double rawSignalStrength = globalPhi > 0 ? clamp01((globalPhi - observed.phi()) / globalPhi) : 0.0;
        boolean significant = pValue < SIGNIFICANCE_THRESHOLD;
        double ri = significant ? rawSignalStrength : 0.0;

        diagnostics.add(k + " informative biallelic sites; " + observed.pairCount() + " pairs within " + windowBp
                + "bp window (of " + (k * (k - 1) / 2) + " total pairs); " + permutations + " permutations run");
        diagnostics.add(String.format("windowed PHI=%.4f, global PHI=%.4f, one-sided permutation p=%.4f", observed.phi(), globalPhi, pValue));
        if (!significant && rawSignalStrength > 0.0) {
            diagnostics.add(String.format(
                    "Permutation p=%.4f is not significant at the %.2f threshold -- ri is reported as 0.0 (no detected "
                            + "recombination signal) rather than the raw windowed/global PHI ratio (%.4f), which is not "
                            + "statistically distinguishable from noise here", pValue, SIGNIFICANCE_THRESHOLD, rawSignalStrength));
        }

        var band = RiReferenceTable.classify(ri);
        String category = band.pathogenContext() + " (" + band.interpretation() + "); " + band.controlImplication();

        return new RiResult(RiMethod.PHI_TEST_BREAKPOINT_SCAN, ri, observed.phi(), globalPhi, pValue, k, permutations, category, diagnostics);
    }

    /** Simpler operational definition (Section 7.1): fraction of isolates an upstream tool already classified as recombinant. */
    public RiResult fromIsolateClassification(int recombinantIsolates, int totalIsolates) {
        if (totalIsolates <= 0) {
            throw new GviInputException("totalIsolates must be positive");
        }
        if (recombinantIsolates < 0 || recombinantIsolates > totalIsolates) {
            throw new GviInputException("recombinantIsolates must be between 0 and totalIsolates");
        }
        double ri = (double) recombinantIsolates / totalIsolates;
        var band = RiReferenceTable.classify(ri);
        String category = band.pathogenContext() + " (" + band.interpretation() + "); " + band.controlImplication();
        List<String> diagnostics = List.of(recombinantIsolates + " of " + totalIsolates + " isolates classified as recombinant by upstream tool");
        return new RiResult(RiMethod.ISOLATE_CLASSIFICATION, ri, null, null, null, null, null, category, diagnostics);
    }

    private record WindowedPhi(double phi, int pairCount) {
    }

    private WindowedPhi windowedPhi(boolean[][] incompatible, int[] positions, int window, int k) {
        long incompatiblePairs = 0, totalPairs = 0;
        for (int i = 0; i < k; i++) {
            for (int j = i + 1; j < k; j++) {
                if (Math.abs(positions[i] - positions[j]) <= window) {
                    totalPairs++;
                    if (incompatible[i][j]) incompatiblePairs++;
                }
            }
        }
        return new WindowedPhi(totalPairs == 0 ? 0.0 : (double) incompatiblePairs / totalPairs, (int) totalPairs);
    }

    private double meanIncompatibility(boolean[][] incompatible, int[] unused, int k) {
        long inc = 0, total = 0;
        for (int i = 0; i < k; i++) {
            for (int j = i + 1; j < k; j++) {
                total++;
                if (incompatible[i][j]) inc++;
            }
        }
        return total == 0 ? 0.0 : (double) inc / total;
    }

    private void shuffle(int[] array, Random rng) {
        for (int i = array.length - 1; i > 0; i--) {
            int j = rng.nextInt(i + 1);
            int tmp = array[i];
            array[i] = array[j];
            array[j] = tmp;
        }
    }

    /** Classic four-gamete test: incompatible iff all 4 combinations of the two biallelic sites' states co-occur in some sequence. */
    boolean fourGameteIncompatible(char[] site1, char[] site2) {
        char a1 = 0, a2 = 0, b1 = 0, b2 = 0;
        boolean sawA1B1 = false, sawA1B2 = false, sawA2B1 = false, sawA2B2 = false;
        for (int s = 0; s < site1.length; s++) {
            char x = site1[s];
            char y = site2[s];
            if (x == 0 || y == 0) continue; // missing/excluded at this sequence for one of the two sites
            if (a1 == 0) a1 = x;
            else if (a1 != x && a2 == 0) a2 = x;
            if (b1 == 0) b1 = y;
            else if (b1 != y && b2 == 0) b2 = y;

            if (x == a1 && y == b1) sawA1B1 = true;
            else if (x == a1 && y == b2) sawA1B2 = true;
            else if (x == a2 && y == b1) sawA2B1 = true;
            else if (x == a2 && y == b2) sawA2B2 = true;
        }
        return sawA1B1 && sawA1B2 && sawA2B1 && sawA2B2;
    }

    private List<InformativeSite> extractBiallelicInformativeSites(SequenceAlignment alignment) {
        List<NucleotideSequence> seqs = alignment.getSequences();
        int n = seqs.size();
        int length = alignment.length();
        List<InformativeSite> sites = new ArrayList<>();

        for (int pos = 0; pos < length; pos++) {
            char[] column = new char[n];
            java.util.Map<Character, Integer> counts = new java.util.HashMap<>();
            for (int s = 0; s < n; s++) {
                char c = Character.toUpperCase(seqs.get(s).getSequence().charAt(pos));
                if (!IupacUtil.isUnambiguousBase(c)) {
                    column[s] = 0; // treated as missing at this site
                    continue;
                }
                column[s] = c;
                counts.merge(c, 1, Integer::sum);
            }
            // parsimony-informative + biallelic: exactly 2 distinct states, each appearing >= 2 times
            if (counts.size() != 2) continue;
            boolean bothAtLeastTwo = counts.values().stream().allMatch(v -> v >= 2);
            if (!bothAtLeastTwo) continue;

            sites.add(new InformativeSite(pos, column));
        }
        return sites;
    }

    private static double clamp01(double v) {
        return Math.max(0.0, Math.min(1.0, v));
    }
}
