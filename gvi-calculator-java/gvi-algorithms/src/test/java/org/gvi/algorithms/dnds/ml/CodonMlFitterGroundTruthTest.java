package org.gvi.algorithms.dnds.ml;

import org.gvi.algorithms.gd.GdMethod;
import org.gvi.algorithms.phylo.PhyloNode;
import org.gvi.algorithms.phylo.PhyloTree;
import org.gvi.algorithms.phylo.PhyloTreeFactory;
import org.gvi.core.model.NucleotideSequence;
import org.gvi.core.model.SequenceAlignment;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.within;

/**
 * Ground-truth recovery test for {@link CodonMlFitter}/{@link CodonModel},
 * same discipline as the project's other ground-truth tests (mu, Re, GD,
 * RI, Nei-Gojobori dN/dS, BDSKY-ML): simulate codon evolution under a KNOWN
 * kappa/omega via an INDEPENDENT Gillespie (Doob-Gillespie SSA) simulator
 * that reimplements the GY94 rate formula directly from first principles
 * (never calling {@link CodonModel}, {@link CodonTreeLikelihood}, or
 * {@link CodonMlFitter} itself), then checks the ML fit -- started from
 * deliberately wrong kappa/omega -- recovers values close to the true
 * generating parameters. This is the non-circular oracle that would have
 * caught a wrong sign/normalization in the rate matrix the way the BDSKY
 * Monte Carlo checks did for the birth-death model.
 * <p>
 * {@code siteCount=1500} (up from an original 300): this is genuine
 * finite-sample ML variance, not a bug (see {@link CodonMlFitter}'s own
 * commit history/the class javadoc on {@link CodonTreeLikelihood} for the
 * distinct performance bug that WAS real and got fixed separately) -- more
 * simulated codons is more information for the same joint kappa/omega/
 * branch-length fit, so error should and does keep shrinking with N here
 * (unlike the Nei-Gojobori purifying-selection ground-truth test elsewhere
 * in this project, which plateaus at a real residual bias no amount of
 * data removes). At 1500 codons/5 taxa this test runs in ~35s/scenario --
 * scaling site count is cheap (linear cost); scaling taxon count is not
 * (more branches to jointly optimize), which is why taxa were left at 5.
 */
class CodonMlFitterGroundTruthTest {

    @Test
    void recoversKnownKappaAndOmegaUnderPurifyingSelectionFromAGillespieSimulation() {
        double trueKappa = 3.0;
        double trueOmega = 0.25;
        int siteCount = 1500;
        long seed = 20260815L;

        runRecoveryCheck(trueKappa, trueOmega, siteCount, seed, 0.5, 0.06);
    }

    @Test
    void recoversKnownKappaAndOmegaUnderPositiveSelectionFromAGillespieSimulation() {
        double trueKappa = 2.0;
        double trueOmega = 2.5;
        int siteCount = 1500;
        long seed = 20260815L + 1;

        runRecoveryCheck(trueKappa, trueOmega, siteCount, seed, 0.3, 0.2);
    }

    private void runRecoveryCheck(double trueKappa, double trueOmega, int siteCount, long seed,
                                   double kappaTolerance, double omegaTolerance) {
        Random random = new Random(seed);
        double[] pi = buildRealisticCodonFrequencies(random);
        NeighborTable table = NeighborTable.build(pi, trueKappa, trueOmega);

        PhyloTree tree = buildFiveTaxonTreeWithKnownBranchLengths();
        Map<String, List<String>> codonsByLabel = simulate(tree, siteCount, pi, table, random);

        double logLAtTrueParamsAndTrueBranches = new CodonTreeLikelihood()
                .logLikelihood(tree, codonsByLabel, siteCount, CodonModel.of(pi, trueKappa, trueOmega));

        // Deliberately wrong starting point, distinct from both true values, so this tests real recovery
        // (not just "didn't move from a lucky start").
        CodonMlFitter.Result result = new CodonMlFitter().optimize(tree, codonsByLabel, siteCount, pi, 1.0, 1.0);

        System.out.println("[CodonMlFitterGroundTruthTest] true kappa=" + trueKappa + " omega=" + trueOmega
                + " -> recovered kappa=" + result.kappa() + " omega=" + result.omega()
                + " logL=" + result.logLikelihood() + " cycles=" + result.cyclesUsed()
                + " | logL at TRUE params+TRUE branch lengths=" + logLAtTrueParamsAndTrueBranches);

        assertThat(result.logLikelihood()).isGreaterThanOrEqualTo(result.initialLogLikelihood());
        assertThat(result.kappa()).isCloseTo(trueKappa, within(kappaTolerance));
        assertThat(result.omega()).isCloseTo(trueOmega, within(omegaTolerance));
    }

    /** F3x4 frequencies from a synthetic, deliberately non-uniform position-wise nucleotide composition (GC-skewed at position 3, a realistic pattern). */
    private double[] buildRealisticCodonFrequencies(Random random) {
        double[][] positionFreq = {
                {0.30, 0.20, 0.30, 0.20}, // position 1: A,C,G,T
                {0.25, 0.25, 0.25, 0.25}, // position 2
                {0.15, 0.35, 0.15, 0.35}, // position 3
        };
        List<String> sample = new ArrayList<>();
        char[] bases = {'A', 'C', 'G', 'T'};
        while (sample.size() < 5000) {
            StringBuilder codon = new StringBuilder(3);
            for (int pos = 0; pos < 3; pos++) {
                codon.append(bases[weightedDraw(positionFreq[pos], random)]);
            }
            String c = codon.toString();
            if (CodonAlphabet.isSenseCodon(c)) sample.add(c);
        }
        return CodonFrequencies.f3x4(sample);
    }

    private int weightedDraw(double[] weights, Random random) {
        double r = random.nextDouble();
        double cum = 0;
        for (int i = 0; i < weights.length; i++) {
            cum += weights[i];
            if (r <= cum) return i;
        }
        return weights.length - 1;
    }

    /** Topology from a real NJ reconstruction (irrelevant which, since branch lengths are overwritten next), rooted at "A". */
    private PhyloTree buildFiveTaxonTreeWithKnownBranchLengths() {
        int length = 5000;
        char[] base = new char[length];
        java.util.Arrays.fill(base, 'A');

        NucleotideSequence a = mutate("A", base, 0, 0);
        NucleotideSequence b = mutate("B", base, 0, 50);
        NucleotideSequence c = mutate("C", base, 50, 150);
        NucleotideSequence d = mutate("D", base, 150, 350);
        NucleotideSequence e = mutate("E", base, 350, 700);

        SequenceAlignment aln = SequenceAlignment.of(List.of(a, b, c, d, e), "A");
        PhyloTree tree = PhyloTreeFactory.build(aln, GdMethod.JUKES_CANTOR).tree();

        double[] knownLengths = {0.08, 0.05, 0.12, 0.03, 0.10, 0.07, 0.09};
        int[] counter = {0};
        assignBranchLengths(tree.root(), knownLengths, counter);
        return tree;
    }

    private void assignBranchLengths(PhyloNode node, double[] values, int[] counter) {
        for (PhyloNode child : node.children()) {
            child.setBranchLength(values[counter[0]++ % values.length]);
            assignBranchLengths(child, values, counter);
        }
    }

    private NucleotideSequence mutate(String id, char[] base, int fromInclusive, int toExclusive) {
        char[] copy = base.clone();
        for (int i = fromInclusive; i < toExclusive; i++) copy[i] = 'G';
        return new NucleotideSequence(id, new String(copy));
    }

    // ---- independent Gillespie (Doob-Gillespie SSA) codon simulator, reimplementing the GY94 rate formula
    //      directly rather than calling any production class -- the non-circular oracle. ----

    private Map<String, List<String>> simulate(PhyloTree tree, int siteCount, double[] pi, NeighborTable table, Random random) {
        Map<String, List<String>> result = new HashMap<>();
        PhyloNode root = tree.root();
        PhyloNode pruningRoot = root.children().get(0);

        int[] ancestral = new int[siteCount];
        for (int s = 0; s < siteCount; s++) ancestral[s] = weightedDraw(pi, random);

        result.put(root.label(), evolveAll(ancestral, pruningRoot.branchLength(), table, random));
        for (PhyloNode child : pruningRoot.children()) {
            simulateSubtree(child, ancestral, table, random, result);
        }
        return result;
    }

    private void simulateSubtree(PhyloNode node, int[] parentCodons, NeighborTable table, Random random, Map<String, List<String>> result) {
        int[] evolved = new int[parentCodons.length];
        for (int s = 0; s < parentCodons.length; s++) {
            evolved[s] = evolveSite(parentCodons[s], node.branchLength(), table, random);
        }
        if (node.isTaxon()) {
            result.put(node.label(), toCodonStrings(evolved));
        }
        for (PhyloNode child : node.children()) {
            simulateSubtree(child, evolved, table, random, result);
        }
    }

    private List<String> evolveAll(int[] ancestral, double branchLength, NeighborTable table, Random random) {
        int[] evolved = new int[ancestral.length];
        for (int s = 0; s < ancestral.length; s++) evolved[s] = evolveSite(ancestral[s], branchLength, table, random);
        return toCodonStrings(evolved);
    }

    private int evolveSite(int startIndex, double branchLength, NeighborTable table, Random random) {
        int current = startIndex;
        double t = 0;
        while (true) {
            double rate = table.rowSum[current] / table.meanRate;
            if (rate <= 0) return current;
            double waitTime = -Math.log(1 - random.nextDouble()) / rate;
            if (t + waitTime > branchLength) return current;
            t += waitTime;
            current = pickNext(current, table, random);
        }
    }

    private int pickNext(int current, NeighborTable table, Random random) {
        double r = random.nextDouble() * table.rowSum[current];
        double cum = 0;
        int[] neighbors = table.neighborIndex[current];
        double[] rates = table.neighborRate[current];
        for (int k = 0; k < neighbors.length; k++) {
            cum += rates[k];
            if (r <= cum) return neighbors[k];
        }
        return neighbors[neighbors.length - 1];
    }

    private List<String> toCodonStrings(int[] indices) {
        List<String> out = new ArrayList<>(indices.length);
        for (int idx : indices) out.add(CodonAlphabet.SENSE_CODONS.get(idx));
        return out;
    }

    /** Precomputed single-nucleotide-difference neighbor rates for every codon, under a fixed (pi, kappa, omega) -- independently reimplements {@link CodonModel}'s rate formula. */
    private static final class NeighborTable {
        final int[][] neighborIndex;
        final double[][] neighborRate;
        final double[] rowSum;
        final double meanRate;

        private NeighborTable(int[][] neighborIndex, double[][] neighborRate, double[] rowSum, double meanRate) {
            this.neighborIndex = neighborIndex;
            this.neighborRate = neighborRate;
            this.rowSum = rowSum;
            this.meanRate = meanRate;
        }

        static NeighborTable build(double[] pi, double kappa, double omega) {
            int n = CodonAlphabet.SIZE;
            int[][] idx = new int[n][];
            double[][] rate = new double[n][];
            double[] rowSum = new double[n];

            for (int i = 0; i < n; i++) {
                String ci = CodonAlphabet.SENSE_CODONS.get(i);
                List<Integer> js = new ArrayList<>();
                List<Double> rs = new ArrayList<>();
                for (int j = 0; j < n; j++) {
                    if (i == j) continue;
                    String cj = CodonAlphabet.SENSE_CODONS.get(j);
                    if (CodonAlphabet.ntDifferences(ci, cj) != 1) continue;
                    int pos = CodonAlphabet.differingPosition(ci, cj);
                    double r = pi[j];
                    if (CodonAlphabet.isTransition(ci.charAt(pos), cj.charAt(pos))) r *= kappa;
                    if (!CodonAlphabet.isSynonymous(ci, cj)) r *= omega;
                    js.add(j);
                    rs.add(r);
                }
                idx[i] = js.stream().mapToInt(Integer::intValue).toArray();
                rate[i] = rs.stream().mapToDouble(Double::doubleValue).toArray();
                double sum = 0;
                for (double v : rate[i]) sum += v;
                rowSum[i] = sum;
            }

            double meanRate = 0;
            for (int i = 0; i < n; i++) meanRate += pi[i] * rowSum[i];

            return new NeighborTable(idx, rate, rowSum, meanRate);
        }
    }
}
