package org.gvi.algorithms.re;

import org.gvi.core.model.NucleotideSequence;
import org.gvi.core.model.SequenceAlignment;
import org.junit.jupiter.api.Test;

import java.time.LocalDate;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Random;

import static org.assertj.core.api.Assertions.assertThat;

/**
 * Ground-truth check on the lineages-through-time fallback, which previously had none.
 * <p>
 * It needed one. The Euler-Lotka conversion discretizes the serial interval over [1, tMax] days and
 * {@link SerialInterval} renormalizes the weights across exactly that window, but tMax was hardcoded
 * to 20 days for every pathogen. A Gamma with a 457-day mean evaluated on a 20-day grid does not
 * become a truncated 457-day distribution -- it becomes a 20-day one, and
 * {@code Re = 1 / sum_u w(u) e^(-r u)} collapses to roughly {@code 1 + r*(20/365)}. Measured before
 * the fix, the estimator returned 1.03, 1.07 and 1.03 for true Re of 1.5, 3.0 and 1.125: it could
 * not tell a tripling epidemic from a static one, which is exactly what the whole 16-dataset corpus
 * showed when every pathogen came back at Re ~ 1.0.
 * <p>
 * The growth rate itself was never the problem -- it read 0.535, 0.689 and 1.496/yr for those three
 * cases, correctly ordered. Only the conversion was broken.
 */
class PhylodynamicReGroundTruthTest {

    private static final long SEED = 20260813L;
    private static final char[] BASES = {'A', 'C', 'G', 'T'};

    /**
     * The estimator must separate epidemics that differ threefold in Re. This is the property the
     * 20-day horizon destroyed, and the one the corpus needs: an index that returns ~1.0 for
     * everything cannot rank pathogens, however precise it looks.
     */
    @Test
    void separatesEpidemicsThatDifferInTrueReproductiveNumber() {
        double slow = reFor(0.9, 0.3, 0.5);    // true Re = 1.125
        double medium = reFor(1.2, 0.3, 0.5);  // true Re = 1.5
        double fast = reFor(2.4, 0.3, 0.5);    // true Re = 3.0

        System.out.printf("[PhylodynamicGroundTruth] true 1.125 -> %.3f | true 1.5 -> %.3f | true 3.0 -> %.3f%n",
                slow, medium, fast);

        assertThat(slow).isLessThan(medium);
        assertThat(medium).isLessThan(fast);
        // Before the fix all three sat inside [1.00, 1.08]; a threefold change in truth must move
        // the estimate by substantially more than that band was wide.
        assertThat(fast - slow).as("spread across a threefold difference in true Re").isGreaterThan(0.5);
    }

    /**
     * A faster epidemic must not be reported as static. The pre-fix estimator put true Re = 3.0 at
     * 1.07, which classifies as "endemic, routine vaccination" -- the wrong control measure.
     */
    @Test
    void aFastEpidemicIsNotReportedAsEndemic() {
        double fast = reFor(2.4, 0.3, 0.5); // true Re = 3.0
        assertThat(fast).as("true Re=3.0 must not read as stable circulation").isGreaterThan(1.5);
    }

    /**
     * The discretization horizon must follow the distribution it is discretizing. A short viral
     * generation time keeps the historic 20-day grid; a slow pathogen gets a grid that actually
     * carries the Gamma's mass.
     */
    @Test
    void theSerialIntervalHorizonScalesWithTheGenerationTime() {
        assertThat(PhylodynamicReEstimator.serialIntervalTMax(5.0, 2.0)).isEqualTo(20);
        assertThat(PhylodynamicReEstimator.serialIntervalTMax(457.0, 182.8))
                .as("a 457-day mean cannot be represented on a 20-day grid")
                .isGreaterThan(1000);
    }

    private double reFor(double lambda, double mu, double psi) {
        SequenceAlignment alignment = simulateAlignment(lambda, mu, psi);
        return new PhylodynamicReEstimator().compute(alignment, 365.25 / (mu + psi)).reMean();
    }

    private SequenceAlignment simulateAlignment(double lambda, double mu, double psi) {
        Random rng = new Random(SEED);
        SimNode simRoot = simulateBirthDeathSampling(lambda, mu, psi, 30, 30.0, rng);
        assertThat(hasSampledDescendant(simRoot)).isTrue();
        SimNode reconstructed = collapseToReconstructedRoot(simRoot);

        char[] ancestral = new char[5000];
        for (int i = 0; i < ancestral.length; i++) ancestral[i] = BASES[rng.nextInt(4)];
        List<NucleotideSequence> tips = new ArrayList<>();
        assignSequences(reconstructed, ancestral, reconstructed.birthTime, 2.0e-3, rng, tips);
        assertThat(tips.size()).isGreaterThanOrEqualTo(10);

        NucleotideSequence reference = tips.stream()
                .min(Comparator.comparingDouble(t -> t.getCollectionDate().get().toEpochDay())).orElseThrow();
        return SequenceAlignment.of(tips, reference.getId());
    }

    static final class SimNode {
        double birthTime;
        double endTime = Double.NaN;
        boolean sampled = false;
        SimNode left, right;
    }

    private SimNode simulateBirthDeathSampling(double lambda, double mu, double psi, int targetSamples, double maxTime, Random rng) {
        List<SimNode> active = new ArrayList<>();
        SimNode origin = new SimNode();
        origin.birthTime = 0;
        active.add(origin);

        double time = 0;
        int sampledCount = 0;
        double totalRatePerLineage = lambda + mu + psi;

        while (!active.isEmpty() && sampledCount < targetSamples) {
            double totalRate = active.size() * totalRatePerLineage;
            double dt = -Math.log(rng.nextDouble()) / totalRate;
            time += dt;
            if (time >= maxTime) break;

            int idx = rng.nextInt(active.size());
            SimNode chosen = active.get(idx);
            double u = rng.nextDouble() * totalRatePerLineage;
            if (u < lambda) {
                SimNode left = new SimNode();
                left.birthTime = time;
                SimNode right = new SimNode();
                right.birthTime = time;
                chosen.left = left;
                chosen.right = right;
                chosen.endTime = time;
                active.remove(idx);
                active.add(left);
                active.add(right);
            } else if (u < lambda + mu) {
                chosen.endTime = time;
                chosen.sampled = false;
                active.remove(idx);
            } else {
                chosen.endTime = time;
                chosen.sampled = true;
                active.remove(idx);
                sampledCount++;
            }
        }
        return origin;
    }

    private boolean hasSampledDescendant(SimNode n) {
        if (n.sampled) return true;
        if (n.left != null) return hasSampledDescendant(n.left) || hasSampledDescendant(n.right);
        return false;
    }

    private SimNode collapseToReconstructedRoot(SimNode n) {
        if (n.sampled || n.left == null) return n;
        boolean leftHas = hasSampledDescendant(n.left);
        boolean rightHas = hasSampledDescendant(n.right);
        if (leftHas && rightHas) return n;
        if (leftHas) return collapseToReconstructedRoot(n.left);
        return collapseToReconstructedRoot(n.right);
    }

    private void assignSequences(SimNode node, char[] parentSequence, double parentTime, double clockRate, Random rng, List<NucleotideSequence> tips) {
        double elapsed = node.endTime - parentTime;
        char[] sequence = elapsed > 0 ? simulateJc69(parentSequence, clockRate, elapsed, rng) : parentSequence;

        if (node.sampled) {
            LocalDate date = LocalDate.of(2000, 1, 1).plusDays(Math.round(node.endTime * 365.25));
            String id = "tip" + tips.size() + "_" + date;
            tips.add(new NucleotideSequence(id, new String(sequence)).withMetadata(date, null, null));
            return;
        }

        boolean leftHas = hasSampledDescendant(node.left);
        boolean rightHas = hasSampledDescendant(node.right);
        SimNode effectiveLeft = leftHas ? collapseToReconstructedRoot(node.left) : null;
        SimNode effectiveRight = rightHas ? collapseToReconstructedRoot(node.right) : null;
        if (effectiveLeft != null) assignSequences(effectiveLeft, sequence, node.endTime, clockRate, rng, tips);
        if (effectiveRight != null) assignSequences(effectiveRight, sequence, node.endTime, clockRate, rng, tips);
    }

    private static char[] simulateJc69(char[] ancestral, double mu, double elapsedYears, Random rng) {
        double pSame = 0.25 + 0.75 * Math.exp(-(4.0 / 3.0) * mu * elapsedYears);
        char[] evolved = new char[ancestral.length];
        for (int i = 0; i < ancestral.length; i++) {
            if (rng.nextDouble() < pSame) {
                evolved[i] = ancestral[i];
            } else {
                char newBase;
                do {
                    newBase = BASES[rng.nextInt(4)];
                } while (newBase == ancestral[i]);
                evolved[i] = newBase;
            }
        }
        return evolved;
    }
}
