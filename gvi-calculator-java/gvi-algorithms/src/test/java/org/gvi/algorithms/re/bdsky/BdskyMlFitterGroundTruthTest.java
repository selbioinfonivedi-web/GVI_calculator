package org.gvi.algorithms.re.bdsky;

import org.gvi.core.model.NucleotideSequence;
import org.gvi.core.model.SequenceAlignment;
import org.gvi.algorithms.re.ReResult;
import org.junit.jupiter.api.Test;

import java.time.LocalDate;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Random;

import static org.assertj.core.api.Assertions.assertThat;

/**
 * End-to-end ground-truth validation of the native BDSKY-equivalent ML
 * fitter: simulates an actual birth-death-sampling branching process
 * (Gillespie algorithm) with a KNOWN true Re, prunes it to the
 * "reconstructed tree" a real phylogenetic analysis would recover, evolves
 * real DNA sequences down that reconstructed topology, then runs the
 * sequences through the REAL pipeline a user would actually exercise --
 * {@code PhyloTreeFactory.build} (NJ reconstruction), {@code LeastSquaresDatingEstimator}
 * (dating), {@code BdskyMlFitter} (this class) -- and checks the recovered
 * Re lands close to the value that generated the data.
 * <p>
 * The first version of this test caught two real bugs in
 * {@code BdskyTreeLikelihood} (see its class javadoc): p/logG evaluated on
 * per-edge-local elapsed time instead of the global remaining-time-to-present
 * axis, and the reference tip incorrectly modeled as both sampled AND
 * having continuing descendants. Both are now fixed and independently
 * verified via direct Monte Carlo simulation (not just this end-to-end
 * test) -- see {@code GPropagatorMonteCarloTest} and {@code Q1VerificationTest}.
 */
class BdskyMlFitterGroundTruthTest {

    private static final long SEED = 20260813L;
    private static final char[] BASES = {'A', 'C', 'G', 'T'};

    /**
     * Re is recovered through the real public entry point, and the reported likelihood-ratio
     * interval must contain the value that generated the data.
     * <p>
     * This replaces an assertion that allowed +/-75% -- by far the loosest in the suite, and loose
     * enough to pass while the estimator was badly wrong. It was: a free three-parameter coordinate
     * ascent over (Re, delta, s) slid along a likelihood ridge to wherever its starting point led,
     * reporting Re = 1.06 for data generated at 1.50 while four other starting points produced
     * 1.15, 1.50 and 1.75 -- all at an identical log-likelihood of -31.72258. The number was a
     * property of the initial guess. delta is now held at the reciprocal of the supplied generation
     * time and Re is profiled, so the answer is a curve over Re rather than one point on a ridge.
     */
    @Test
    void theProfileIntervalCoversTheTrueReproductiveNumber() {
        double[][] cases = {{1.2, 0.3, 0.5}, {2.4, 0.3, 0.5}, {0.9, 0.3, 0.5}};
        for (double[] c : cases) {
            double lambdaTrue = c[0], muTrue = c[1], psiTrue = c[2];
            double deltaTrue = muTrue + psiTrue;
            double trueRe = lambdaTrue / deltaTrue;
            double generationTimeDays = 365.25 / deltaTrue;

            SequenceAlignment alignment = simulateAlignment(lambdaTrue, muTrue, psiTrue, SEED);
            ReResult result = new BdskyReEstimator().compute(alignment, generationTimeDays);

            System.out.printf("[BdskyGroundTruth] true Re=%.3f -> Re=%.4f CI=[%.2f, %.2f]%n",
                    trueRe, result.reMean(), result.interval().lower(), result.interval().upper());

            assertThat(result.interval()).as("Re interval for true Re=%s", trueRe).isNotNull();
            assertThat(result.interval().lower())
                    .as("95%% interval must contain the true Re=%s (got [%s, %s])",
                            trueRe, result.interval().lower(), result.interval().upper())
                    .isLessThanOrEqualTo(trueRe);
            assertThat(result.interval().upper()).isGreaterThanOrEqualTo(trueRe);
            assertThat(result.interval().lower()).isLessThan(result.interval().upper());
        }
    }

    /**
     * Re must respond to the epidemic that generated the tree. The point estimates need not be
     * exact at these tree sizes -- the intervals above carry that -- but a threefold difference in
     * the truth must show up as a difference in the estimate, or the index cannot rank pathogens.
     */
    @Test
    void reRanksASlowEpidemicBelowAFastOne() {
        double slowRe = reFor(0.9, 0.3, 0.5);
        double mediumRe = reFor(1.2, 0.3, 0.5);
        double fastRe = reFor(2.4, 0.3, 0.5);

        System.out.printf("[BdskyGroundTruth] ranking: slow(1.125)=%.3f medium(1.5)=%.3f fast(3.0)=%.3f%n",
                slowRe, mediumRe, fastRe);
        assertThat(slowRe).isLessThan(mediumRe);
        assertThat(mediumRe).isLessThan(fastRe);
    }

    private double reFor(double lambda, double mu, double psi) {
        SequenceAlignment alignment = simulateAlignment(lambda, mu, psi, SEED);
        return new BdskyReEstimator().compute(alignment, 365.25 / (mu + psi)).reMean();
    }

    /** Simulates a birth-death-sampling tree and evolves real sequences down its reconstructed topology. */
    private SequenceAlignment simulateAlignment(double lambda, double mu, double psi, long seed) {
        Random rng = new Random(seed);
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
