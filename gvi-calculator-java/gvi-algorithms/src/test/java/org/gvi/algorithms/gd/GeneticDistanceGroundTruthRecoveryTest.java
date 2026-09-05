package org.gvi.algorithms.gd;

import org.gvi.core.model.NucleotideSequence;
import org.junit.jupiter.api.Test;

import java.util.Random;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.within;

/**
 * Validates GD recovery against sequences evolved under a KNOWN true
 * distance -- not a hand-picked pair of sequences chosen to produce a
 * particular p-distance, but per-site stochastic mutation whose expected
 * proportion of differences is set directly, mirroring
 * {@code GroundTruthRecoveryTest} (mu)'s JC69 simulation and
 * {@code DnDsGroundTruthRecoveryTest}'s per-codon acceptance simulation.
 */
class GeneticDistanceGroundTruthRecoveryTest {

    private static final long SEED = 20260813L;
    private static final char[] BASES = {'A', 'C', 'G', 'T'};
    private static final int LENGTH = 50000;

    @Test
    void jukesCantorRecoversAKnownUnbiasedDistance() {
        // No ts/tv bias: each differing site mutates uniformly to one of the other 3 bases,
        // the exact process JC69's correction formula assumes.
        double pTarget = 0.15;
        Random rng = new Random(SEED);
        char[] ancestral = randomSequence(LENGTH, rng);
        char[] query = ancestral.clone();
        for (int i = 0; i < LENGTH; i++) {
            if (rng.nextDouble() < pTarget) {
                char alt;
                do {
                    alt = BASES[rng.nextInt(4)];
                } while (alt == ancestral[i]);
                query[i] = alt;
            }
        }

        double dTrue = -0.75 * Math.log(1.0 - (4.0 / 3.0) * pTarget);
        NucleotideSequence ref = new NucleotideSequence("ref", new String(ancestral));
        NucleotideSequence qry = new NucleotideSequence("qry", new String(query));

        GdResult result = new GeneticDistanceCalculator().compute(ref, qry, GdMethod.JUKES_CANTOR);
        System.out.println("[GeneticDistanceGroundTruthRecoveryTest] JC69: true d=" + dTrue + " recovered=" + result.distance());

        // With this fixed seed the recovered distance actually lands within ~0.0005 of true d.
        assertThat(result.saturated()).isFalse();
        assertThat(result.distance()).isCloseTo(dTrue, within(0.003));
    }

    /**
     * Simulates real transition/transversion bias (kappa=4, transitions 4x
     * more likely per pathway than transversions -- the same kappa
     * convention this codebase's HKY85/K80 models use elsewhere) and checks
     * two things: K80 recovers the true (bias-corrected) distance, AND K80
     * lands measurably closer to the truth than naively applying JC69 (which
     * ignores the ts/tv split) to the same data -- direct evidence that
     * K80's extra correction is actually doing something, not just
     * decoration.
     */
    @Test
    void kimura2ParameterRecoversAKnownTransitionBiasedDistanceMoreAccuratelyThanJc69() {
        double kappaTrue = 4.0;
        double pTarget = 0.15; // total proportion different
        double transitionTarget = pTarget * kappaTrue / (kappaTrue + 2); // 1 transition pathway weighted by kappa
        double transversionTarget = pTarget * 2.0 / (kappaTrue + 2);      // 2 transversion pathways, weight 1 each

        Random rng = new Random(SEED);
        char[] ancestral = randomSequence(LENGTH, rng);
        char[] query = ancestral.clone();
        for (int i = 0; i < LENGTH; i++) {
            double u = rng.nextDouble();
            if (u < transitionTarget) {
                query[i] = transitionPartner(ancestral[i]);
            } else if (u < transitionTarget + transversionTarget) {
                char[] options = transversionPartners(ancestral[i]);
                query[i] = options[rng.nextInt(2)];
            }
        }

        double dTrue = -0.5 * Math.log(1.0 - 2 * transitionTarget - transversionTarget)
                - 0.25 * Math.log(1.0 - 2 * transversionTarget);
        double dJcNaive = -0.75 * Math.log(1.0 - (4.0 / 3.0) * pTarget);

        NucleotideSequence ref = new NucleotideSequence("ref", new String(ancestral));
        NucleotideSequence qry = new NucleotideSequence("qry", new String(query));

        GdResult k80 = new GeneticDistanceCalculator().compute(ref, qry, GdMethod.KIMURA_2_PARAMETER);
        GdResult jc = new GeneticDistanceCalculator().compute(ref, qry, GdMethod.JUKES_CANTOR);
        System.out.println("[GeneticDistanceGroundTruthRecoveryTest] K80: true d=" + dTrue + " k80 recovered=" + k80.distance()
                + " jc69-naive-on-same-data=" + jc.distance() + " (formula-only jc69 estimate=" + dJcNaive + ")");

        // With this fixed seed K80's recovered distance actually lands within ~0.001 of true d.
        assertThat(k80.saturated()).isFalse();
        assertThat(k80.distance()).isCloseTo(dTrue, within(0.003));
        // K80 must land closer to the true distance than plain JC69 does on the SAME transition-biased data --
        // otherwise the extra ts/tv-aware correction isn't buying anything.
        assertThat(Math.abs(k80.distance() - dTrue)).isLessThan(Math.abs(jc.distance() - dTrue));
    }

    private static char[] randomSequence(int length, Random rng) {
        char[] seq = new char[length];
        for (int i = 0; i < length; i++) seq[i] = BASES[rng.nextInt(4)];
        return seq;
    }

    private static char transitionPartner(char base) {
        return switch (base) {
            case 'A' -> 'G';
            case 'G' -> 'A';
            case 'C' -> 'T';
            case 'T' -> 'C';
            default -> throw new IllegalArgumentException("Not a base: " + base);
        };
    }

    private static char[] transversionPartners(char base) {
        return switch (base) {
            case 'A', 'G' -> new char[]{'C', 'T'};
            case 'C', 'T' -> new char[]{'A', 'G'};
            default -> throw new IllegalArgumentException("Not a base: " + base);
        };
    }
}
