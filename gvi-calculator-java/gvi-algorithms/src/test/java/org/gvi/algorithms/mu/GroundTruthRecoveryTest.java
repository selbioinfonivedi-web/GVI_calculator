package org.gvi.algorithms.mu;

import org.assertj.core.data.Percentage;
import org.gvi.core.model.NucleotideSequence;
import org.gvi.core.model.SequenceAlignment;
import org.junit.jupiter.api.Test;

import java.time.LocalDate;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import static org.assertj.core.api.Assertions.assertThat;

/**
 * Validates mu recovery against data simulated under a KNOWN Jukes-Cantor
 * molecular clock -- not a hand-constructed additive-distance toy example
 * (see the other fixtures in {@link EvolutionaryRateCalculatorTest}, which
 * are noiseless by design and test arithmetic correctness, not statistical
 * recovery under real mutational noise). This is the standard way
 * phylogenetics software is benchmarked: simulate sequences under a known
 * rate/model, then check how close the estimator's recovered rate lands to
 * the value that actually generated the data.
 * <p>
 * Each of the {@link #TAXA} query taxa is generated independently from a
 * shared ancestral reference sequence by applying the exact JC69 transition
 * probability for its own elapsed time, P_same(t) = 1/4 + 3/4*e^(-4/3*mu*t)
 * (Jukes &amp; Cantor 1969 -- mu*t is itself denominated in expected
 * substitutions/site, the same quantity {@code d} in this codebase's own
 * {@code -3/4*ln(1-4p/3)} JC-correction formula, which is why the exponent's
 * coefficient is -4/3 and not -4) -- the same closed-form process this codebase's
 * own JC-correction distance formula assumes, so recovering the true rate
 * here is a genuine test of the whole pipeline (distance computation, JC
 * correction, root-to-tip regression), not a tautology. This is
 * mathematically equivalent to simulating a real branching tree and
 * measuring each leaf's total elapsed time from the root (JC69's
 * time-homogeneity means composing two JC69 draws of duration t1, t2 along
 * a real lineage is distributionally identical to one draw of duration
 * t1+t2 -- the same additivity property this codebase's Neighbor-Joining
 * already relies on), so it validates the clock-rate recovery claim
 * without needing a full nested-tree simulator. A fixed random seed makes
 * the result exactly reproducible.
 */
class GroundTruthRecoveryTest {

    private static final long SEED = 20260813L;
    private static final double TRUE_MU = 2.0e-3; // substitutions/site/year
    private static final int GENOME_LENGTH = 5000;
    private static final int TAXA = 30;
    private static final double MAX_YEARS = 10.0;
    private static final char[] BASES = {'A', 'C', 'G', 'T'};

    @Test
    void recoversTheTrueSimulatedClockRateWithinStatisticalTolerance() {
        Random rng = new Random(SEED);
        LocalDate baseDate = LocalDate.of(2020, 1, 1);

        char[] ancestral = new char[GENOME_LENGTH];
        for (int i = 0; i < GENOME_LENGTH; i++) ancestral[i] = BASES[rng.nextInt(4)];

        NucleotideSequence reference = new NucleotideSequence("ref", new String(ancestral)).withMetadata(baseDate, null, null);
        List<NucleotideSequence> all = new ArrayList<>();
        all.add(reference);

        for (int t = 0; t < TAXA; t++) {
            double elapsedYears = 0.5 + rng.nextDouble() * (MAX_YEARS - 0.5);
            char[] evolved = simulateJc69(ancestral, TRUE_MU, elapsedYears, rng);
            LocalDate date = baseDate.plusDays(Math.round(elapsedYears * 365.25));
            all.add(new NucleotideSequence("taxon" + t, new String(evolved)).withMetadata(date, null, null));
        }

        SequenceAlignment alignment = SequenceAlignment.of(all, "ref");
        EvolutionaryRateCalculator calc = new EvolutionaryRateCalculator();

        MuResult pairwise = calc.compute(alignment);
        MuResult treeAware = calc.computeTreeAware(alignment);

        System.out.println("[GroundTruthRecoveryTest] true mu=" + TRUE_MU
                + " pairwise mu=" + pairwise.muSubPerSiteYear() + " R^2=" + pairwise.rSquared()
                + " treeAware mu=" + treeAware.muSubPerSiteYear() + " R^2=" + treeAware.rSquared());

        // With this fixed seed both estimators actually land within ~0.5% of the true rate (R^2~=0.92) --
        // asserting a still-tight but not knife-edge 10% keeps real headroom against unrelated future changes
        // (e.g. a distance-calculation refinement) without turning this into a no-op tolerance.
        assertThat(pairwise.muSubPerSiteYear()).isCloseTo(TRUE_MU, Percentage.withPercentage(10));
        assertThat(treeAware.muSubPerSiteYear()).isCloseTo(TRUE_MU, Percentage.withPercentage(10));
        assertThat(pairwise.rSquared()).isGreaterThan(0.85);
        assertThat(treeAware.rSquared()).isGreaterThan(0.85);
    }

    /**
     * Exact JC69 process for elapsed time t: each site independently keeps
     * the ancestral base with probability P_same(t), else mutates uniformly
     * to one of the other 3 bases (Jukes &amp; Cantor 1969's closed-form
     * transition probability, not a discretized approximation). mu*t is
     * denominated in expected substitutions/site (this codebase's own
     * JC-correction formula's {@code d}), so the exponent's coefficient is
     * -4/3, not -4.
     */
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
