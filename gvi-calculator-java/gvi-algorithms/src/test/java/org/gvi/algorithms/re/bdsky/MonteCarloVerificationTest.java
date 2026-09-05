package org.gvi.algorithms.re.bdsky;

import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.within;

/**
 * Independent Monte Carlo verification of p(tau) for GENERIC (non-degenerate)
 * rate combinations, via direct Gillespie simulation of the birth-death-
 * sampling process -- not just the closed-form special cases already
 * checked in {@link BirthDeathSamplingPropagatorTest}. This is the same
 * simulator technique used to investigate the BDSKY-ML ground-truth
 * failure; kept as a permanent regression check once a normalization fix
 * is verified this way.
 */
class MonteCarloVerificationTest {

    @Test
    void pMatchesDirectSimulationForGenericRates() {
        double lambda = 1.2, mu = 0.3, psi = 0.5, tau = 2.0;
        Random rng = new Random(20260814L);
        int trials = 300_000;
        int noSampleCount = 0;
        for (int i = 0; i < trials; i++) {
            if (!anySampleWithinTime(lambda, mu, psi, tau, rng)) noSampleCount++;
        }
        double empiricalP = (double) noSampleCount / trials;
        double analyticalP = new BirthDeathSamplingPropagator(lambda, mu, psi).integrate(tau).p();

        System.out.println("[MonteCarloVerificationTest] p(tau): empirical=" + empiricalP + " analytical=" + analyticalP);
        assertThat(analyticalP).isCloseTo(empiricalP, within(0.01)); // 300k trials -> binomial SE well under 0.01
    }

    /** Simulates a single lineage's birth-death-sampling process for exactly {@code tau} time; returns true if any sampling event occurred. */
    private boolean anySampleWithinTime(double lambda, double mu, double psi, double tau, Random rng) {
        List<Double> active = new ArrayList<>(); // just need a count really, but keep it simple
        active.add(0.0);
        double time = 0;
        double totalRatePerLineage = lambda + mu + psi;
        while (!active.isEmpty()) {
            double totalRate = active.size() * totalRatePerLineage;
            double dt = -Math.log(rng.nextDouble()) / totalRate;
            time += dt;
            if (time >= tau) return false;
            int idx = rng.nextInt(active.size());
            double u = rng.nextDouble() * totalRatePerLineage;
            if (u < lambda) {
                active.add(0.0); // birth: one more active lineage
            } else if (u < lambda + mu) {
                active.remove(idx); // death
            } else {
                return true; // sampled
            }
        }
        return false; // all lineages died out without ever being sampled
    }
}
