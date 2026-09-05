package org.gvi.algorithms.re.bdsky;

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.nonstiff.DormandPrince853Integrator;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.within;

/**
 * Self-contained Monte Carlo verification of q1(tau) = P(exactly one
 * sampling event by tau) for a single-lineage birth-death-sampling
 * process, derived via the INHOMOGENEOUS linear ODE
 * {@code dq1/dtau = -(lambda+mu+psi-2*lambda*p)*q1 + psi, q1(0)=0} (same
 * homogeneous part as {@code BirthDeathSamplingPropagator}'s logG, but
 * with an extra "+psi" source term). This derivation is what led to
 * finding and fixing the two real bugs documented in
 * {@link BdskyTreeLikelihood}'s class javadoc -- an earlier, WRONG
 * assumption (that P(exactly one sample by T) equals the simple integral
 * of psi*g(t) from 0 to T) was off by nearly 2x against direct simulation;
 * this correctly-derived inhomogeneous ODE matches simulation almost
 * exactly. Kept as a permanent regression check on the underlying math,
 * independent of {@link BdskyTreeLikelihood} itself.
 */
class Q1VerificationTest {

    @Test
    void q1MatchesDirectSimulation() {
        double lambda = 1.2, mu = 0.3, psi = 0.5, T = 2.0;

        FirstOrderDifferentialEquations ode = new FirstOrderDifferentialEquations() {
            @Override
            public int getDimension() {
                return 2; // y[0]=p, y[1]=q1
            }

            @Override
            public void computeDerivatives(double t, double[] y, double[] yDot) {
                double p = y[0];
                double q1 = y[1];
                yDot[0] = mu - (lambda + mu + psi) * p + lambda * p * p;
                yDot[1] = -(lambda + mu + psi - 2 * lambda * p) * q1 + psi;
            }
        };
        double[] y = {1.0, 0.0};
        DormandPrince853Integrator integrator = new DormandPrince853Integrator(1e-12, T, 1e-12, 1e-12);
        integrator.integrate(ode, 0.0, y, T, y);
        double analyticalQ1 = y[1];

        Random rng = new Random(20260814L);
        int trials = 1_000_000;
        int exactlyOne = 0;
        for (int i = 0; i < trials; i++) {
            if (countSamples(lambda, mu, psi, T, rng) == 1) exactlyOne++;
        }
        double empiricalQ1 = (double) exactlyOne / trials;

        System.out.println("[Q1VerificationTest] q1(T): analytical=" + analyticalQ1 + " empirical=" + empiricalQ1);
        assertThat(analyticalQ1).isCloseTo(empiricalQ1, within(0.01)); // 1M trials -> binomial SE well under 0.01
    }

    private int countSamples(double lambda, double mu, double psi, double T, Random rng) {
        List<Double> active = new ArrayList<>();
        active.add(0.0);
        double time = 0;
        int sampled = 0;
        double totalRatePerLineage = lambda + mu + psi;
        while (!active.isEmpty()) {
            double totalRate = active.size() * totalRatePerLineage;
            double dt = -Math.log(rng.nextDouble()) / totalRate;
            time += dt;
            if (time >= T) break;
            int idx = rng.nextInt(active.size());
            double u = rng.nextDouble() * totalRatePerLineage;
            if (u < lambda) {
                active.add(0.0);
            } else if (u < lambda + mu) {
                active.remove(idx);
            } else {
                active.remove(idx);
                sampled++;
            }
        }
        return sampled;
    }
}
