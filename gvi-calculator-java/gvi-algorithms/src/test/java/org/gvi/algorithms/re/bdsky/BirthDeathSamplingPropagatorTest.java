package org.gvi.algorithms.re.bdsky;

import org.junit.jupiter.api.Test;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.within;

/**
 * Verifies the numerically-integrated p(tau)/logG(tau) ODEs against two
 * independently-derivable closed-form special cases BEFORE anything is
 * built on top of them (the tree likelihood, the ML fitter) -- see
 * {@link BirthDeathSamplingPropagator}'s class javadoc for the derivations.
 */
class BirthDeathSamplingPropagatorTest {

    @Test
    void pIsExactlyOneForAnyTauWhenSamplingRateIsZero() {
        // psi=0 means sampling never happens -- a lineage can never leave a sampled descendant,
        // so p(tau) (prob of leaving NO sampled descendants) must be exactly 1 for every tau, any lambda/mu.
        for (double lambda : new double[]{0.1, 1.0, 5.0}) {
            for (double mu : new double[]{0.0, 0.5, 2.0}) {
                BirthDeathSamplingPropagator propagator = new BirthDeathSamplingPropagator(lambda, mu, 0.0);
                for (double tau : new double[]{0.1, 1.0, 5.0, 20.0}) {
                    BirthDeathSamplingPropagator.State state = propagator.integrate(tau);
                    assertThat(state.p()).as("lambda=%s mu=%s tau=%s", lambda, mu, tau).isCloseTo(1.0, within(1e-8));
                }
            }
        }
    }

    @Test
    void pMatchesExactPoissonSurvivalWhenOnlySamplingIsActive() {
        // lambda=mu=0: p(tau) reduces to dp/dtau=-psi*p, i.e. the exact closed form p(tau)=e^(-psi*tau)
        // (textbook Poisson "not yet sampled by tau" survival probability).
        for (double psi : new double[]{0.1, 0.5, 2.0}) {
            BirthDeathSamplingPropagator propagator = new BirthDeathSamplingPropagator(0.0, 0.0, psi);
            for (double tau : new double[]{0.01, 0.5, 1.0, 3.0, 10.0}) {
                BirthDeathSamplingPropagator.State state = propagator.integrate(tau);
                double expected = Math.exp(-psi * tau);
                assertThat(state.p()).as("psi=%s tau=%s", psi, tau).isCloseTo(expected, within(1e-8));
            }
        }
    }

    @Test
    void logGIsZeroAtTauZero() {
        BirthDeathSamplingPropagator propagator = new BirthDeathSamplingPropagator(1.0, 0.3, 0.2);
        BirthDeathSamplingPropagator.State state = propagator.integrate(0.0);
        assertThat(state.p()).isEqualTo(1.0);
        assertThat(state.logG()).isEqualTo(0.0);
    }

    @Test
    void logGMatchesExactConstantRateWhenSamplingIsZero() {
        // With psi=0, p(tau)=1 for all tau (check above), so dlogG/dtau = -(lambda+mu+0-2*lambda*1) = lambda-mu
        // exactly (constant) -- giving the exact closed form logG(tau) = (lambda-mu)*tau.
        double lambda = 0.8, mu = 0.3;
        BirthDeathSamplingPropagator propagator = new BirthDeathSamplingPropagator(lambda, mu, 0.0);
        for (double tau : new double[]{0.5, 2.0, 10.0}) {
            BirthDeathSamplingPropagator.State state = propagator.integrate(tau);
            assertThat(state.logG()).isCloseTo((lambda - mu) * tau, within(1e-6));
        }
    }
}
