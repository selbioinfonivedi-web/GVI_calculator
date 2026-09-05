package org.gvi.algorithms.re.bdsky;

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.nonstiff.DormandPrince853Integrator;

/**
 * Numerically integrates the two ODEs behind the birth-death-sampling
 * process tree likelihood (Stadler 2010) -- the same generative model
 * BEAST2's BDSKY package uses for Re -- derived here directly from the
 * underlying continuous-time Markov process via a small-time-interval
 * expansion (not recalled from a memorized closed-form algebraic
 * solution, which has several sign/parameterization conventions across
 * the literature and was judged too risky to reproduce from memory alone;
 * see the project README's "Validation" section for the reasoning).
 * <p>
 * Consider a single lineage, lambda = birth (transmission) rate, mu =
 * death/become-uninfectious-without-sampling rate, psi = sampling (with
 * removal) rate. Over an infinitesimal interval dt, exactly one of four
 * things happens: birth (prob lambda*dt), death (prob mu*dt), sampling
 * (prob psi*dt), or nothing (prob 1-(lambda+mu+psi)*dt).
 * <p>
 * <b>p(tau)</b> = probability a lineage alive tau time units before the
 * present leaves NO sampled descendants by the present. Conditioning on
 * the first event: birth -> both resulting lineages must independently
 * leave no sampled descendants, contributing p(tau-dt)^2; death -> no
 * descendants at all, automatically satisfies "no sampled descendants",
 * contributing 1; sampling -> this lineage IS now sampled, contributing 0.
 * Expanding p(tau) = (1-(lambda+mu+psi)dt)p(tau-dt) + lambda*dt*p(tau-dt)^2
 * + mu*dt and taking dt-&gt;0 gives the Riccati ODE:
 * <pre>  dp/dtau = mu - (lambda+mu+psi)*p + lambda*p^2,  p(0) = 1</pre>
 * (p(0)=1: this project's data is serially sampled through time, not a
 * single-timepoint sample, so there is no separate present-day
 * rho-sampling pulse to condition on).
 * <p>
 * Sanity checks used to verify this ODE (see {@code BirthDeathSamplingPropagatorTest}):
 * if psi=0 (sampling never happens), p=1 is a fixed point of the ODE for
 * any lambda/mu (dp/dtau=0 at p=1 exactly when psi=0) -- consistent with
 * "can never leave a sampled descendant if sampling never occurs". If
 * lambda=mu=0 (pure sampling, no birth/death), the ODE reduces to
 * dp/dtau=-psi*p, i.e. p(tau)=e^(-psi*tau) exactly -- the textbook Poisson
 * survival probability for "not yet sampled by tau", independently
 * checkable in closed form.
 * <p>
 * <b>logG(tau)</b> = log of the branch-likelihood propagator: the
 * (log) density that a specific observed lineage survives tau time units
 * without producing any OTHER sampled descendants beyond what is already
 * accounted for elsewhere in the tree. Conditioning similarly: with rate
 * 2*lambda*p(tau)*dt a "hidden" (unsampled) side-birth occurs and is
 * consistent with continuing to observe only this one lineage (factor 2
 * because either of the two resulting lineages could be "the" continuing
 * one, and the other must itself leave no sampled descendants, probability
 * p(tau)); with rate (lambda+mu+psi)*dt something happens that is
 * inconsistent with this exact continuation. This gives the linear ODE
 * (linear once p(tau) is already known):
 * <pre>  dlogG/dtau = -(lambda+mu+psi - 2*lambda*p(tau)),  logG(0) = 0</pre>
 * Solved in log space throughout (not g(tau) directly) to avoid underflow
 * over long branches / high rates -- the same numerical-stability
 * motivation behind working in log-likelihood space elsewhere in this
 * codebase (e.g. {@code TreeLikelihoodCalculator}).
 */
final class BirthDeathSamplingPropagator {

    private final double lambda, mu, psi;

    BirthDeathSamplingPropagator(double lambda, double mu, double psi) {
        this.lambda = lambda;
        this.mu = mu;
        this.psi = psi;
    }

    record State(double p, double logG) {
    }

    State integrate(double tau) {
        if (tau <= 0) return new State(1.0, 0.0);

        FirstOrderDifferentialEquations ode = new FirstOrderDifferentialEquations() {
            @Override
            public int getDimension() {
                return 2;
            }

            @Override
            public void computeDerivatives(double t, double[] y, double[] yDot) {
                double p = y[0];
                yDot[0] = mu - (lambda + mu + psi) * p + lambda * p * p;
                yDot[1] = -(lambda + mu + psi - 2 * lambda * p);
            }
        };

        double[] y = {1.0, 0.0};
        double minStep = Math.min(1e-10, tau * 1e-8);
        DormandPrince853Integrator integrator = new DormandPrince853Integrator(minStep, tau, 1e-12, 1e-12);
        integrator.integrate(ode, 0.0, y, tau, y);
        return new State(y[0], y[1]);
    }
}
