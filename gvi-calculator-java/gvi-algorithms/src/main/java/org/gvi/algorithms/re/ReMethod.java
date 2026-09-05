package org.gvi.algorithms.re;

/** Which estimator produced a {@link ReResult}, per Section 5.2's data-availability rule. */
public enum ReMethod {
    CORI_INCIDENCE,
    PHYLODYNAMIC_FALLBACK,
    /** Native birth-death-sampling maximum-likelihood fit (the same generative model BEAST2's BDSKY package uses) -- see {@code org.gvi.algorithms.re.bdsky.BdskyReEstimator}. */
    BDSKY_ML
}
