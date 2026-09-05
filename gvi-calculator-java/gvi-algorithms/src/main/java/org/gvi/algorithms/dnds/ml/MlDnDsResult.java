package org.gvi.algorithms.dnds.ml;

/**
 * Result of an ML (GY94/{@code codeml} M0-equivalent) dN/dS fit for one
 * gene -- {@code omega} is the maximum-likelihood dN/dS estimate,
 * {@code kappa} the jointly-fit transition/transversion ratio, alongside
 * the fitted log-likelihood for reference (e.g. comparing across genes or
 * against a null model, should that be added later).
 */
public record MlDnDsResult(String geneName, double omega, double kappa, double logLikelihood,
                            int taxonCount, int codonsUsed, int cyclesUsed) {
}
