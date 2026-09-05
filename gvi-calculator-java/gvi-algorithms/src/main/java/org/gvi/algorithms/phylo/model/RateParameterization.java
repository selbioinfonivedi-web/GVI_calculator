package org.gvi.algorithms.phylo.model;

/**
 * Maps a small vector of free rate parameters (e.g. a single
 * transition/transversion ratio kappa) onto the 6-element GTR
 * exchangeability rate array (order: AC, AG, AT, CG, CT, GT), and records
 * how many free parameters that costs for AIC/BIC model comparison.
 * <p>
 * This is what makes JC69, F81, K80, HKY85, and TN93 available as named
 * models without separate substitution-model math: they are all special
 * cases of GTR (equal frequencies and/or constrained/tied rates), so each
 * is just a different, smaller parameterization feeding the same
 * {@link GtrModel}. See {@link NamedSubstitutionModels} for the concrete
 * instances.
 */
public record RateParameterization(
        String modelName,
        boolean useEmpiricalFrequencies,
        double[] initialTheta,
        double[] lowerBounds,
        double[] upperBounds,
        java.util.function.Function<double[], double[]> thetaToRates
) {
    public int freeRateParameterCount() {
        return initialTheta.length;
    }

    /** Total free parameters for AIC/BIC: rate params + (3 empirical frequency params, "+F" convention, if used). */
    public int totalSubstitutionParameterCount() {
        return freeRateParameterCount() + (useEmpiricalFrequencies ? 3 : 0);
    }

    public GtrModel buildModel(double[] theta, double[] empiricalFrequencies) {
        double[] freq = useEmpiricalFrequencies ? empiricalFrequencies : new double[]{0.25, 0.25, 0.25, 0.25};
        double[] rates = thetaToRates.apply(theta);
        return GtrModel.of(freq, rates);
    }
}
