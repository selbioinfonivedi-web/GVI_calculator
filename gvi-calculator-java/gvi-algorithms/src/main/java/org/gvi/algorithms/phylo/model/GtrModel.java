package org.gvi.algorithms.phylo.model;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.gvi.core.exception.GviInputException;

/**
 * The General Time-Reversible (GTR) nucleotide substitution model -- the
 * model RAxML/IQ-TREE default to, and a strict generalization of
 * Jukes-Cantor (equal frequencies, equal rates) and Kimura 2-Parameter
 * (equal frequencies, transition/transversion rate split) used elsewhere in
 * this codebase for simple pairwise distances. GTR allows each of the 4
 * base frequencies and each of the 6 pairwise exchangeability rates
 * (AC, AG, AT, CG, CT, GT) to differ independently, which matters for real
 * sequence data that rarely has exactly equal composition or rates.
 * <p>
 * Transition probabilities P(t) = exp(Qt) are computed via the standard
 * symmetric-eigendecomposition trick for time-reversible rate matrices
 * (Q is similar to a symmetric matrix B = Pi^(1/2) Q Pi^(-1/2), which has
 * real eigenvalues/orthogonal eigenvectors by construction, avoiding the
 * numerical fragility of a general complex eigendecomposition):
 * <pre>
 *   P(t)[i][j] = sqrt(pi_j / pi_i) * sum_k V[i][k] * V[j][k] * exp(lambda_k * t)
 * </pre>
 * where V, lambda are the eigenvectors/eigenvalues of the symmetrized B.
 */
public final class GtrModel {

    /** exchangeability rate order: AC, AG, AT, CG, CT, GT */
    private final double[] baseFrequencies; // [piA, piC, piG, piT]
    private final RealMatrix eigenvectors;  // V, columns are eigenvectors of the symmetrized matrix
    private final double[] eigenvalues;     // lambda_k

    private GtrModel(double[] baseFrequencies, RealMatrix eigenvectors, double[] eigenvalues) {
        this.baseFrequencies = baseFrequencies;
        this.eigenvectors = eigenvectors;
        this.eigenvalues = eigenvalues;
    }

    public static GtrModel of(double[] baseFrequencies, double[] exchangeabilityRates) {
        if (baseFrequencies.length != 4) {
            throw new GviInputException("GTR base frequencies must have exactly 4 entries (A,C,G,T), got " + baseFrequencies.length);
        }
        if (exchangeabilityRates.length != 6) {
            throw new GviInputException("GTR exchangeability rates must have exactly 6 entries (AC,AG,AT,CG,CT,GT), got " + exchangeabilityRates.length);
        }
        double freqSum = 0;
        for (double f : baseFrequencies) {
            if (f <= 0) throw new GviInputException("GTR base frequencies must all be positive");
            freqSum += f;
        }
        if (Math.abs(freqSum - 1.0) > 1e-6) {
            throw new GviInputException("GTR base frequencies must sum to 1.0, got " + freqSum);
        }
        for (double r : exchangeabilityRates) {
            if (r <= 0) throw new GviInputException("GTR exchangeability rates must all be positive");
        }

        double[][] r = symmetricRateMatrix(exchangeabilityRates);
        double[][] q = buildRateMatrix(baseFrequencies, r);
        normalizeToUnitMeanRate(q, baseFrequencies);

        double[] sqrtPi = new double[4];
        for (int i = 0; i < 4; i++) sqrtPi[i] = Math.sqrt(baseFrequencies[i]);

        double[][] b = new double[4][4];
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                b[i][j] = sqrtPi[i] * q[i][j] / sqrtPi[j];
            }
        }
        // b should be exactly symmetric by construction; average with transpose to kill floating-point asymmetry
        for (int i = 0; i < 4; i++) {
            for (int j = i + 1; j < 4; j++) {
                double avg = (b[i][j] + b[j][i]) / 2.0;
                b[i][j] = avg;
                b[j][i] = avg;
            }
        }

        EigenDecomposition eig = new EigenDecomposition(new Array2DRowRealMatrix(b));
        double[] lambda = eig.getRealEigenvalues();
        RealMatrix v = eig.getV();

        return new GtrModel(baseFrequencies.clone(), v, lambda);
    }

    /** JC69 as a degenerate GTR: equal frequencies, equal rates. Useful as a sanity-check baseline. */
    public static GtrModel jukesCantor() {
        return of(new double[]{0.25, 0.25, 0.25, 0.25}, new double[]{1, 1, 1, 1, 1, 1});
    }

    private static double[][] symmetricRateMatrix(double[] rates) {
        double[][] r = new double[4][4];
        // order: AC, AG, AT, CG, CT, GT
        r[Nucleotide.A][Nucleotide.C] = r[Nucleotide.C][Nucleotide.A] = rates[0];
        r[Nucleotide.A][Nucleotide.G] = r[Nucleotide.G][Nucleotide.A] = rates[1];
        r[Nucleotide.A][Nucleotide.T] = r[Nucleotide.T][Nucleotide.A] = rates[2];
        r[Nucleotide.C][Nucleotide.G] = r[Nucleotide.G][Nucleotide.C] = rates[3];
        r[Nucleotide.C][Nucleotide.T] = r[Nucleotide.T][Nucleotide.C] = rates[4];
        r[Nucleotide.G][Nucleotide.T] = r[Nucleotide.T][Nucleotide.G] = rates[5];
        return r;
    }

    private static double[][] buildRateMatrix(double[] freq, double[][] r) {
        double[][] q = new double[4][4];
        for (int i = 0; i < 4; i++) {
            double rowSum = 0;
            for (int j = 0; j < 4; j++) {
                if (i == j) continue;
                q[i][j] = r[i][j] * freq[j];
                rowSum += q[i][j];
            }
            q[i][i] = -rowSum;
        }
        return q;
    }

    private static void normalizeToUnitMeanRate(double[][] q, double[] freq) {
        double meanRate = 0;
        for (int i = 0; i < 4; i++) meanRate += -freq[i] * q[i][i];
        if (meanRate <= 0) throw new GviInputException("GTR rate matrix has non-positive mean rate; check inputs");
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                q[i][j] /= meanRate;
            }
        }
    }

    /** P(t): probability of ending in state j after time t, starting from state i. Branch length t in expected substitutions/site. */
    public double[][] transitionProbabilities(double t) {
        double[][] p = new double[4][4];
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                double sum = 0;
                for (int k = 0; k < 4; k++) {
                    sum += eigenvectors.getEntry(i, k) * eigenvectors.getEntry(j, k) * Math.exp(eigenvalues[k] * t);
                }
                p[i][j] = Math.sqrt(baseFrequencies[j] / baseFrequencies[i]) * sum;
            }
        }
        return p;
    }

    public double[] baseFrequencies() {
        return baseFrequencies.clone();
    }
}
