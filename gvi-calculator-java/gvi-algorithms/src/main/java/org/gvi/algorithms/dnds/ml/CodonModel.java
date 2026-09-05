package org.gvi.algorithms.dnds.ml;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.gvi.core.exception.GviInputException;

/**
 * The Goldman &amp; Yang (1994) codon substitution model (GY94) -- the same
 * core continuous-time Markov model PAML's {@code codeml} builds its M0
 * "one-ratio" analysis on. A single instantaneous-rate matrix Q over the 61
 * sense codons, parameterized by:
 * <ul>
 *   <li>{@code kappa}: the transition/transversion rate ratio</li>
 *   <li>{@code omega}: dN/dS itself -- the parameter this whole model exists
 *       to estimate</li>
 *   <li>codon frequencies {@code pi} (F3x4, {@link CodonFrequencies})</li>
 * </ul>
 * Only codon pairs differing at exactly one nucleotide position get a
 * nonzero instantaneous rate (multi-step substitutions in one instant have
 * probability zero under a CTMC and are the standard GY94/codeml
 * convention); for such a pair i -&gt; j:
 * <pre>
 *   Q[i][j] = pi[j] * (kappa if the difference is a transition else 1) * (omega if nonsynonymous else 1)
 * </pre>
 * scaled so the mean instantaneous rate is 1 (branch lengths are then in
 * expected substitutions/codon, same convention as {@code codeml}).
 * <p>
 * This construction is time-reversible by inspection (Q[i][j]*pi[i] ==
 * Q[j][i]*pi[j], since the transition/transversion and synonymous/
 * nonsynonymous classification of the single differing position is the
 * same in both directions), so P(t)=exp(Qt) is computed via the identical
 * symmetric-eigendecomposition trick {@link org.gvi.algorithms.phylo.model.GtrModel}
 * uses for the 4-state nucleotide case, just generalized to N=61 states.
 */
public final class CodonModel {

    private final double[] pi;
    private final double kappa;
    private final double omega;
    private final RealMatrix eigenvectors;
    private final double[] eigenvalues;

    private CodonModel(double[] pi, double kappa, double omega, RealMatrix eigenvectors, double[] eigenvalues) {
        this.pi = pi;
        this.kappa = kappa;
        this.omega = omega;
        this.eigenvectors = eigenvectors;
        this.eigenvalues = eigenvalues;
    }

    public static CodonModel of(double[] pi, double kappa, double omega) {
        int n = CodonAlphabet.SIZE;
        if (pi.length != n) {
            throw new GviInputException("Codon frequency vector must have " + n + " entries, got " + pi.length);
        }
        if (kappa <= 0 || omega <= 0) {
            throw new GviInputException("kappa and omega must be positive (kappa=" + kappa + ", omega=" + omega + ")");
        }

        double[][] q = new double[n][n];
        for (int i = 0; i < n; i++) {
            String ci = CodonAlphabet.SENSE_CODONS.get(i);
            double rowSum = 0;
            for (int j = 0; j < n; j++) {
                if (i == j) continue;
                String cj = CodonAlphabet.SENSE_CODONS.get(j);
                if (CodonAlphabet.ntDifferences(ci, cj) != 1) continue;
                int pos = CodonAlphabet.differingPosition(ci, cj);
                double rate = pi[j];
                if (CodonAlphabet.isTransition(ci.charAt(pos), cj.charAt(pos))) rate *= kappa;
                if (!CodonAlphabet.isSynonymous(ci, cj)) rate *= omega;
                q[i][j] = rate;
                rowSum += rate;
            }
            q[i][i] = -rowSum;
        }

        double meanRate = 0;
        for (int i = 0; i < n; i++) meanRate += -pi[i] * q[i][i];
        if (meanRate <= 0) {
            throw new GviInputException("Codon rate matrix has non-positive mean rate; check frequency/kappa/omega inputs");
        }
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                q[i][j] /= meanRate;
            }
        }

        double[] sqrtPi = new double[n];
        for (int i = 0; i < n; i++) sqrtPi[i] = Math.sqrt(pi[i]);

        double[][] b = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                b[i][j] = sqrtPi[i] * q[i][j] / sqrtPi[j];
            }
        }
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                double avg = (b[i][j] + b[j][i]) / 2.0;
                b[i][j] = avg;
                b[j][i] = avg;
            }
        }

        EigenDecomposition eig = new EigenDecomposition(new Array2DRowRealMatrix(b));
        return new CodonModel(pi.clone(), kappa, omega, eig.getV(), eig.getRealEigenvalues());
    }

    /** P(t): probability of ending in codon j after time t (expected substitutions/codon), starting from codon i. */
    public double[][] transitionProbabilities(double t) {
        int n = CodonAlphabet.SIZE;
        double[][] p = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                double sum = 0;
                for (int k = 0; k < n; k++) {
                    sum += eigenvectors.getEntry(i, k) * eigenvectors.getEntry(j, k) * Math.exp(eigenvalues[k] * t);
                }
                p[i][j] = Math.sqrt(pi[j] / pi[i]) * sum;
            }
        }
        return p;
    }

    public double[] frequencies() {
        return pi.clone();
    }

    public double kappa() {
        return kappa;
    }

    public double omega() {
        return omega;
    }
}
