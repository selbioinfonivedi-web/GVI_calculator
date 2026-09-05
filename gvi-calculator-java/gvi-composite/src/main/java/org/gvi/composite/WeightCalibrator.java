package org.gvi.composite;

import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.NelderMeadSimplex;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.SimplexOptimizer;
import org.gvi.core.exception.GviComputationException;
import org.gvi.core.exception.GviInputException;

import java.util.ArrayList;
import java.util.Collections;
import java.util.EnumMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

/**
 * Weight-calibration module (Section 5.9 / 8.5 of the build spec, marked
 * "optional/advanced" -- fits composite GVI weights against an observed
 * epidemic signal instead of using the spec's fixed default midpoints).
 * <p>
 * Fits weights on the probability simplex (each w_i &gt;= 0, sum(w) = 1) by
 * minimizing sum-of-squared-residuals between the weighted-normalized-index
 * dot product and an observed target (e.g. an independently measured growth
 * rate or attack-rate proxy for each historical sample/timepoint). The
 * simplex constraint is enforced structurally via a softmax
 * reparameterization (w_i = exp(theta_i) / sum_j exp(theta_j)) rather than
 * a constrained optimizer, so plain derivative-free Nelder-Mead
 * (Apache Commons Math's {@link SimplexOptimizer}) can be used directly on
 * the unconstrained theta.
 */
public final class WeightCalibrator {

    public record Observation(Map<IndexKey, Double> normalizedValues, double observedTarget) {
    }

    public record CalibrationResult(CompositeWeights weights, double trainRmse, Double holdoutRmse, int observationsUsed) {
    }

    /** Fits weights against all observations, reporting in-sample RMSE. */
    public CalibrationResult calibrate(List<Observation> observations, List<IndexKey> keysInPlay) {
        double[] theta = fit(observations, keysInPlay);
        double[] w = softmax(theta);
        double rmse = Math.sqrt(sumSquaredResiduals(w, observations, keysInPlay) / observations.size());
        return new CalibrationResult(toWeights(keysInPlay, w), rmse, null, observations.size());
    }

    /**
     * Hold-out validated fit (Section 8.5: "Hold-out validation: predict
     * 20% of data; compare predicted vs. observed"). Splits observations
     * (seeded shuffle for reproducibility), fits on the training split,
     * reports RMSE on both splits.
     */
    /**
     * Result of leave-one-group-out cross-validation: how well weights fitted on every OTHER group
     * predict the group that was held out.
     *
     * @param meanHeldOutRmse averaged across folds -- the headline generalisation figure
     * @param perGroupRmse    each held-out group's own error, so a single pathological group is visible
     * @param weights         weights refitted on ALL observations, for actual use
     */
    public record GroupCvResult(CompositeWeights weights, double trainRmse, double meanHeldOutRmse,
                                Map<String, Double> perGroupRmse) {
    }

    /**
     * Leave-one-group-out cross-validation, where a group is normally one pathogen.
     * <p>
     * A random holdout answers "can these weights interpolate among observations I have already
     * seen?" -- which, on a corpus of a few observations per pathogen, is nearly free: the same
     * pathogen appears on both sides of the split, so the model can succeed by memorising it. The
     * question that actually matters for a general-purpose index is "do weights fitted on other
     * pathogens predict a pathogen never seen during fitting?", and only holding out whole groups
     * asks it.
     * <p>
     * The gap between {@code trainRmse} and {@code meanHeldOutRmse} is the diagnostic: close
     * together means the weights carry transferable signal; a large gap means they have been fitted
     * to the idiosyncrasies of the training pathogens and should not be applied to a new one.
     *
     * @param groupOf group label per observation, positionally aligned with {@code observations}
     */
    public GroupCvResult calibrateWithGroupHoldout(List<Observation> observations, List<IndexKey> keysInPlay,
                                                   List<String> groupOf) {
        if (observations.size() != groupOf.size()) {
            throw new GviInputException("groupOf must have one label per observation (" + observations.size()
                    + " observations, " + groupOf.size() + " labels)");
        }
        java.util.LinkedHashSet<String> groups = new java.util.LinkedHashSet<>(groupOf);
        if (groups.size() < 2) {
            throw new GviInputException("Leave-one-group-out needs at least 2 distinct groups, found " + groups.size()
                    + ". With one group there is nothing to generalize to -- use calibrateWithHoldout instead.");
        }

        Map<String, Double> perGroup = new java.util.LinkedHashMap<>();
        for (String held : groups) {
            List<Observation> train = new ArrayList<>();
            List<Observation> test = new ArrayList<>();
            for (int i = 0; i < observations.size(); i++) {
                (groupOf.get(i).equals(held) ? test : train).add(observations.get(i));
            }
            if (train.isEmpty() || test.isEmpty()) continue;
            double[] w = softmax(fit(train, keysInPlay));
            perGroup.put(held, Math.sqrt(sumSquaredResiduals(w, test, keysInPlay) / test.size()));
        }
        if (perGroup.isEmpty()) {
            throw new GviInputException("No usable leave-one-group-out fold could be built");
        }

        double[] all = softmax(fit(observations, keysInPlay));
        double trainRmse = Math.sqrt(sumSquaredResiduals(all, observations, keysInPlay) / observations.size());
        double meanHeldOut = perGroup.values().stream().mapToDouble(Double::doubleValue).average().orElse(Double.NaN);
        return new GroupCvResult(toWeights(keysInPlay, all), trainRmse, meanHeldOut, Map.copyOf(perGroup));
    }

    public CalibrationResult calibrateWithHoldout(List<Observation> observations, List<IndexKey> keysInPlay,
                                                   double holdoutFraction, long seed) {
        if (holdoutFraction <= 0 || holdoutFraction >= 1) {
            throw new GviInputException("holdoutFraction must be in (0,1), got " + holdoutFraction);
        }
        List<Observation> shuffled = new ArrayList<>(observations);
        Collections.shuffle(shuffled, new Random(seed));
        int holdoutSize = Math.max(1, (int) Math.round(shuffled.size() * holdoutFraction));
        if (holdoutSize >= shuffled.size()) {
            throw new GviInputException("Not enough observations (" + shuffled.size() + ") to hold out " + holdoutFraction + " and still have a training set");
        }
        List<Observation> holdout = shuffled.subList(0, holdoutSize);
        List<Observation> train = shuffled.subList(holdoutSize, shuffled.size());

        double[] theta = fit(train, keysInPlay);
        double[] w = softmax(theta);
        double trainRmse = Math.sqrt(sumSquaredResiduals(w, train, keysInPlay) / train.size());
        double holdoutRmse = Math.sqrt(sumSquaredResiduals(w, holdout, keysInPlay) / holdout.size());

        return new CalibrationResult(toWeights(keysInPlay, w), trainRmse, holdoutRmse, train.size());
    }

    private double[] fit(List<Observation> observations, List<IndexKey> keysInPlay) {
        if (observations == null || observations.size() < 2) {
            throw new GviInputException("Need at least 2 observations to calibrate weights, got "
                    + (observations == null ? 0 : observations.size()));
        }
        if (keysInPlay == null || keysInPlay.isEmpty()) {
            throw new GviInputException("Need at least one index key to calibrate weights against");
        }
        int k = keysInPlay.size();

        MultivariateFunction objective = theta -> sumSquaredResiduals(softmax(theta), observations, keysInPlay);

        SimplexOptimizer optimizer = new SimplexOptimizer(1e-10, 1e-12);
        NelderMeadSimplex simplex = new NelderMeadSimplex(k);
        double[] initialTheta = new double[k]; // all-zero theta -> uniform weights as the starting point

        try {
            PointValuePair result = optimizer.optimize(
                    new MaxEval(20_000),
                    new ObjectiveFunction(objective),
                    GoalType.MINIMIZE,
                    new InitialGuess(initialTheta),
                    simplex);
            return result.getPoint();
        } catch (RuntimeException e) {
            throw new GviComputationException("Weight calibration failed to converge: " + e.getMessage(), e);
        }
    }

    private double sumSquaredResiduals(double[] weights, List<Observation> observations, List<IndexKey> keysInPlay) {
        double sse = 0;
        for (Observation obs : observations) {
            double predicted = 0;
            for (int i = 0; i < keysInPlay.size(); i++) {
                predicted += weights[i] * obs.normalizedValues().getOrDefault(keysInPlay.get(i), 0.0);
            }
            double residual = predicted - obs.observedTarget();
            sse += residual * residual;
        }
        return sse;
    }

    private double[] softmax(double[] theta) {
        double max = Double.NEGATIVE_INFINITY;
        for (double t : theta) max = Math.max(max, t);
        double[] exp = new double[theta.length];
        double sum = 0;
        for (int i = 0; i < theta.length; i++) {
            exp[i] = Math.exp(theta[i] - max); // subtract max for numerical stability
            sum += exp[i];
        }
        double[] w = new double[theta.length];
        for (int i = 0; i < theta.length; i++) w[i] = exp[i] / sum;
        return w;
    }

    private CompositeWeights toWeights(List<IndexKey> keysInPlay, double[] w) {
        Map<IndexKey, Double> map = new EnumMap<>(IndexKey.class);
        for (int i = 0; i < keysInPlay.size(); i++) map.put(keysInPlay.get(i), w[i]);
        return CompositeWeights.of(map);
    }
}
