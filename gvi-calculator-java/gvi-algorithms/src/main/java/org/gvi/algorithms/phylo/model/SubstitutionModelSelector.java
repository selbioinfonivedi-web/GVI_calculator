package org.gvi.algorithms.phylo.model;

import org.gvi.algorithms.gd.GdMethod;
import org.gvi.algorithms.phylo.PhyloTree;
import org.gvi.algorithms.phylo.PhyloTreeFactory;
import org.gvi.core.exception.GviComputationException;
import org.gvi.core.model.NucleotideSequence;
import org.gvi.core.model.SequenceAlignment;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Fits each candidate substitution model (see {@link NamedSubstitutionModels})
 * to the data via {@link MlPhylogeneticOptimizer} and picks the best one by
 * AIC (Akaike Information Criterion) -- the same model-selection approach
 * jModelTest/ModelTest-NG/IQ-TREE's ModelFinder use, scoped down to a
 * deterministic, dependency-free Java implementation. This answers "what if
 * the alignment doesn't actually fit a GTR model" directly: don't assume,
 * fit several models and let the data decide.
 * <p>
 * AIC = -2*logLikelihood + 2*k; BIC = -2*logLikelihood + k*ln(sites), where
 * k = substitution-model free parameters (frequencies "+F" convention +
 * rate parameters + alpha if Gamma is enabled) + branch lengths. Since
 * every candidate shares the same topology and branch count, the branch
 * term is a constant across models and does not affect which model wins --
 * it is included for textbook-correct absolute AIC/BIC values.
 * <p>
 * This is, unavoidably, expensive: a full joint branch-length +
 * substitution-parameter ML optimization for each of up to 6 candidate
 * models. It inherits {@link EvolutionaryRateCalculator}'s size caps via
 * whatever caller wires it in.
 */
public final class SubstitutionModelSelector {

    public record ModelFit(String modelName, double logLikelihood, int freeParameters, double aic, double bic,
                            double[] optimizedTheta, Double optimizedAlpha) {
    }

    public record SelectionResult(List<ModelFit> allFits, ModelFit best, PhyloTree treeUnderBestModel,
                                   GtrModel bestModel, double[] bestGammaRates) {
    }

    public SelectionResult select(SequenceAlignment alignment, List<RateParameterization> candidates, int gammaCategoryCount) {
        return select(alignment, candidates, gammaCategoryCount, alignment.getReference().getId());
    }

    /**
     * Same as {@link #select(SequenceAlignment, List, int)}, but rooting every candidate's tree at an explicit
     * taxon instead of always the pipeline's {@code --reference-id}. Callers that go on to regress root-to-tip
     * distance against sampling date (mu's {@code --substitution-model auto} path) should pass
     * {@link PhyloTreeFactory#temporalAnchorId} here -- rooting is irrelevant to the ML fit itself (branch-length
     * optimization under a time-reversible model is root-invariant), so this only changes what the final
     * regression sees, not which model wins or how well it fits.
     */
    public SelectionResult select(SequenceAlignment alignment, List<RateParameterization> candidates, int gammaCategoryCount, String rootSequenceId) {
        if (candidates == null || candidates.isEmpty()) {
            throw new GviComputationException("Cannot select a substitution model: no candidate models supplied");
        }
        GtrParameterEstimator.Estimate estimate = GtrParameterEstimator.estimate(alignment);
        Map<String, String> sequences = new HashMap<>();
        for (NucleotideSequence seq : alignment.getSequences()) sequences.put(seq.getId(), seq.getSequence());

        List<ModelFit> fits = new ArrayList<>();
        ModelFit best = null;
        PhyloTree bestTree = null;
        GtrModel bestModel = null;
        double[] bestGammaRates = null;
        int branchCount = -1;

        for (RateParameterization candidate : candidates) {
            // NJ topology is identical every time (same data), but branch lengths get mutated in place by the
            // optimizer -- rebuilding fresh per candidate is simpler and cheap relative to the ML optimization itself.
            PhyloTreeFactory.BuiltTree built = PhyloTreeFactory.build(alignment, GdMethod.JUKES_CANTOR, rootSequenceId);
            if (branchCount < 0) branchCount = MlBranchLengthOptimizer.collectBranches(built.tree()).size();

            MlPhylogeneticOptimizer.Result result = new MlPhylogeneticOptimizer().optimize(
                    built.tree(), sequences, alignment.length(), candidate, estimate.baseFrequencies(), gammaCategoryCount, 1.0);

            int freeParams = candidate.totalSubstitutionParameterCount() + (gammaCategoryCount >= 2 ? 1 : 0) + branchCount;
            double aic = -2 * result.finalLogLikelihood() + 2 * freeParams;
            double bic = -2 * result.finalLogLikelihood() + freeParams * Math.log(alignment.length());
            ModelFit fit = new ModelFit(candidate.modelName(), result.finalLogLikelihood(), freeParams, aic, bic,
                    result.optimizedTheta(), result.optimizedAlpha());
            fits.add(fit);

            if (best == null || aic < best.aic()) {
                best = fit;
                bestTree = built.tree();
                bestModel = result.finalModel();
                bestGammaRates = result.finalGammaRates();
            }
        }

        fits.sort(Comparator.comparingDouble(ModelFit::aic));
        return new SelectionResult(fits, best, bestTree, bestModel, bestGammaRates);
    }
}
