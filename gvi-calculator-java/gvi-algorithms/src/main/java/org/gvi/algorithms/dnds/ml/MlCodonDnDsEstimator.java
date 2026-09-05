package org.gvi.algorithms.dnds.ml;

import org.gvi.algorithms.common.CodonUtil;
import org.gvi.algorithms.gd.GdMethod;
import org.gvi.algorithms.phylo.PhyloTree;
import org.gvi.algorithms.phylo.PhyloTreeFactory;
import org.gvi.core.exception.GviComputationException;
import org.gvi.core.model.GeneAnnotation;
import org.gvi.core.model.NucleotideSequence;
import org.gvi.core.model.SequenceAlignment;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Public entry point for the native GY94/{@code codeml} M0-equivalent ML
 * dN/dS estimator (see {@link CodonModel}, {@link CodonTreeLikelihood},
 * {@link CodonMlFitter}): builds a Neighbor-Joining tree, extracts each
 * taxon's codon sequence for one gene, ML-fits kappa and omega jointly with
 * every branch length under the GY94 codon substitution model, and returns
 * omega as a genuine maximum-likelihood dN/dS estimate -- as distinct from
 * this project's default Nei-Gojobori counting-based dN/dS
 * ({@link org.gvi.algorithms.dnds.DnDsCalculator}), which remains the
 * default composite-facing method (already validated, much cheaper).
 * <p>
 * Capped much more tightly than the nucleotide ML paths
 * ({@code --high-accuracy-mu}'s 40-taxon/50,000-site cap): the codon model
 * has 61 states instead of 4, so the same per-site pruning cost is roughly
 * (61/4)^2 =~ 232x higher, and every Brent evaluation during branch-length
 * optimization still recomputes the full-tree likelihood (the same
 * established, if computationally simple, pattern
 * {@link org.gvi.algorithms.phylo.model.MlBranchLengthOptimizer} uses for
 * nucleotides). These caps keep runtime bounded on typical desktop
 * hardware; genes/datasets that exceed them should fall back to
 * Nei-Gojobori (the caller in {@code GviPipeline} does this automatically).
 */
public final class MlCodonDnDsEstimator {

    public static final int MAX_TAXA_FOR_ML_DNDS = 15;
    public static final int MAX_CODONS_FOR_ML_DNDS = 300;

    private static final double INITIAL_KAPPA = 2.0;
    private static final double INITIAL_OMEGA = 0.5;

    public MlDnDsResult compute(SequenceAlignment alignment, GeneAnnotation gene) {
        if (alignment.size() > MAX_TAXA_FOR_ML_DNDS) {
            throw new GviComputationException(alignment.size() + " taxa exceeds the " + MAX_TAXA_FOR_ML_DNDS
                    + "-taxon cap for ML dN/dS (61-state codon likelihood is expensive; see MlCodonDnDsEstimator's class javadoc)");
        }

        PhyloTree tree = PhyloTreeFactory.build(alignment, GdMethod.JUKES_CANTOR).tree();

        Map<String, List<String>> codonsByLabel = new HashMap<>();
        int siteCount = -1;
        List<String> allObservedCodons = new ArrayList<>();
        for (NucleotideSequence seq : alignment.getSequences()) {
            String cds = CodonUtil.extractCds(seq.getSequence(), gene.start(), gene.end(), gene.strand());
            CodonUtil.SplitResult split = CodonUtil.splitIntoCodons(cds);
            List<String> codons = split.codons();
            if (siteCount == -1) {
                siteCount = codons.size();
            } else if (codons.size() != siteCount) {
                throw new GviComputationException("Gene '" + gene.geneName() + "': codon count mismatch across taxa ("
                        + seq.getId() + " has " + codons.size() + ", expected " + siteCount + ")");
            }
            codonsByLabel.put(seq.getId(), codons);
            for (String c : codons) {
                if (CodonAlphabet.isSenseCodon(c)) allObservedCodons.add(c);
            }
        }

        if (siteCount > MAX_CODONS_FOR_ML_DNDS) {
            throw new GviComputationException("Gene '" + gene.geneName() + "' has " + siteCount + " codons, exceeding the "
                    + MAX_CODONS_FOR_ML_DNDS + "-codon cap for ML dN/dS (61-state codon likelihood is expensive; see "
                    + "MlCodonDnDsEstimator's class javadoc)");
        }
        if (allObservedCodons.size() < 30) {
            throw new GviComputationException("Gene '" + gene.geneName() + "' has too few usable (sense) codons ("
                    + allObservedCodons.size() + ") to reliably estimate F3x4 codon frequencies");
        }

        double[] pi = CodonFrequencies.f3x4(allObservedCodons);
        CodonMlFitter.Result fit = new CodonMlFitter().optimize(tree, codonsByLabel, siteCount, pi, INITIAL_KAPPA, INITIAL_OMEGA);

        return new MlDnDsResult(gene.geneName(), fit.omega(), fit.kappa(), fit.logLikelihood(),
                alignment.size(), siteCount, fit.cyclesUsed());
    }
}
