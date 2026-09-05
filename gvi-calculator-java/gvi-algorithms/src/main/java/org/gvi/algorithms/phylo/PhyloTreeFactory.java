package org.gvi.algorithms.phylo;

import org.gvi.algorithms.gd.GdMethod;
import org.gvi.algorithms.gd.GdResult;
import org.gvi.algorithms.gd.GeneticDistanceCalculator;
import org.gvi.core.exception.GviComputationException;
import org.gvi.core.model.NucleotideSequence;
import org.gvi.core.model.SequenceAlignment;
import org.gvi.core.util.TemporalUtil;

import java.util.ArrayList;
import java.util.List;

/**
 * Shared entrypoint for building a Neighbor-Joining tree from a
 * {@link SequenceAlignment}, rooted at its designated reference -- used by
 * both tree-aware mu (root-to-tip regression) and the lineage-through-time
 * Re estimator so the (fairly expensive, O(n^3)) tree is only built once
 * per caller and the size guard/diagnostics are consistent across both.
 */
public final class PhyloTreeFactory {

    /** NJ is O(n^3); above this taxon count building a tree is not attempted. */
    public static final int MAX_TAXA_FOR_TREE = 300;

    private PhyloTreeFactory() {
    }

    public record BuiltTree(PhyloTree tree, List<String> diagnostics) {
    }

    public static BuiltTree build(SequenceAlignment alignment, GdMethod method) {
        return build(alignment, method, alignment.getReference().getId());
    }

    /**
     * Same as {@link #build(SequenceAlignment, GdMethod)}, but rooted at an
     * explicit taxon instead of always the pipeline's {@code --reference-id}.
     * Callers that regress root-to-tip distance against sampling date (mu,
     * Re) should root at {@link #temporalAnchorId} instead of the declared
     * reference: rooting at an arbitrary tip -- especially a recent one --
     * makes root-to-tip distance reflect that tip's own position in the tree
     * rather than accumulated divergence since a temporally basal point, and
     * can flip a real molecular-clock signal's sign or discard it entirely.
     * Callers that don't regress against date (dN/dS's ML tree, bootstrap
     * support, model selection) are unaffected by rooting choice and can
     * keep using the declared reference.
     */
    public static BuiltTree build(SequenceAlignment alignment, GdMethod method, String rootSequenceId) {
        if (alignment.size() > MAX_TAXA_FOR_TREE) {
            throw new GviComputationException(alignment.size() + " taxa exceeds the " + MAX_TAXA_FOR_TREE
                    + "-taxon cap for Neighbor-Joining (O(n^3) cost); use a non-tree-based estimator instead");
        }
        if (alignment.size() < 3) {
            throw new GviComputationException("Building a tree needs at least 3 sequences, got " + alignment.size());
        }

        List<NucleotideSequence> sequences = alignment.getSequences();
        int n = sequences.size();
        String[] labels = new String[n];
        double[][] distances = new double[n][n];
        List<String> diagnostics = new ArrayList<>();
        GeneticDistanceCalculator distanceCalculator = new GeneticDistanceCalculator();
        int saturatedPairs = 0;

        for (int i = 0; i < n; i++) labels[i] = sequences.get(i).getId();
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                GdResult gd = distanceCalculator.compute(sequences.get(i), sequences.get(j), method);
                double d = (gd.saturated() || gd.distance() == null) ? gd.pDistance() : gd.distance();
                if (gd.saturated()) saturatedPairs++;
                distances[i][j] = d;
                distances[j][i] = d;
            }
        }
        if (saturatedPairs > 0) {
            diagnostics.add(saturatedPairs + " pair(s) had saturated distance; raw p-distance used for those pairs when building the tree");
        }

        NeighborJoining.Result njResult = NeighborJoining.build(labels, distances);
        diagnostics.addAll(njResult.diagnostics());
        PhyloTree tree = PhyloTree.rootAt(njResult, rootSequenceId);
        diagnostics.add("Built a " + n + "-taxon Neighbor-Joining tree rooted at '" + rootSequenceId + "'");

        return new BuiltTree(tree, diagnostics);
    }

    /**
     * The taxon that should root a root-to-tip date regression: the
     * earliest-dated sequence in the alignment, since only a temporally
     * basal root makes "distance from root" track "time since divergence".
     * Falls back to the pipeline reference when no sequence is dated, or
     * when the reference already is the earliest-dated sequence (the common
     * case, and a no-op). Adds a diagnostic explaining the substitution
     * whenever it actually changes which taxon roots the tree, so a user
     * relying on {@code --reference-id} for GD/GC/CAI isn't surprised that
     * mu/Re used a different taxon to root their own regression.
     */
    public static String temporalAnchorId(SequenceAlignment alignment, List<String> diagnostics) {
        NucleotideSequence declaredReference = alignment.getReference();
        NucleotideSequence earliest = null;
        double earliestYear = Double.POSITIVE_INFINITY;
        for (NucleotideSequence seq : alignment.getSequences()) {
            if (seq.getCollectionDate().isEmpty()) continue;
            double year = TemporalUtil.toDecimalYear(seq.getCollectionDate().get());
            if (year < earliestYear) {
                earliestYear = year;
                earliest = seq;
            }
        }
        if (earliest == null || earliest.getId().equals(declaredReference.getId())) {
            return declaredReference.getId();
        }
        diagnostics.add(String.format(
                "Rooting the temporal (root-to-tip) regression at '%s' (earliest collection date, %.2f) instead of the "
                        + "pipeline reference '%s' -- an arbitrary or recent --reference-id can otherwise invert or discard "
                        + "a real molecular-clock signal. GD, GC deviation, and CAI are unaffected and still use '%s'.",
                earliest.getId(), earliestYear, declaredReference.getId(), declaredReference.getId()));
        return earliest.getId();
    }
}
