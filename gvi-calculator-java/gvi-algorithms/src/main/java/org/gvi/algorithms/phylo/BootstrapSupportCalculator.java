package org.gvi.algorithms.phylo;

import org.gvi.algorithms.gd.GdMethod;
import org.gvi.core.exception.GviComputationException;
import org.gvi.core.exception.GviException;
import org.gvi.core.model.NucleotideSequence;
import org.gvi.core.model.SequenceAlignment;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

/**
 * Felsenstein (1985) nonparametric bootstrap support for a Neighbor-Joining
 * tree's internal splits: resample alignment columns with replacement B
 * times, rebuild the NJ tree from each resampled alignment, and report what
 * fraction of replicates recover each of the original tree's bipartitions
 * (clades) -- the same method RAxML/IQ-TREE/PAUP* report as "bootstrap
 * support %" on a tree figure, and something this project's trees didn't
 * carry before.
 * <p>
 * Deliberately bootstraps only the (fast, O(n^3)) NJ topology, not the much
 * slower ML branch-length/substitution-parameter fit (`--high-accuracy-mu`)
 * -- bootstrapping that too would multiply an already-expensive
 * optimization by B replicates. This is exactly why RAxML's own "rapid
 * bootstrap" mode is a separate, cheaper pass from its full ML search
 * rather than literally rerunning ML search B times.
 */
public final class BootstrapSupportCalculator {

    public static final int DEFAULT_REPLICATES = 200;
    /** Bootstrapping multiplies NJ's own tree-building cost (dominated in practice by O(n^2 x alignment length) pairwise distances) by the replicate count -- tighter than the general {@link PhyloTreeFactory#MAX_TAXA_FOR_TREE} cap so a default-replicate run stays tractable. */
    public static final int MAX_TAXA_FOR_BOOTSTRAP = 60;

    private final int replicates;
    private final long seed;

    public BootstrapSupportCalculator() {
        this(DEFAULT_REPLICATES, 42L);
    }

    public BootstrapSupportCalculator(int replicates, long seed) {
        this.replicates = replicates;
        this.seed = seed;
    }

    /** {@code supportByClade} maps each of the original tree's informative splits (as a taxon-label set) to its bootstrap support percentage. */
    public record Support(Map<Set<String>, Double> supportByClade, int successfulReplicates, List<String> diagnostics) {

        public Double forClade(Set<String> taxonLabels) {
            return supportByClade.get(taxonLabels);
        }
    }

    public Support compute(SequenceAlignment alignment, GdMethod method, PhyloTree originalTree) {
        int n = alignment.size();
        List<Set<String>> originalClades = Bipartitions.of(originalTree, n);
        if (originalClades.isEmpty()) {
            return new Support(Map.of(), 0, List.of("No informative internal splits to test (tree too small)"));
        }

        List<NucleotideSequence> sequences = alignment.getSequences();
        int length = alignment.length();
        Random rng = new Random(seed);

        Map<Set<String>, Integer> hits = new LinkedHashMap<>();
        for (Set<String> c : originalClades) hits.put(c, 0);

        int successfulReplicates = 0;
        int failedReplicates = 0;
        for (int b = 0; b < replicates; b++) {
            int[] columns = new int[length];
            for (int i = 0; i < length; i++) columns[i] = rng.nextInt(length);

            List<NucleotideSequence> resampled = new ArrayList<>(sequences.size());
            for (NucleotideSequence s : sequences) {
                String seq = s.getSequence();
                char[] chars = new char[length];
                for (int i = 0; i < length; i++) chars[i] = seq.charAt(columns[i]);
                resampled.add(new NucleotideSequence(s.getId(), new String(chars)));
            }

            try {
                SequenceAlignment resampledAlignment = SequenceAlignment.of(resampled, alignment.getReference().getId());
                PhyloTreeFactory.BuiltTree built = PhyloTreeFactory.build(resampledAlignment, method);
                Set<Set<String>> replicateClades = new HashSet<>(Bipartitions.of(built.tree(), n));
                for (Set<String> c : originalClades) {
                    if (replicateClades.contains(c)) hits.merge(c, 1, Integer::sum);
                }
                successfulReplicates++;
            } catch (GviException e) {
                failedReplicates++; // a pathological resample (e.g. a degenerate all-invariant column set) shouldn't abort the whole bootstrap
            }
        }

        if (successfulReplicates == 0) {
            throw new GviComputationException("All " + replicates + " bootstrap replicates failed to build a tree; cannot compute support");
        }

        Map<Set<String>, Double> support = new LinkedHashMap<>();
        for (var entry : hits.entrySet()) {
            support.put(entry.getKey(), 100.0 * entry.getValue() / successfulReplicates);
        }

        List<String> diagnostics = new ArrayList<>();
        diagnostics.add(successfulReplicates + " of " + replicates + " bootstrap replicate(s) succeeded"
                + (failedReplicates > 0 ? " (" + failedReplicates + " failed to build a tree and were excluded)" : ""));

        return new Support(support, successfulReplicates, diagnostics);
    }
}
