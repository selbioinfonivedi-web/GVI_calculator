package org.gvi.algorithms.phylo;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Extracts the informative bipartitions (clades) of a reference-rooted
 * binary {@link PhyloTree}, as a taxon-label set per internal split -- the
 * standard representation for comparing whether two independently built
 * trees over the same taxon set agree on a split, regardless of drawing
 * order/orientation (used by {@link BootstrapSupportCalculator}).
 */
final class Bipartitions {

    private Bipartitions() {
    }

    /**
     * One entry per internal branch point, excluding the trivial split that
     * always exists by construction (root's own child's subtree = every
     * taxon except the rooting reference) -- that split carries no
     * topological information since every replicate is rooted at the same
     * reference.
     */
    static List<Set<String>> of(PhyloTree tree, int totalTaxa) {
        List<Set<String>> result = new ArrayList<>();
        for (PhyloNode branchPoint : tree.branchPoints()) {
            Set<String> taxa = new HashSet<>();
            collectTaxa(branchPoint, taxa);
            if (taxa.size() >= 2 && taxa.size() <= totalTaxa - 2) {
                result.add(taxa);
            }
        }
        return result;
    }

    private static void collectTaxa(PhyloNode node, Set<String> out) {
        if (node.isTaxon()) {
            out.add(node.label());
            return;
        }
        for (PhyloNode child : node.children()) {
            collectTaxa(child, out);
        }
    }
}
