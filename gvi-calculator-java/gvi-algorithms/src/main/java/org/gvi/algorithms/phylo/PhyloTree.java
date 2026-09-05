package org.gvi.algorithms.phylo;

import org.gvi.core.exception.GviComputationException;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Deque;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * A Neighbor-Joining result, rooted at a chosen taxon (typically the
 * designated reference sequence -- a defensible, simple "outgroup rooting"
 * choice when a known reference exists, standard practice when there's no
 * independent outgroup available).
 */
public final class PhyloTree {

    private final PhyloNode root;
    private final Map<String, PhyloNode> taxonByLabel = new HashMap<>();
    private final List<PhyloNode> branchPoints = new ArrayList<>(); // internal (non-taxon) nodes, i.e. true bifurcation events

    private PhyloTree(PhyloNode root) {
        this.root = root;
        collect(root);
    }

    private void collect(PhyloNode node) {
        if (node.isTaxon()) {
            taxonByLabel.put(node.label(), node);
        } else if (!node.isRoot()) {
            branchPoints.add(node);
        }
        for (PhyloNode child : node.children()) collect(child);
    }

    public static PhyloTree rootAt(NeighborJoining.Result njResult, String rootLabel) {
        Map<Integer, List<double[]>> adjacency = new HashMap<>(); // id -> list of [neighborId, length]
        for (PhyloEdge e : njResult.edges()) {
            adjacency.computeIfAbsent(e.a(), k -> new ArrayList<>()).add(new double[]{e.b(), e.length()});
            adjacency.computeIfAbsent(e.b(), k -> new ArrayList<>()).add(new double[]{e.a(), e.length()});
        }

        Integer rootId = null;
        for (var entry : njResult.leafLabels().entrySet()) {
            if (entry.getValue().equals(rootLabel)) {
                rootId = entry.getKey();
                break;
            }
        }
        if (rootId == null) {
            throw new GviComputationException("Cannot root tree: label '" + rootLabel + "' not found among the taxa");
        }

        Map<Integer, PhyloNode> built = new HashMap<>();
        PhyloNode rootNode = new PhyloNode(rootId, njResult.leafLabels().get(rootId));
        built.put(rootId, rootNode);

        Deque<Integer> queue = new ArrayDeque<>();
        queue.add(rootId);
        while (!queue.isEmpty()) {
            int current = queue.poll();
            PhyloNode currentNode = built.get(current);
            for (double[] edge : adjacency.getOrDefault(current, List.of())) {
                int neighborId = (int) edge[0];
                double length = edge[1];
                if (built.containsKey(neighborId)) continue; // already visited (that's where we came from)
                PhyloNode neighborNode = new PhyloNode(neighborId, njResult.leafLabels().get(neighborId));
                neighborNode.setParent(currentNode, length);
                currentNode.addChild(neighborNode);
                built.put(neighborId, neighborNode);
                queue.add(neighborId);
            }
        }

        return new PhyloTree(rootNode);
    }

    public PhyloNode root() {
        return root;
    }

    public double rootToTip(String taxonLabel) {
        PhyloNode node = taxonByLabel.get(taxonLabel);
        if (node == null) {
            throw new GviComputationException("Taxon '" + taxonLabel + "' not found in tree");
        }
        return node.rootToTipDistance();
    }

    public Map<String, Double> allRootToTipDistances() {
        Map<String, Double> result = new HashMap<>();
        for (var entry : taxonByLabel.entrySet()) {
            result.put(entry.getKey(), entry.getValue().rootToTipDistance());
        }
        return result;
    }

    /** Neighbor-Joining's internal branching points (true bifurcation events), excluding the root itself. */
    public List<PhyloNode> branchPoints() {
        return branchPoints;
    }

    public int taxonCount() {
        return taxonByLabel.size();
    }
}
