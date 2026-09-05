package org.gvi.algorithms.phylo;

import org.gvi.core.exception.GviComputationException;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Neighbor-Joining tree construction (Saitou &amp; Nei, 1987). Builds an
 * unrooted bifurcating tree from a pairwise distance matrix -- this is what
 * backs tree-aware mu (root-to-tip regression along real branches instead
 * of raw distance-to-a-single-reference) and the lineage-through-time Re
 * estimator.
 * <p>
 * NJ is provably exact when the input distance matrix is perfectly
 * additive (i.e. genuinely tree-derived): it recovers both the true
 * topology and the true branch lengths. That property is what the unit
 * test for this class exploits -- a distance matrix is constructed by hand
 * from a known tree, and NJ is checked to reconstruct it exactly.
 * <p>
 * Negative branch lengths (a well-known NJ artifact on noisy/non-additive
 * real-world distances) are clamped to zero, matching standard practice in
 * reference implementations (e.g. PHYLIP's neighbor).
 */
public final class NeighborJoining {

    private NeighborJoining() {
    }

    public record Result(List<PhyloEdge> edges, Map<Integer, String> leafLabels, List<String> diagnostics) {
    }

    public static Result build(String[] labels, double[][] distances) {
        int n = labels.length;
        if (n < 3) {
            throw new GviComputationException("Neighbor-Joining needs at least 3 taxa, got " + n);
        }
        validateMatrix(labels, distances);

        List<String> diagnostics = new ArrayList<>();
        List<PhyloEdge> edges = new ArrayList<>();
        Map<Integer, String> leafLabels = new HashMap<>();

        // active node id -> (other active node id -> distance)
        Map<Integer, Map<Integer, Double>> d = new HashMap<>();
        int nextId = n;
        for (int i = 0; i < n; i++) {
            leafLabels.put(i, labels[i]);
            Map<Integer, Double> row = new HashMap<>();
            for (int j = 0; j < n; j++) {
                if (i != j) row.put(j, distances[i][j]);
            }
            d.put(i, row);
        }

        while (d.size() > 2) {
            int activeCount = d.size();
            Map<Integer, Double> r = new HashMap<>();
            for (var entry : d.entrySet()) {
                r.put(entry.getKey(), entry.getValue().values().stream().mapToDouble(Double::doubleValue).sum());
            }

            int bestI = -1, bestJ = -1;
            double bestQ = Double.POSITIVE_INFINITY;
            List<Integer> ids = new ArrayList<>(d.keySet());
            for (int a = 0; a < ids.size(); a++) {
                for (int b = a + 1; b < ids.size(); b++) {
                    int i = ids.get(a), j = ids.get(b);
                    double q = (activeCount - 2) * d.get(i).get(j) - r.get(i) - r.get(j);
                    if (q < bestQ) {
                        bestQ = q;
                        bestI = i;
                        bestJ = j;
                    }
                }
            }

            double dij = d.get(bestI).get(bestJ);
            double di = 0.5 * dij + (activeCount > 2 ? (r.get(bestI) - r.get(bestJ)) / (2.0 * (activeCount - 2)) : 0.0);
            double dj = dij - di;
            if (di < 0) {
                diagnostics.add("Negative branch length (" + di + ") clamped to 0 for leaf/clade '" + labelOf(bestI, leafLabels) + "'");
                di = 0;
            }
            if (dj < 0) {
                diagnostics.add("Negative branch length (" + dj + ") clamped to 0 for leaf/clade '" + labelOf(bestJ, leafLabels) + "'");
                dj = 0;
            }

            int u = nextId++;
            edges.add(new PhyloEdge(bestI, u, di));
            edges.add(new PhyloEdge(bestJ, u, dj));

            Map<Integer, Double> newRow = new HashMap<>();
            for (int k : ids) {
                if (k == bestI || k == bestJ) continue;
                double dk = (d.get(bestI).get(k) + d.get(bestJ).get(k) - dij) / 2.0;
                newRow.put(k, dk);
                d.get(k).put(u, dk);
                d.get(k).remove(bestI);
                d.get(k).remove(bestJ);
            }
            d.remove(bestI);
            d.remove(bestJ);
            d.put(u, newRow);
        }

        // exactly 2 active nodes left -- connect them with the one remaining edge
        List<Integer> last = new ArrayList<>(d.keySet());
        double finalLength = d.get(last.get(0)).get(last.get(1));
        if (finalLength < 0) {
            diagnostics.add("Negative final branch length (" + finalLength + ") clamped to 0");
            finalLength = 0;
        }
        edges.add(new PhyloEdge(last.get(0), last.get(1), finalLength));

        return new Result(edges, leafLabels, diagnostics);
    }

    private static String labelOf(int id, Map<Integer, String> leafLabels) {
        return leafLabels.getOrDefault(id, "internal#" + id);
    }

    private static void validateMatrix(String[] labels, double[][] distances) {
        int n = labels.length;
        if (distances.length != n) {
            throw new GviComputationException("Distance matrix row count (" + distances.length + ") doesn't match label count (" + n + ")");
        }
        for (int i = 0; i < n; i++) {
            if (distances[i].length != n) {
                throw new GviComputationException("Distance matrix row " + i + " has " + distances[i].length + " columns, expected " + n);
            }
            for (int j = 0; j < n; j++) {
                if (distances[i][j] < 0) {
                    throw new GviComputationException("Distance matrix contains a negative entry at (" + i + "," + j + ")");
                }
            }
        }
    }
}
