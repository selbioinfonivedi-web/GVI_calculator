package org.gvi.core.util;

import org.gvi.core.model.NucleotideSequence;
import org.gvi.core.model.SequenceAlignment;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Detects whether an alignment pools several genetically distinct lineages.
 * <p>
 * <b>Why this exists.</b> A molecular clock rate is only meaningful when every sequence descends
 * from one continuously evolving lineage. When an alignment instead contains, say, two serotypes
 * introduced independently, root-to-tip regression measures the distance <em>between</em> those
 * serotypes rather than substitution accumulated over time -- and it does so with a high R2, because
 * the points really do sit tightly on a line. That is precisely the failure the date-randomization
 * test flags, and the standard remedy is to split by lineage and estimate a rate within each.
 * <p>
 * The remedy is hard to act on without knowing where the split is, which is what this reports:
 * single-linkage clustering over pairwise p-distance at a supplied cut-off. Single linkage is the
 * right choice here because the question is whether the set separates into groups with a gap between
 * them, not whether the groups are individually tight.
 * <p>
 * This is deliberately a <em>diagnostic</em>, not an automatic re-analysis. Choosing a lineage
 * boundary is a judgement about the organism -- the threshold that separates FMDV serotypes is not
 * the one that separates SARS-CoV-2 variants -- so the tool reports the partition it found and lets
 * the analyst decide, rather than silently estimating rates on groups it invented.
 */
public final class LineagePartitioner {

    /**
     * Default p-distance cut-off. Sequences closer than this join the same putative lineage.
     * 5% is where this project's own genetic-distance reference table places "different strain /
     * serotype", so it is the same boundary the report already uses when describing a distance --
     * not a new constant invented here.
     */
    public static final double DEFAULT_CUTOFF = 0.05;

    /** Fewer clusters than this means no meaningful partition was found. */
    private static final int MIN_CLUSTERS_TO_REPORT = 2;

    public record Cluster(int index, List<String> sequenceIds) {
    }

    public record Result(List<Cluster> clusters, double cutoff, boolean partitioned, String description) {
    }

    private LineagePartitioner() {
    }

    public static Result detect(SequenceAlignment alignment) {
        return detect(alignment, DEFAULT_CUTOFF);
    }

    public static Result detect(SequenceAlignment alignment, double cutoff) {
        List<NucleotideSequence> seqs = alignment.getSequences();
        int n = seqs.size();
        if (n < 3) {
            return new Result(List.of(), cutoff, false, "Too few sequences to look for lineage structure.");
        }

        // Single-linkage agglomeration: start with every sequence in its own group, then merge any
        // two groups holding a pair closer than the cut-off.
        int[] group = new int[n];
        for (int i = 0; i < n; i++) group[i] = i;

        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                if (group[i] == group[j]) continue;
                if (pDistance(seqs.get(i).getSequence(), seqs.get(j).getSequence()) < cutoff) {
                    int from = group[j];
                    int to = group[i];
                    for (int k = 0; k < n; k++) if (group[k] == from) group[k] = to;
                }
            }
        }

        Map<Integer, List<String>> byGroup = new LinkedHashMap<>();
        for (int i = 0; i < n; i++) {
            byGroup.computeIfAbsent(group[i], g -> new ArrayList<>()).add(seqs.get(i).getId());
        }

        List<Cluster> clusters = new ArrayList<>();
        int idx = 1;
        for (var e : byGroup.entrySet()) clusters.add(new Cluster(idx++, List.copyOf(e.getValue())));

        if (clusters.size() < MIN_CLUSTERS_TO_REPORT) {
            return new Result(clusters, cutoff, false, String.format(
                    "No lineage structure detected: all %d sequences fall within one group at a %.0f%% p-distance "
                            + "cut-off.", n, cutoff * 100));
        }

        StringBuilder sb = new StringBuilder(String.format(
                "Lineage structure detected: the %d sequences separate into %d groups at a %.0f%% p-distance cut-off, "
                        + "which is where this tool's own distance table places \"different strain / serotype\". A single "
                        + "molecular-clock rate across groups measures the distance BETWEEN them, not substitution "
                        + "accumulated over time. Estimate a rate within each group separately: ",
                n, clusters.size(), cutoff * 100));
        for (Cluster c : clusters) {
            sb.append(String.format("[group %d, %d seq: %s] ", c.index(), c.sequenceIds().size(),
                    String.join(", ", c.sequenceIds())));
        }
        return new Result(clusters, cutoff, true, sb.toString().trim());
    }

    /** Pairwise-deletion proportion of differing comparable sites. */
    private static double pDistance(String a, String b) {
        int len = Math.min(a.length(), b.length());
        long comparable = 0;
        long diffs = 0;
        for (int i = 0; i < len; i++) {
            char x = a.charAt(i);
            char y = b.charAt(i);
            if (!IupacUtil.isComparable(x) || !IupacUtil.isComparable(y)) continue;
            comparable++;
            if (x != y) diffs++;
        }
        return comparable == 0 ? 0.0 : (double) diffs / comparable;
    }
}
