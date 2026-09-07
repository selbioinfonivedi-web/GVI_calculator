package org.gvi.core.util;

import org.gvi.core.model.NucleotideSequence;
import org.gvi.core.model.SequenceAlignment;

import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

/**
 * Two cheap checks on an alignment, run before anything is computed from it.
 * <p>
 * Both exist because the corpus contained datasets that ran to completion and produced a
 * confident-looking GVI from inputs that could not support one. Neither failure was visible until
 * the indices had been computed and gated, by which point the diagnosis lived scattered across
 * per-index exclusion messages. Stating it once, up front, in the alignment's own terms is a
 * better place to learn it.
 *
 * <h2>What the two checks caught</h2>
 * <table>
 *   <caption>Measured on the project corpus</caption>
 *   <tr><th>Dataset</th><th>Identity floor</th><th>Fully-covered columns</th></tr>
 *   <tr><td>foot_and_mouth_disease</td><td>87%</td><td>100%</td></tr>
 *   <tr><td>enterotoxaemia</td><td>49%</td><td>47%</td></tr>
 *   <tr><td>haemorrhagic_septicaemia</td><td>46%</td><td>43%</td></tr>
 * </table>
 * <p>
 * The two broken sets separate cleanly from the healthy one on both axes, which is what makes a
 * threshold worth having at all.
 * <p>
 * These are <em>warnings</em>, not gates. The pipeline already gates individual indices on their own
 * evidence; this only makes the shape of the input legible before that happens. A user analysing a
 * deliberately divergent panel should not be blocked, only told what they are holding.
 */
public final class AlignmentPreflight {

    /**
     * Pairwise identity below which two sequences are unlikely to be the same locus from the same
     * organism group. Sequences this far apart are usually a pooled alignment: several unrelated
     * groups concatenated into one file, where every position-wise index then compares things that
     * are not homologous.
     * <p>
     * 0.70 sits well below the healthy corpus dataset (0.87 floor) and well above both broken ones
     * (0.49 and 0.46), so it separates the observed cases without sitting on top of either.
     */
    public static final double MIN_PAIRWISE_IDENTITY = 0.70;

    /**
     * Fraction of columns every sequence must actually cover. A low value means the records span
     * different windows of the locus, so position-wise indices are comparing observed bases against
     * absent data rather than against differences.
     * <p>
     * This is the condition {@code --trim-to-covered} exists to fix, so the finding names it.
     */
    public static final double MIN_FULLY_COVERED_FRACTION = 0.60;

    private static final String BASES = "ACGT";

    private AlignmentPreflight() {
    }

    /** What the checks found. {@link #findings()} is empty when the alignment looks usable. */
    public record Report(int sequences, int length, double minPairwiseIdentity,
                         double fullyCoveredFraction, int fullyCoveredColumns,
                         List<String> findings) {

        public boolean clean() {
            return findings.isEmpty();
        }
    }

    public static Report check(SequenceAlignment alignment) {
        List<NucleotideSequence> seqs = alignment.getSequences();
        int n = seqs.size();
        int length = alignment.length();
        List<String> findings = new ArrayList<>();

        double minIdentity = minPairwiseIdentity(seqs);
        int covered = fullyCoveredColumns(seqs, length);
        double coveredFraction = length == 0 ? 0.0 : covered / (double) length;

        if (n >= 2 && minIdentity < MIN_PAIRWISE_IDENTITY) {
            findings.add(String.format(Locale.ROOT,
                    "Pairwise identity falls to %.0f%% between the two most divergent sequences, below the %.0f%% "
                            + "expected for one locus from one organism group. This usually means the file pools "
                            + "several unrelated groups, in which case every position-wise index (pi, MB, GD, dN/dS) "
                            + "is comparing positions that are not homologous. Split the file by group and score each "
                            + "separately, or confirm the divergence is real for this locus.",
                    minIdentity * 100, MIN_PAIRWISE_IDENTITY * 100));
        }

        if (length > 0 && coveredFraction < MIN_FULLY_COVERED_FRACTION) {
            findings.add(String.format(Locale.ROOT,
                    "Only %d of %d columns (%.0f%%) are covered by every sequence, below %.0f%%. The records span "
                            + "different windows of the locus, so position-wise indices compare observed bases "
                            + "against absent data. Use --trim-to-covered to restrict the alignment to the shared "
                            + "window, and check how much of it survives.",
                    covered, length, coveredFraction * 100, MIN_FULLY_COVERED_FRACTION * 100));
        }

        return new Report(n, length, minIdentity, coveredFraction, covered, List.copyOf(findings));
    }

    /** Lowest identity over all pairs, counting only positions where both sequences have a base. */
    private static double minPairwiseIdentity(List<NucleotideSequence> seqs) {
        double min = 1.0;
        for (int i = 0; i < seqs.size(); i++) {
            String a = seqs.get(i).getSequence();
            for (int j = i + 1; j < seqs.size(); j++) {
                String b = seqs.get(j).getSequence();
                int comparable = 0, same = 0;
                int len = Math.min(a.length(), b.length());
                for (int k = 0; k < len; k++) {
                    char x = a.charAt(k), y = b.charAt(k);
                    if (BASES.indexOf(x) >= 0 && BASES.indexOf(y) >= 0) {
                        comparable++;
                        if (x == y) same++;
                    }
                }
                // A pair with no comparable position at all says nothing about identity; excluding
                // it stops disjoint coverage windows being reported as 0% identity, which is a
                // coverage finding rather than a divergence one and is reported separately below.
                if (comparable > 0) {
                    min = Math.min(min, same / (double) comparable);
                }
            }
        }
        return min;
    }

    private static int fullyCoveredColumns(List<NucleotideSequence> seqs, int length) {
        int covered = 0;
        outer:
        for (int j = 0; j < length; j++) {
            for (NucleotideSequence s : seqs) {
                String v = s.getSequence();
                if (j >= v.length() || BASES.indexOf(v.charAt(j)) < 0) {
                    continue outer;
                }
            }
            covered++;
        }
        return covered;
    }
}
