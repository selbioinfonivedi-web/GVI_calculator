package org.gvi.core.util;

import org.gvi.core.model.NucleotideSequence;
import org.gvi.core.model.SequenceAlignment;
import org.junit.jupiter.api.Test;

import java.util.List;
import java.util.Random;

import static org.assertj.core.api.Assertions.assertThat;

/**
 * The thresholds here were set from the project corpus, where two datasets ran to completion and
 * produced a confident-looking GVI from inputs that could not support one. These tests pin the
 * separation that made the thresholds defensible: the healthy dataset sits at 87% identity and
 * 100% coverage, the two broken ones at 46-49% and 43-47%.
 */
class AlignmentPreflightTest {

    private static final char[] BASES = {'A', 'C', 'G', 'T'};

    @Test
    void aCleanAlignmentProducesNoFindings() {
        SequenceAlignment alignment = alignmentOf(
                "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG",
                "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG",
                "ATGGCCATTGTAATGGCCCGCTGAAAGGGTGCCCGATAC",
                "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCACGATAG");

        AlignmentPreflight.Report report = AlignmentPreflight.check(alignment);

        assertThat(report.clean()).as("findings: %s", report.findings()).isTrue();
        assertThat(report.fullyCoveredFraction()).isEqualTo(1.0);
        assertThat(report.minPairwiseIdentity()).isGreaterThan(AlignmentPreflight.MIN_PAIRWISE_IDENTITY);
    }

    /**
     * Two unrelated groups concatenated into one file. This is enterotoxaemia's actual shape:
     * pairwise identity spans 49% to 100%, because within a group the sequences match and across
     * groups they do not. Every position-wise index then compares positions that are not homologous.
     */
    @Test
    void pooledUnrelatedGroupsAreReported() {
        Random rng = new Random(42);
        String groupA = randomSequence(600, rng);
        String groupB = randomSequence(600, rng); // independent, so ~25% identity to groupA

        SequenceAlignment alignment = alignmentOf(groupA, mutate(groupA, 5, rng),
                                                  groupB, mutate(groupB, 5, rng));

        AlignmentPreflight.Report report = AlignmentPreflight.check(alignment);

        assertThat(report.minPairwiseIdentity()).isLessThan(AlignmentPreflight.MIN_PAIRWISE_IDENTITY);
        assertThat(report.findings()).anyMatch(f -> f.contains("Pairwise identity falls to"));
        assertThat(report.findings()).anyMatch(f -> f.contains("pools"));
    }

    /**
     * Records spanning different windows of the locus. This is haemorrhagic septicaemia's shape:
     * 43% of columns covered by every sequence, so the rest is absent data being compared against
     * observed bases.
     */
    @Test
    void differingCoverageWindowsAreReported() {
        // One sequence covers the first half, one the second; the overlap is small.
        String full = "ACGT".repeat(100);
        String firstHalf = full.substring(0, 220) + "-".repeat(180);
        String secondHalf = "-".repeat(180) + full.substring(180);

        SequenceAlignment alignment = alignmentOf(full, firstHalf, secondHalf);

        AlignmentPreflight.Report report = AlignmentPreflight.check(alignment);

        assertThat(report.fullyCoveredFraction()).isLessThan(AlignmentPreflight.MIN_FULLY_COVERED_FRACTION);
        assertThat(report.findings()).anyMatch(f -> f.contains("covered by every sequence"));
        assertThat(report.findings()).anyMatch(f -> f.contains("--trim-to-covered"));
    }

    /**
     * A pair sharing no comparable position says nothing about divergence. Counting it as 0%
     * identity would report a coverage problem as a divergence one; the coverage check reports it.
     */
    @Test
    void disjointCoverageIsNotReportedAsZeroIdentity() {
        String a = "ACGTACGTAC" + "-".repeat(10);
        String b = "-".repeat(10) + "ACGTACGTAC";

        AlignmentPreflight.Report report = AlignmentPreflight.check(alignmentOf(a, b));

        assertThat(report.minPairwiseIdentity())
                .as("no shared position means no identity evidence, not zero identity")
                .isEqualTo(1.0);
        assertThat(report.findings()).noneMatch(f -> f.contains("Pairwise identity falls to"));
        assertThat(report.findings()).anyMatch(f -> f.contains("covered by every sequence"));
    }

    @Test
    void aSingleSequenceCannotBeDivergentFromAnything() {
        AlignmentPreflight.Report report = AlignmentPreflight.check(alignmentOf("ACGTACGTAC"));
        assertThat(report.findings()).noneMatch(f -> f.contains("Pairwise identity falls to"));
    }

    // ── helpers ──────────────────────────────────────────────────────────────
    private static SequenceAlignment alignmentOf(String... seqs) {
        List<NucleotideSequence> list = new java.util.ArrayList<>();
        for (int i = 0; i < seqs.length; i++) {
            list.add(new NucleotideSequence("s" + i, seqs[i]));
        }
        return SequenceAlignment.of(list, "s0");
    }

    private static String randomSequence(int length, Random rng) {
        StringBuilder sb = new StringBuilder(length);
        for (int i = 0; i < length; i++) sb.append(BASES[rng.nextInt(4)]);
        return sb.toString();
    }

    private static String mutate(String seq, int changes, Random rng) {
        char[] c = seq.toCharArray();
        for (int i = 0; i < changes; i++) c[rng.nextInt(c.length)] = BASES[rng.nextInt(4)];
        return new String(c);
    }
}
