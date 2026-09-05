package org.gvi.cli;

import org.junit.jupiter.api.Test;

import java.util.List;

import static org.assertj.core.api.Assertions.assertThat;

class ExclusionExplainerTest {

    /** Verbatim shapes the pipeline actually emits, so the mapping is pinned to real strings. */
    private static final String GAPS =
            "MB (excluded from composite GVI): the alignment is 38% gap characters and the mean mutation burden is "
                    + "21% of its length, so the sequences are not the closely-related, genuinely-aligned set...";
    private static final String GAPS_GD = GAPS.replace("MB (", "GD (");
    private static final String GAPS_PI = GAPS.replace("MB (", "pi (");
    private static final String DRT =
            "mu (excluded from composite GVI): the date-randomization test failed -- randomly reshuffled collection "
                    + "dates reproduce this rate...";
    private static final String ZTEST =
            "dN/dS (excluded from composite GVI): the pooled dN/dS is 2.536 but the Nei-Gojobori Z-test cannot "
                    + "distinguish it from neutrality...";
    private static final String LTT =
            "Re (excluded from composite GVI): Re came from the lineages-through-time fallback, which on real data "
                    + "returns ~1.0 for every pathogen...";

    @Test
    void groupsOneUnderlyingProblemAsOneProblem() {
        // Five indices, one bad alignment. The user has one thing to fix, not five.
        var groups = ExclusionExplainer.explain(List.of(GAPS, GAPS_GD, GAPS_PI));
        assertThat(groups).hasSize(1);
        assertThat(groups.get(0).affectedIndices()).containsExactly("MB", "GD", "pi");
    }

    @Test
    void separatesDistinctCauses() {
        var groups = ExclusionExplainer.explain(List.of(GAPS, DRT, ZTEST, LTT));
        assertThat(groups).hasSize(4);
        assertThat(groups).extracting(ExclusionExplainer.Group::cause)
                .doesNotHaveDuplicates();
    }

    @Test
    void ignoresNonExclusionSkippedEntries() {
        var groups = ExclusionExplainer.explain(List.of(
                "CAI: no --codon-usage table or --codon-usage-species supplied",
                "RI: PHI test needs at least 4 informative biallelic sites, found 0.",
                GAPS));
        assertThat(groups).hasSize(1);
        assertThat(groups.get(0).affectedIndices()).containsExactly("MB");
    }

    @Test
    void everyRecognisedCauseCarriesAConcreteRemedy() {
        for (String msg : List.of(GAPS, DRT, ZTEST, LTT)) {
            String remedy = ExclusionExplainer.remedyFor(ExclusionExplainer.causeOf(msg));
            assertThat(remedy).isNotBlank();
            assertThat(remedy).as("remedy for %s should be actionable, not a pointer", msg)
                    .doesNotContain("See the full report");
        }
    }

    @Test
    void theAlignmentRemedyWarnsAgainstTheWrongInstinct() {
        // Re-running MAFFT is the instinctive fix and it does not work -- the aligner already placed
        // these correctly, they simply have different coverage windows.
        String remedy = ExclusionExplainer.remedyFor(ExclusionExplainer.causeOf(GAPS));
        assertThat(remedy).contains("MAFFT").contains("will not help");
    }

    @Test
    void unrecognisedMessageFallsBackWithoutThrowing() {
        var groups = ExclusionExplainer.explain(List.of("XX (excluded from composite GVI): something new"));
        assertThat(groups).hasSize(1);
        assertThat(groups.get(0).cause()).isEqualTo("Other quality-gate failure");
    }

    @Test
    void handlesEmptyAndNullSafely() {
        assertThat(ExclusionExplainer.explain(List.of())).isEmpty();
        assertThat(ExclusionExplainer.causeOf(null)).isEqualTo("Other quality-gate failure");
    }
}
