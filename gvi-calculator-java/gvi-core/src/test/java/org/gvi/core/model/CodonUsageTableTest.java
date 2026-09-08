package org.gvi.core.model;

import org.gvi.core.exception.GviInputException;
import org.gvi.core.util.GeneticCode;
import org.junit.jupiter.api.Test;

import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.assertThatThrownBy;
import static org.assertj.core.api.Assertions.within;

/**
 * CAI is the geometric mean of per-codon weights from this table (Sharp &amp; Li 1987), so every
 * property here multiplies straight into the reported index. It was at 0% coverage.
 * <p>
 * The weight is w(codon) = f(codon) / f(most frequent synonymous codon), which makes two things
 * load-bearing and easy to get wrong: the maximum is taken <em>within an amino acid's family</em>
 * and not across the table, and a zero weight is floored, because a geometric mean takes logs and
 * log(0) is -Infinity -- one never-observed codon would otherwise drive the whole CAI to zero.
 */
class CodonUsageTableTest {

    /** Frequencies for one family; every other codon defaults to absent. */
    private static Map<String, Double> freqs(Object... pairs) {
        Map<String, Double> m = new LinkedHashMap<>();
        for (int i = 0; i < pairs.length; i += 2) {
            m.put((String) pairs[i], ((Number) pairs[i + 1]).doubleValue());
        }
        return m;
    }

    @Test
    void theMostFrequentCodonInAFamilyGetsWeightOne() {
        // Phe: TTT and TTC. TTC is three times as common.
        CodonUsageTable t = CodonUsageTable.fromFrequencies("host", freqs("TTT", 10.0, "TTC", 30.0));

        assertThat(t.weightOf("TTC")).isEqualTo(1.0);
        assertThat(t.weightOf("TTT")).isCloseTo(10.0 / 30.0, within(1e-12));
    }

    /**
     * The maximum is per amino acid. If it were taken across the whole table, a codon that is common
     * within a rare family would be scored as though it were rare, and CAI would fall for sequences
     * that are in fact well adapted.
     */
    @Test
    void theMaximumIsTakenWithinEachFamilyNotAcrossTheWholeTable() {
        // Phe is rare overall; Gly is common. Each family's own best codon must still score 1.0.
        CodonUsageTable t = CodonUsageTable.fromFrequencies("host", freqs(
                "TTT", 1.0, "TTC", 2.0,
                "GGT", 50.0, "GGC", 100.0, "GGA", 25.0, "GGG", 25.0));

        assertThat(t.weightOf("TTC")).as("best Phe codon").isEqualTo(1.0);
        assertThat(t.weightOf("GGC")).as("best Gly codon").isEqualTo(1.0);
        assertThat(t.weightOf("TTT")).isCloseTo(0.5, within(1e-12));
        assertThat(t.weightOf("GGT")).isCloseTo(0.5, within(1e-12));
    }

    /**
     * CAI is a geometric mean, so it takes the log of every weight. An exact zero would make that
     * -Infinity and collapse the whole index on a single unobserved codon.
     */
    @Test
    void anUnobservedCodonIsFlooredRatherThanZeroedSoTheGeometricMeanSurvives() {
        CodonUsageTable t = CodonUsageTable.fromFrequencies("host", freqs("TTT", 0.0, "TTC", 30.0));

        assertThat(t.weightOf("TTT")).isGreaterThan(0.0);
        assertThat(Math.log(t.weightOf("TTT"))).isFinite();
        assertThat(t.weightOf("TTT")).isCloseTo(0.0001, within(1e-12));
    }

    @Test
    void aFamilyWithNoObservationsAtAllIsFlooredRatherThanNaN() {
        // Only Phe supplied; every other family has max = 0, which would be 0/0.
        CodonUsageTable t = CodonUsageTable.fromFrequencies("host", freqs("TTT", 10.0, "TTC", 30.0));

        for (String gly : GeneticCode.SYNONYMOUS_CODONS.get('G')) {
            assertThat(t.weightOf(gly)).as("weight of %s", gly).isNotNaN().isGreaterThan(0.0);
        }
    }

    /** Stops are excluded from CAI by convention, so they are not in the table at all. */
    @Test
    void stopCodonsAreNotWeighted() {
        CodonUsageTable t = CodonUsageTable.fromFrequencies("host", freqs("TTT", 10.0, "TAA", 5.0));

        for (String stop : new String[]{"TAA", "TAG", "TGA"}) {
            assertThatThrownBy(() -> t.weightOf(stop))
                    .as("stop codon %s", stop)
                    .isInstanceOf(GviInputException.class);
        }
    }

    /** Every sense codon must have a weight, or CAI would throw partway through a real sequence. */
    @Test
    void everySenseCodonIsWeighted() {
        CodonUsageTable t = CodonUsageTable.fromFrequencies("host", freqs("TTT", 1.0));

        int weighted = 0;
        for (Map.Entry<String, Character> e : GeneticCode.CODON_TABLE.entrySet()) {
            if (e.getValue() == GeneticCode.STOP) continue;
            assertThat(t.weightOf(e.getKey())).as("weight of %s", e.getKey()).isBetween(0.0, 1.0);
            weighted++;
        }
        assertThat(weighted).isEqualTo(61);
    }

    @Test
    void codonLookupIsCaseInsensitive() {
        CodonUsageTable t = CodonUsageTable.fromFrequencies("host", freqs("TTT", 10.0, "TTC", 30.0));
        assertThat(t.weightOf("ttc")).isEqualTo(t.weightOf("TTC"));
    }

    @Test
    void anUnknownCodonIsRejectedNamingTheTable() {
        CodonUsageTable t = CodonUsageTable.fromFrequencies("e_coli", freqs("TTT", 1.0));

        assertThatThrownBy(() -> t.weightOf("XYZ"))
                .isInstanceOf(GviInputException.class)
                .hasMessageContaining("XYZ")
                .hasMessageContaining("e_coli");
    }

    @Test
    void anEmptyOrNullTableIsRejectedRatherThanProducingAllFlooredWeights() {
        assertThatThrownBy(() -> CodonUsageTable.fromFrequencies("host", new HashMap<>()))
                .isInstanceOf(GviInputException.class)
                .hasMessageContaining("is empty");
        assertThatThrownBy(() -> CodonUsageTable.fromFrequencies("host", null))
                .isInstanceOf(GviInputException.class)
                .hasMessageContaining("is empty");
    }

    @Test
    void theNameIsCarriedForReporting() {
        assertThat(CodonUsageTable.fromFrequencies("Bos taurus", freqs("TTT", 1.0)).getName())
                .isEqualTo("Bos taurus");
    }

    /**
     * Weights are ratios, so scaling every count by a constant -- counts versus per-thousand
     * frequencies, the two forms real tables ship in -- must not move CAI.
     */
    @Test
    void weightsAreInvariantToTheScaleOfTheInputCounts() {
        CodonUsageTable counts = CodonUsageTable.fromFrequencies("a",
                freqs("TTT", 10.0, "TTC", 30.0, "GGT", 5.0, "GGC", 20.0));
        CodonUsageTable perThousand = CodonUsageTable.fromFrequencies("b",
                freqs("TTT", 0.10, "TTC", 0.30, "GGT", 0.05, "GGC", 0.20));

        for (String c : new String[]{"TTT", "TTC", "GGT", "GGC"}) {
            assertThat(perThousand.weightOf(c)).as("weight of %s", c)
                    .isCloseTo(counts.weightOf(c), within(1e-12));
        }
    }
}
