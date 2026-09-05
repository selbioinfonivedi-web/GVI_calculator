package org.gvi.algorithms.re;

import org.gvi.core.exception.GviComputationException;
import org.gvi.core.model.NucleotideSequence;
import org.gvi.core.model.SequenceAlignment;
import org.junit.jupiter.api.Test;

import java.time.LocalDate;
import java.util.ArrayList;
import java.util.List;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.assertThatThrownBy;
import static org.assertj.core.api.Assertions.within;

class PhylodynamicReEstimatorTest {

    /**
     * Direct tests of the Euler-Lotka (Wallinga &amp; Lipsitch 2007)
     * growth-rate-to-Re conversion, independent of the full tree-building
     * pipeline -- verifies the formula itself against its own defining
     * property (Re * sum(w(u)*e^-ru) == 1) and known boundary behavior,
     * rather than only checking the end-to-end result is directionally
     * plausible (which {@link #detectsPositiveGrowthRateFromAGenuinelyBranchingOutbreak}
     * already does).
     */
    @Test
    void zeroGrowthRateGivesReExactlyOne() {
        SerialInterval si = new SerialInterval(5.0, 2.0, 20);
        assertThat(PhylodynamicReEstimator.reFromGrowthRate(0.0, si)).isCloseTo(1.0, within(1e-12));
    }

    @Test
    void positiveGrowthRateGivesReGreaterThanOneAndNegativeGivesLessThanOne() {
        SerialInterval si = new SerialInterval(5.0, 2.0, 20);
        assertThat(PhylodynamicReEstimator.reFromGrowthRate(50.0, si)).isGreaterThan(1.0); // strong growth, per-year units
        assertThat(PhylodynamicReEstimator.reFromGrowthRate(-50.0, si)).isLessThan(1.0);
    }

    @Test
    void satisfiesItsOwnDefiningRenewalEquationForArbitraryGrowthRates() {
        SerialInterval si = new SerialInterval(5.0, 2.0, 20);
        for (double rPerYear : new double[]{-100.0, -10.0, 5.0, 30.0, 200.0}) {
            double re = PhylodynamicReEstimator.reFromGrowthRate(rPerYear, si);
            double rPerDay = rPerYear / 365.25;
            double sum = 0.0;
            for (int u = 1; u <= si.tMax(); u++) sum += si.weight(u) * Math.exp(-rPerDay * u);
            assertThat(re * sum).isCloseTo(1.0, within(1e-9));
        }
    }

    @Test
    void reducesToTheOldLinearApproximationForSmallGrowthRates() {
        SerialInterval si = new SerialInterval(5.0, 2.0, 20);
        double meanDays = 5.0;
        double rPerYearSmall = 0.001; // small enough that the first-order Taylor expansion should dominate
        double exact = PhylodynamicReEstimator.reFromGrowthRate(rPerYearSmall, si);
        double linearApprox = 1.0 + rPerYearSmall * (meanDays / 365.25);
        assertThat(exact).isCloseTo(linearApprox, within(1e-6));
    }

    private final PhylodynamicReEstimator estimator = new PhylodynamicReEstimator();

    /**
     * Simulates a genuinely branching outbreak: a lineage that splits in two
     * every generation (1 -> 2 -> 4 -> 8), with samples taken at every
     * generation (not just the final one, mirroring real surveillance that
     * samples throughout an outbreak's growth). Each node gets a unique
     * "private mutation" position, and its sequence carries its own private
     * mutation plus every ancestor's along its path from the root -- so the
     * resulting Neighbor-Joining tree genuinely reflects this branching
     * structure instead of being inferred from arbitrary data.
     * <p>
     * This is a real exponentially-doubling lineage count, so the
     * lineages-through-time regression should detect clear positive growth
     * (Re > 1) with a reasonably strong signal -- checked directionally
     * (not against a hand-solved exact value, since the full pipeline's
     * closed-form answer isn't practical to re-derive independently) the
     * same way the underlying algorithm's own literature standing is
     * validated: by simulation, not by trusting a single point value.
     */
    @Test
    void detectsPositiveGrowthRateFromAGenuinelyBranchingOutbreak() {
        int length = 2000;
        char[] base = new char[length];
        java.util.Arrays.fill(base, 'A');
        LocalDate rootDate = LocalDate.of(2020, 1, 1);

        List<NucleotideSequence> seqs = new ArrayList<>();
        seqs.add(new NucleotideSequence("ref", new String(base)).withMetadata(rootDate, null, null));

        // gen1: 2 nodes, each with 1 private mutation, dated root+1yr
        int[][] gen1Positions = {{0}, {1}};
        List<int[]> gen1Muts = new ArrayList<>();
        int idx = 0;
        for (int i = 0; i < 2; i++) {
            int[] muts = {gen1Positions[i][0]};
            gen1Muts.add(muts);
            seqs.add(taxon("gen1_" + i, base, muts, rootDate.plusYears(1)));
        }

        // gen2: 4 nodes (2 children per gen1 parent), each inherits parent's mutation + 1 new private one
        int nextPos = 2;
        List<int[]> gen2Muts = new ArrayList<>();
        for (int parent = 0; parent < 2; parent++) {
            for (int child = 0; child < 2; child++) {
                int[] muts = append(gen1Muts.get(parent), nextPos++);
                gen2Muts.add(muts);
                seqs.add(taxon("gen2_" + parent + "_" + child, base, muts, rootDate.plusYears(2)));
            }
        }

        // gen3: 8 nodes (2 children per gen2 parent), same pattern
        for (int parent = 0; parent < 4; parent++) {
            for (int child = 0; child < 2; child++) {
                int[] muts = append(gen2Muts.get(parent), nextPos++);
                seqs.add(taxon("gen3_" + parent + "_" + child, base, muts, rootDate.plusYears(3)));
            }
        }

        SequenceAlignment aln = SequenceAlignment.of(seqs, "ref");
        ReResult r = estimator.compute(aln, 5.0);

        assertThat(r.method()).isEqualTo(ReMethod.PHYLODYNAMIC_FALLBACK);
        assertThat(r.reMean()).isGreaterThan(1.0);
        assertThat(r.doublingTimeDays()).isNotNull();
        assertThat(r.doublingTimeDays()).isPositive();
        assertThat(r.diagnostics()).anyMatch(s -> s.contains("lineages-through-time") || s.contains("Pybus"));
    }

    private static int[] append(int[] existing, int extra) {
        int[] result = java.util.Arrays.copyOf(existing, existing.length + 1);
        result[existing.length] = extra;
        return result;
    }

    private static NucleotideSequence taxon(String id, char[] base, int[] mutationPositions, LocalDate date) {
        char[] copy = base.clone();
        for (int pos : mutationPositions) copy[pos] = 'G';
        return new NucleotideSequence(id, new String(copy)).withMetadata(date, null, null);
    }

    @Test
    void throwsWhenTooFewDatedSequences() {
        List<NucleotideSequence> seqs = List.of(
                new NucleotideSequence("s1", "ACGT").withMetadata(LocalDate.of(2024, 1, 1), null, null),
                new NucleotideSequence("s2", "ACGT").withMetadata(LocalDate.of(2024, 1, 2), null, null)
        );
        SequenceAlignment aln = SequenceAlignment.of(seqs);
        assertThatThrownBy(() -> estimator.compute(aln, 5.0)).isInstanceOf(GviComputationException.class);
    }
}
