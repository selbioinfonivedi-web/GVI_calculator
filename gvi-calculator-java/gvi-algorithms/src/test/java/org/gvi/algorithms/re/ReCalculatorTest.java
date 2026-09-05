package org.gvi.algorithms.re;

import org.gvi.core.exception.GviInputException;
import org.gvi.core.model.IncidencePoint;
import org.gvi.core.model.NucleotideSequence;
import org.gvi.core.model.SequenceAlignment;
import org.junit.jupiter.api.Test;

import java.time.LocalDate;
import java.util.ArrayList;
import java.util.List;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.assertThatThrownBy;

class ReCalculatorTest {

    private final ReCalculator calculator = new ReCalculator();

    @Test
    void prefersIncidenceBasedEstimatorWhenAvailable() {
        SerialInterval si = SerialInterval.defaultProfile();
        LocalDate start = LocalDate.of(2024, 1, 1);
        double[] incidence = new double[80];
        for (int i = 0; i < si.tMax(); i++) incidence[i] = 100.0;
        for (int t = si.tMax(); t < incidence.length; t++) {
            double lambda = 0;
            for (int u = 1; u <= si.tMax(); u++) lambda += incidence[t - u] * si.weight(u);
            incidence[t] = 1.2 * lambda;
        }
        List<IncidencePoint> points = new ArrayList<>();
        for (int i = 0; i < incidence.length; i++) points.add(new IncidencePoint(start.plusDays(i), incidence[i]));

        ReResult r = calculator.compute(points, null);
        assertThat(r.method()).isEqualTo(ReMethod.CORI_INCIDENCE);
    }

    @Test
    void fallsBackToPhylodynamicWhenOnlyAlignmentAvailable() {
        ReResult r = calculator.compute(List.of(), buildBranchingAlignment());
        assertThat(r.method()).isEqualTo(ReMethod.PHYLODYNAMIC_FALLBACK);
    }

    @Test
    void throwsWhenNeitherIncidenceNorAlignmentAvailable() {
        assertThatThrownBy(() -> calculator.compute(List.of(), null)).isInstanceOf(GviInputException.class);
    }

    @Test
    void fallsBackWhenIncidenceSeriesTooShort() {
        List<IncidencePoint> tooShort = List.of(new IncidencePoint(LocalDate.of(2024, 1, 1), 10));
        ReResult r = calculator.compute(tooShort, buildBranchingAlignment());
        assertThat(r.method()).isEqualTo(ReMethod.PHYLODYNAMIC_FALLBACK);
        assertThat(r.diagnostics().get(0)).contains("Incidence-based estimation failed");
    }

    /**
     * A genuinely branching (1->2->4->8), genuinely mutating outbreak
     * simulation -- the tree-based phylodynamic estimator needs real
     * genetic + temporal signal to time-scale a tree, unlike the old
     * sample-count-only heuristic it replaced. See
     * PhylodynamicReEstimatorTest for the full annotated version of this
     * construction.
     */
    private static SequenceAlignment buildBranchingAlignment() {
        int length = 2000;
        char[] base = new char[length];
        java.util.Arrays.fill(base, 'A');
        LocalDate rootDate = LocalDate.of(2020, 1, 1);

        List<NucleotideSequence> seqs = new ArrayList<>();
        seqs.add(new NucleotideSequence("ref", new String(base)).withMetadata(rootDate, null, null));

        List<int[]> gen1Muts = new ArrayList<>();
        for (int i = 0; i < 2; i++) {
            int[] muts = {i};
            gen1Muts.add(muts);
            seqs.add(taxon("gen1_" + i, base, muts, rootDate.plusYears(1)));
        }

        int nextPos = 2;
        List<int[]> gen2Muts = new ArrayList<>();
        for (int parent = 0; parent < 2; parent++) {
            for (int child = 0; child < 2; child++) {
                int[] muts = append(gen1Muts.get(parent), nextPos++);
                gen2Muts.add(muts);
                seqs.add(taxon("gen2_" + parent + "_" + child, base, muts, rootDate.plusYears(2)));
            }
        }

        for (int parent = 0; parent < 4; parent++) {
            for (int child = 0; child < 2; child++) {
                int[] muts = append(gen2Muts.get(parent), nextPos++);
                seqs.add(taxon("gen3_" + parent + "_" + child, base, muts, rootDate.plusYears(3)));
            }
        }

        return SequenceAlignment.of(seqs, "ref");
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
}
