package org.gvi.cli;

import org.gvi.algorithms.dnds.NeiGojoboriSelectionTest;
import org.gvi.composite.IndexKey;
import org.gvi.core.model.NucleotideSequence;
import org.gvi.core.model.SequenceAlignment;
import org.gvi.core.spi.IndexResult;
import org.gvi.core.spi.SimpleIndexResult;
import org.junit.jupiter.api.Test;

import java.util.EnumMap;
import java.util.List;
import java.util.Map;

import static org.assertj.core.api.Assertions.assertThat;

/**
 * dN/dS above 1 only means positive selection when it is statistically distinguishable from 1.
 * These pin the gate that decides whether such a ratio is allowed to carry composite weight.
 */
class SelectionSignificanceGateTest {

    /** A clean, well-aligned 4-sequence alignment, so only the selection gate can fire. */
    private static SequenceAlignment cleanAlignment() {
        return SequenceAlignment.of(List.of(
                new NucleotideSequence("ref", "ATGCATGCATGCATGCATGCATGCATGCATGC"),
                new NucleotideSequence("s1", "ATGCATGCATGCATGCATGCATGCATGCATGA"),
                new NucleotideSequence("s2", "ATGCATGCATGCATGCATGCATGCATGCATGT"),
                new NucleotideSequence("s3", "ATGCATGCATGCATGCATGCATGCATGCATGG")), null);
    }

    private static Map<IndexKey, IndexResult> withDnDs(double omega) {
        Map<IndexKey, IndexResult> m = new EnumMap<>(IndexKey.class);
        m.put(IndexKey.DNDS, new SimpleIndexResult("dN/dS", omega, "category", List.of()));
        return m;
    }

    private static NeiGojoboriSelectionTest.Result verdict(NeiGojoboriSelectionTest.Verdict v) {
        return new NeiGojoboriSelectionTest.Result(1.1, 0.27, v, "test fixture");
    }

    @Test
    void omegaAboveOneWithoutSignificanceIsExcluded() {
        // The real Trypanosomiasis case: omega = 2.536, not distinguishable from neutrality.
        var gate = QualityGate.evaluate(cleanAlignment(), Map.of(),
                withDnDs(2.536), verdict(NeiGojoboriSelectionTest.Verdict.NOT_SIGNIFICANT));

        assertThat(gate.isExcluded(IndexKey.DNDS)).isTrue();
        assertThat(gate.exclusions()).singleElement().satisfies(e ->
                assertThat(e.reason()).contains("cannot distinguish it from neutrality").contains("frameshift"));
    }

    @Test
    void omegaAboveOneWithSignificanceIsKept() {
        var gate = QualityGate.evaluate(cleanAlignment(), Map.of(),
                withDnDs(2.536), verdict(NeiGojoboriSelectionTest.Verdict.POSITIVE_SELECTION));

        assertThat(gate.isExcluded(IndexKey.DNDS)).isFalse();
    }

    @Test
    void omegaBelowOneIsNeverGatedOnSignificance() {
        // Purifying selection, or simply no detectable selection, does not overstate anything --
        // gating it would discard the common, unremarkable case.
        for (var v : NeiGojoboriSelectionTest.Verdict.values()) {
            var gate = QualityGate.evaluate(cleanAlignment(), Map.of(), withDnDs(0.08), verdict(v));
            assertThat(gate.isExcluded(IndexKey.DNDS)).as("omega 0.08 with verdict %s", v).isFalse();
        }
    }

    @Test
    void alignmentIntegrityFailureStillTakesPrecedence() {
        // A gap-riddled alignment invalidates dN/dS regardless of the selection verdict, and the
        // reason reported should be the alignment, not the Z-test.
        SequenceAlignment gappy = SequenceAlignment.of(List.of(
                new NucleotideSequence("ref", "ATGC--------ATGCATGC--------ATGC"),
                new NucleotideSequence("s1", "ATGC--------ATGCATGC--------ATGA")), null);

        var gate = QualityGate.evaluate(gappy, Map.of(),
                withDnDs(2.536), verdict(NeiGojoboriSelectionTest.Verdict.POSITIVE_SELECTION));

        assertThat(gate.isExcluded(IndexKey.DNDS)).isTrue();
        assertThat(gate.exclusions()).anySatisfy(e -> {
            if (e.key() == IndexKey.DNDS) assertThat(e.reason()).contains("gap characters");
        });
    }
}
