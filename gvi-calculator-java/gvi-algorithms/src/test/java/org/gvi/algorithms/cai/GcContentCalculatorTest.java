package org.gvi.algorithms.cai;

import org.gvi.core.model.NucleotideSequence;
import org.junit.jupiter.api.Test;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.within;

class GcContentCalculatorTest {

    private final GcContentCalculator calc = new GcContentCalculator();

    @Test
    void computesExactGcPercentAndDeviation() {
        // 10 bases: 4 G/C -> 40% GC
        NucleotideSequence seq = new NucleotideSequence("s1", "GCGCAAAAAA");
        GcResult r = calc.compute(seq, 38.0);
        assertThat(r.observedGcPercent()).isCloseTo(40.0, within(1e-9));
        assertThat(r.deviationPercent()).isCloseTo(2.0, within(1e-9));
        assertThat(r.category()).contains("Minor drift");
    }

    @Test
    void classifiesLargeDeviationAsAnomaly() {
        NucleotideSequence seq = new NucleotideSequence("s1", "GGGGGGGGGG"); // 100% GC
        GcResult r = calc.compute(seq, 40.0);
        assertThat(r.deviationPercent()).isCloseTo(60.0, within(1e-9));
        assertThat(r.category()).contains("Anomaly");
    }

    @Test
    void excludesAmbiguousBasesFromCalculation() {
        // 8 real bases (4 GC) + 2 N's -> GC% computed over 8, not 10
        NucleotideSequence seq = new NucleotideSequence("s1", "GCGCAAAANN");
        GcResult r = calc.compute(seq, 50.0);
        assertThat(r.observedGcPercent()).isCloseTo(50.0, within(1e-9));
    }

    @Test
    void threeArgOverloadHasNoSignatureFieldWithoutAReference() {
        NucleotideSequence seq = new NucleotideSequence("s1", "GCGCAAAAAA");
        assertThat(calc.compute(seq, 38.0).mutationalSignature()).isNull();
    }

    @Test
    void flagsApobecLikeSignatureInCategoryAndDiagnosticsWhenReferenceSuppliedAndEnriched() {
        StringBuilder refB = new StringBuilder();
        StringBuilder qryB = new StringBuilder();
        for (int i = 0; i < 15; i++) {
            refB.append("TC");
            qryB.append(i < 12 ? "TT" : "TC");
        }
        for (int i = 0; i < 15; i++) {
            refB.append("GC");
            qryB.append("GC");
        }
        NucleotideSequence reference = new NucleotideSequence("ref", refB.toString());
        NucleotideSequence query = new NucleotideSequence("qry", qryB.toString());

        GcResult r = calc.compute(query, 50.0, reference);

        assertThat(r.mutationalSignature()).isNotNull();
        assertThat(r.mutationalSignature().anySignatureDetected()).isTrue();
        assertThat(r.category()).contains("APOBEC3-like C->T").contains("not necessarily HGT");
        assertThat(r.diagnostics()).anyMatch(d -> d.contains("SIGNIFICANT enrichment"));
    }

    @Test
    void fallsBackWithoutSignatureScanWhenReferenceLengthMismatches() {
        NucleotideSequence query = new NucleotideSequence("qry", "GCGCAAAAAA");
        NucleotideSequence reference = new NucleotideSequence("ref", "GCGCAAAA"); // shorter -- not aligned

        GcResult r = calc.compute(query, 38.0, reference);

        assertThat(r.mutationalSignature()).isNull();
        assertThat(r.diagnostics()).anyMatch(d -> d.contains("Mutational-signature scan skipped"));
    }
}
