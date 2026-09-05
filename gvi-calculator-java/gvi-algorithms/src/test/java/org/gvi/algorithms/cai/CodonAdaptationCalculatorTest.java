package org.gvi.algorithms.cai;

import org.gvi.core.model.CodonUsageTable;
import org.junit.jupiter.api.Test;

import java.util.Map;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.within;

class CodonAdaptationCalculatorTest {

    private final CodonAdaptationCalculator calc = new CodonAdaptationCalculator();

    @Test
    void allOptimalCodonsGiveCaiOfOne() {
        // TTT (F) is the only codon used for F in this fake table (freq 1.0), TTC freq 0
        Map<String, Double> freq = Map.of("TTT", 1.0, "TTC", 0.0);
        CodonUsageTable table = CodonUsageTable.fromFrequencies("test", freq);
        // gene made entirely of the optimal codon TTT
        String cds = "TTT".repeat(10);
        CaiResult r = calc.compute("s1", "geneX", cds, table);
        assertThat(r.cai()).isCloseTo(1.0, within(1e-9));
        assertThat(r.codonsUsed()).isEqualTo(10);
    }

    @Test
    void matchesHandCalculatedGeometricMeanForMixedCodons() {
        // Phe: TTT freq=6 (max), TTC freq=2 -> w(TTT)=1.0, w(TTC)=2/6=0.3333...
        Map<String, Double> freq = Map.of("TTT", 6.0, "TTC", 2.0);
        CodonUsageTable table = CodonUsageTable.fromFrequencies("test", freq);
        String cds = "TTT" + "TTC"; // 2 codons: w=1.0, w=0.3333
        CaiResult r = calc.compute("s1", "geneX", cds, table);
        double expected = Math.exp((Math.log(1.0) + Math.log(2.0 / 6.0)) / 2.0);
        assertThat(r.cai()).isCloseTo(expected, within(1e-9));
    }

    @Test
    void excludesStopCodonsAndSingleCodonAminoAcidsFromL() {
        Map<String, Double> freq = Map.of("TTT", 1.0, "TTC", 1.0, "ATG", 1.0, "TGG", 1.0, "TAA", 1.0);
        CodonUsageTable table = CodonUsageTable.fromFrequencies("test", freq);
        // ATG (Met, single-codon) + TGG (Trp, single-codon) + TAA (stop) + TTT (Phe, counted)
        String cds = "ATG" + "TGG" + "TTT" + "TAA";
        CaiResult r = calc.compute("s1", "geneX", cds, table);
        assertThat(r.codonsUsed()).isEqualTo(1); // only TTT counted
        assertThat(r.excludedStop()).isEqualTo(1);
        assertThat(r.excludedSingleAminoAcid()).isEqualTo(2);
    }

    @Test
    void dropsTrailingPartialCodon() {
        Map<String, Double> freq = Map.of("TTT", 1.0, "TTC", 1.0);
        CodonUsageTable table = CodonUsageTable.fromFrequencies("test", freq);
        String cds = "TTT" + "AC"; // trailing 2 bases dropped
        CaiResult r = calc.compute("s1", "geneX", cds, table);
        assertThat(r.codonsUsed()).isEqualTo(1);
        assertThat(r.diagnostics()).anyMatch(s -> s.contains("trailing"));
    }
}
