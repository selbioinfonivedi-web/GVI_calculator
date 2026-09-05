package org.gvi.algorithms.phylo;

import org.gvi.algorithms.gd.GdMethod;
import org.gvi.core.model.NucleotideSequence;
import org.gvi.core.model.SequenceAlignment;
import org.junit.jupiter.api.Test;

import java.util.List;
import java.util.Set;

import static org.assertj.core.api.Assertions.assertThat;

class BootstrapSupportCalculatorTest {

    /**
     * 5 taxa built so the true topology is ((a1,a2),(b1,b2)) rooted at ref:
     * a1/a2 share 10 identical columns (block2) found nowhere else, b1/b2
     * share 10 different identical columns (block3) found nowhere else, and
     * a 10bp block of small private per-taxon variation adds terminal branch
     * length without threatening either clade. With 20 of 40 columns
     * strongly and redundantly supporting the two true clades, a resampled
     * alignment would need to miss ALL 10 copies of a block simultaneously
     * to lose that clade's signal (~0.75^40, vanishingly unlikely) -- so a
     * real bootstrap implementation must recover both clades in nearly
     * every replicate.
     */
    @Test
    void recoversHighSupportForCladesWithStrongRedundantSignal() {
        String block1 = "A".repeat(10); // invariant
        NucleotideSequence ref = new NucleotideSequence("ref", block1 + "A".repeat(10) + "A".repeat(10) + "A".repeat(10));
        NucleotideSequence a1 = new NucleotideSequence("a1", block1 + "C".repeat(10) + "A".repeat(10) + "AAAAAAAAGA");
        NucleotideSequence a2 = new NucleotideSequence("a2", block1 + "C".repeat(10) + "A".repeat(10) + "AAAAATAAAA");
        NucleotideSequence b1 = new NucleotideSequence("b1", block1 + "A".repeat(10) + "G".repeat(10) + "AAACAAAAAA");
        NucleotideSequence b2 = new NucleotideSequence("b2", block1 + "A".repeat(10) + "G".repeat(10) + "AAAAAAATAA");

        SequenceAlignment alignment = SequenceAlignment.of(List.of(ref, a1, a2, b1, b2), "ref");
        PhyloTreeFactory.BuiltTree built = PhyloTreeFactory.build(alignment, GdMethod.JUKES_CANTOR);

        BootstrapSupportCalculator.Support support = new BootstrapSupportCalculator(200, 42L)
                .compute(alignment, GdMethod.JUKES_CANTOR, built.tree());

        assertThat(support.successfulReplicates()).isGreaterThan(0);
        assertThat(support.supportByClade()).hasSize(2); // exactly the 2 informative internal splits for 5 taxa

        Double cladeA = support.forClade(Set.of("a1", "a2"));
        Double cladeB = support.forClade(Set.of("b1", "b2"));
        assertThat(cladeA).isNotNull();
        assertThat(cladeB).isNotNull();
        assertThat(cladeA).isGreaterThanOrEqualTo(90.0);
        assertThat(cladeB).isGreaterThanOrEqualTo(90.0);
    }

    /** For 3 taxa, the only branch point is root's own child (the trivial, always-present split) -- nothing informative to test. */
    @Test
    void reportsNoInformativeSplitsForThreeTaxa() {
        NucleotideSequence ref = new NucleotideSequence("ref", "ACGTACGTAC");
        NucleotideSequence q1 = new NucleotideSequence("q1", "ACGTACGTAG");
        NucleotideSequence q2 = new NucleotideSequence("q2", "ACGTACGTCC");
        SequenceAlignment alignment = SequenceAlignment.of(List.of(ref, q1, q2), "ref");
        PhyloTreeFactory.BuiltTree built = PhyloTreeFactory.build(alignment, GdMethod.JUKES_CANTOR);

        BootstrapSupportCalculator.Support support = new BootstrapSupportCalculator(50, 42L)
                .compute(alignment, GdMethod.JUKES_CANTOR, built.tree());

        assertThat(support.supportByClade()).isEmpty();
        assertThat(support.diagnostics()).anyMatch(d -> d.contains("No informative internal splits"));
    }
}
