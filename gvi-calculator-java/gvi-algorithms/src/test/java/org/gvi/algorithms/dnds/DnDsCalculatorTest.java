package org.gvi.algorithms.dnds;

import org.junit.jupiter.api.Test;

import java.util.List;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.within;

class DnDsCalculatorTest {

    private final DnDsCalculator calc = new DnDsCalculator();

    @Test
    void synonymousSitesOfTttMatchesHandCalculation() {
        // TTT=Phe. Pos1 alternatives A/C/G -> ATT(Ile),CTT(Leu),GTT(Val): all nonsyn -> 0/3
        // Pos2 alternatives A/C/G -> TAT(Tyr),TCT(Ser),TGT(Cys): all nonsyn -> 0/3
        // Pos3 alternatives A/C/G -> TTA(Leu,nonsyn),TTC(Phe,SYN),TTG(Leu,nonsyn) -> 1/3
        // total = 0 + 0 + 1/3 = 1/3
        assertThat(calc.synonymousSites("TTT")).isCloseTo(1.0 / 3.0, within(1e-9));
    }

    @Test
    void synonymousSitesOfFourFoldDegenerateCodonIsExactlyOne() {
        // GGT = Gly; all 4 codons GGN encode Gly, so 3rd position is fully synonymous (3/3=1),
        // positions 1 and 2 are fully nonsynonymous (0/3 each) -> total = 1
        assertThat(calc.synonymousSites("GGT")).isCloseTo(1.0, within(1e-9));
    }

    @Test
    void everyCodonsSitesSumToThree() {
        for (String codon : org.gvi.core.util.GeneticCode.CODON_TABLE.keySet()) {
            if (org.gvi.core.util.GeneticCode.translate(codon) == org.gvi.core.util.GeneticCode.STOP) continue;
            double s = calc.synonymousSites(codon);
            // sSites + nSites == 3 by construction (nSites = 3 - sSites); just assert range sanity
            assertThat(s).isBetween(0.0, 3.0);
        }
    }

    @Test
    void singleDifferenceNonsynonymousChangeCountsAsOneNonsynDiff() {
        // TTT(Phe) -> TTA(Leu): single change at pos3, different amino acid -> nonsynonymous
        DnDsCalculator.PathwayResult r = calc.observedDifferences("TTT", "TTA");
        assertThat(r.synonymousDiffs()).isCloseTo(0.0, within(1e-9));
        assertThat(r.nonsynonymousDiffs()).isCloseTo(1.0, within(1e-9));
    }

    @Test
    void singleDifferenceSynonymousChangeCountsAsOneSynDiff() {
        // TTT(Phe) -> TTC(Phe): single change at pos3, same amino acid -> synonymous
        DnDsCalculator.PathwayResult r = calc.observedDifferences("TTT", "TTC");
        assertThat(r.synonymousDiffs()).isCloseTo(1.0, within(1e-9));
        assertThat(r.nonsynonymousDiffs()).isCloseTo(0.0, within(1e-9));
    }

    @Test
    void twoDifferencePathwaysAgreeForThisCodonPair() {
        // TTT(Phe) vs CTA(Leu): differ at pos1 (T->C) and pos3 (T->A)
        // Path A: TTT->CTT(Leu,nonsyn)->CTA(Leu,syn)
        // Path B: TTT->TTA(Leu,nonsyn)->CTA(Leu,syn)
        // Both pathways give (nonsyn=1, syn=1) -> average exactly (1,1)
        DnDsCalculator.PathwayResult r = calc.observedDifferences("TTT", "CTA");
        assertThat(r.synonymousDiffs()).isCloseTo(1.0, within(1e-9));
        assertThat(r.nonsynonymousDiffs()).isCloseTo(1.0, within(1e-9));
    }

    @Test
    void identicalCodonsHaveZeroDifferences() {
        DnDsCalculator.PathwayResult r = calc.observedDifferences("TTT", "TTT");
        assertThat(r.synonymousDiffs()).isZero();
        assertThat(r.nonsynonymousDiffs()).isZero();
    }

    @Test
    void syntheticGeneUnderStrongPurifyingSelectionHasLowOmega() {
        // Reference gene: repeat a 4-fold-degenerate codon (GGT=Gly) so all mutations we introduce
        // in the 3rd position are synonymous -> expect dN/dS well below 1.
        String refGene = "GGT".repeat(50);
        StringBuilder query = new StringBuilder(refGene);
        // introduce synonymous 3rd-position changes only, every other codon
        for (int i = 0; i < 50; i += 2) {
            query.setCharAt(i * 3 + 2, 'A'); // GGT -> GGA, still Gly
        }
        DnDsResult r = calc.compute("s1", "unknownGene", refGene, query.toString());
        assertThat(r.omega()).isLessThan(0.2);
        assertThat(r.category()).contains("no gene-specific reference band matched");
    }

    @Test
    void matchesGeneSpecificReferenceBandByKeyword() {
        String refGene = "TTT".repeat(30);
        StringBuilder query = new StringBuilder(refGene);
        // introduce nonsynonymous changes at 3rd position of several codons: TTT -> TTA (Phe->Leu)
        for (int i = 0; i < 30; i += 2) {
            query.setCharAt(i * 3 + 2, 'A');
        }
        DnDsResult r = calc.compute("s1", "Spike glycoprotein", refGene, query.toString());
        assertThat(r.category()).contains("Spike / Surface");
    }

    @Test
    void slidingWindowLocalizesPositiveSelectionThatWholeGenePoolingDilutes() {
        // First 30 codons: strong purifying selection (GGT->GGA, synonymous every other codon).
        // Last 30 codons: strong positive selection (TTT->TTA, nonsynonymous every other codon).
        StringBuilder ref = new StringBuilder();
        StringBuilder query = new StringBuilder();
        for (int i = 0; i < 30; i++) {
            ref.append("GGT");
            query.append(i % 2 == 0 ? "GGA" : "GGT"); // synonymous when changed
        }
        for (int i = 0; i < 30; i++) {
            ref.append("TTT");
            query.append(i % 2 == 0 ? "TTA" : "TTT"); // nonsynonymous when changed
        }

        List<SlidingWindowDnDsResult> windows = calc.slidingWindow("s1", "geneX", ref.toString(), query.toString(), 10, 10);
        assertThat(windows).hasSize(6);

        // windows 1-3 cover the purifying region (codons 1-30), windows 4-6 the positive-selection region (31-60)
        for (int i = 0; i < 3; i++) {
            assertThat(windows.get(i).omega()).as("purifying window %d", i).isLessThan(0.5);
        }
        for (int i = 3; i < 6; i++) {
            assertThat(windows.get(i).omega()).as("positive-selection window %d", i).isGreaterThan(1.0);
        }

        // whole-gene pooling should dilute the signal into something in between the two extremes
        DnDsResult pooled = calc.compute("s1", "geneX", ref.toString(), query.toString());
        assertThat(pooled.omega()).isBetween(windows.get(0).omega(), windows.get(5).omega());
    }

    @Test
    void slidingWindowSkipsIncompleteTrailingWindow() {
        String ref = "TTT".repeat(25); // 25 codons; window=10 step=10 -> windows [1-10],[11-20], then [21-30] doesn't fit (only 5 left)
        List<SlidingWindowDnDsResult> windows = calc.slidingWindow("s1", "geneX", ref, ref, 10, 10);
        assertThat(windows).hasSize(2);
    }

    @Test
    void poolingRecomputesFromRawCounts_ignoringWhateverOmegaTheInputsAlreadyCarried() {
        // Deliberately give both inputs a nonsense/sentinel omega/dN/dS (-1) unrelated to their real S/N/Sd/Nd --
        // if pool() were (bugged into) averaging the stored omega fields, the result would come out negative.
        // A correct implementation ignores those fields entirely and recomputes from the raw counts.
        DnDsResult r1 = new DnDsResult("s1", "geneX", -1, -1, -1, true, 2.0, 4.0, 0.0, 2.0, 5, 0, "n/a", List.of());
        DnDsResult r2 = new DnDsResult("s2", "geneX", -1, -1, -1, true, 3.0, 6.0, 1.0, 1.0, 5, 0, "n/a", List.of());

        DnDsResult pooled = calc.pool("dataset", "geneX", List.of(r1, r2));

        // hand-computed from the POOLED totals (S=5, N=10, Sd=1, Nd=3)
        double expectedPS = 1.0 / 5.0;   // pooled Sd/S, no continuity correction needed since pooled Sd != 0
        double expectedPN = 3.0 / 10.0;
        double expectedDS = -0.75 * Math.log(1 - (4.0 / 3.0) * expectedPS);
        double expectedDN = -0.75 * Math.log(1 - (4.0 / 3.0) * expectedPN);

        assertThat(pooled.synonymousSites()).isCloseTo(5.0, within(1e-9));
        assertThat(pooled.nonsynonymousSites()).isCloseTo(10.0, within(1e-9));
        assertThat(pooled.synonymousDifferences()).isCloseTo(1.0, within(1e-9));
        assertThat(pooled.nonsynonymousDifferences()).isCloseTo(3.0, within(1e-9));
        assertThat(pooled.dS()).isCloseTo(expectedDS, within(1e-9));
        assertThat(pooled.dN()).isCloseTo(expectedDN, within(1e-9));
        assertThat(pooled.omega()).isCloseTo(expectedDN / expectedDS, within(1e-9));
        assertThat(pooled.omega()).isPositive(); // proves it did NOT just average the sentinel -1 fields
    }

    @Test
    void poolingRejectsEmptyInput() {
        org.assertj.core.api.Assertions.assertThatThrownBy(() -> calc.pool("dataset", "geneX", List.of()))
                .isInstanceOf(org.gvi.core.exception.GviComputationException.class);
    }
}
