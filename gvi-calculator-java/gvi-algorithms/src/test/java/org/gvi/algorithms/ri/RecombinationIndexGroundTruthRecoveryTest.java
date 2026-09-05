package org.gvi.algorithms.ri;

import org.gvi.core.model.NucleotideSequence;
import org.gvi.core.model.SequenceAlignment;
import org.junit.jupiter.api.Test;

import java.util.List;
import java.util.Random;

import static org.assertj.core.api.Assertions.assertThat;

/**
 * Validates RI/the PHI test against alignments with a KNOWN recombination
 * status -- the same style of validation Bruen, Bryant &amp; Poss (2006)
 * used for the PHI test itself: simulated recombinant vs. clonal
 * alignments, not real data whose true recombination history is unknown.
 * <p>
 * Construction: two "pure" lineages (A and B) each carry a private
 * diagnostic allele at ~150 biallelic sites spread across a 900bp
 * alignment (two representatives of each, so every diagnostic site is
 * genuinely polymorphic in the sample -- required for the PHI test's
 * four-gamete counting). Two more taxa are added that are either:
 * <ul>
 *   <li>a complementary pair of mosaics ("R1"/"R2"): R1 carries clade-A's
 *   allele in the first half and clade-B's in the second half; R2 carries
 *   the opposite. A single mosaic alone only ever produces ONE of the two
 *   "cross" allele combinations at a pair of sites spanning the breakpoint
 *   (e.g. A-then-B), and the four-gamete test needs BOTH cross
 *   combinations present to declare incompatibility -- exactly the real
 *   biological point that a lone novel haplotype can still be explained by
 *   a compatible tree, but two complementary recombination products
 *   cannot; or</li>
 *   <li>two more non-mosaic clade-A representatives ("C1"/"C2") -- every
 *   diagnostic site supports the identical bipartition throughout, so
 *   there is nothing for the PHI test to find.</li>
 * </ul>
 */
class RecombinationIndexGroundTruthRecoveryTest {

    private static final long SEED = 20260813L;
    private static final char[] BASES = {'A', 'C', 'G', 'T'};
    private static final int LENGTH = 900;
    private static final int DIAGNOSTIC_STEP = 6; // ~150 diagnostic sites across 900bp

    @Test
    void flagsComplementaryMosaicRecombinantsAsSignificant() {
        Random rng = new Random(SEED);
        char[] ancestral = randomNonDiagnosticBackground(rng);

        char[] pureA1 = ancestral.clone();
        char[] pureA2 = ancestral.clone();
        char[] pureB1 = ancestral.clone();
        char[] pureB2 = ancestral.clone();
        char[] recombinant1 = ancestral.clone(); // A then B
        char[] recombinant2 = ancestral.clone(); // B then A

        for (int pos = 0; pos < LENGTH; pos += DIAGNOSTIC_STEP) {
            char alleleA = ancestral[pos];
            char alleleB = alternateBase(alleleA);
            pureA1[pos] = alleleA;
            pureA2[pos] = alleleA;
            pureB1[pos] = alleleB;
            pureB2[pos] = alleleB;
            boolean firstHalf = pos < LENGTH / 2;
            recombinant1[pos] = firstHalf ? alleleA : alleleB;
            recombinant2[pos] = firstHalf ? alleleB : alleleA;
        }

        SequenceAlignment alignment = SequenceAlignment.of(List.of(
                seq("pureA1", pureA1), seq("pureA2", pureA2),
                seq("pureB1", pureB1), seq("pureB2", pureB2),
                seq("recombinant1", recombinant1), seq("recombinant2", recombinant2)));

        RiResult result = new RecombinationIndexCalculator().computeFromAlignment(alignment);
        System.out.println("[RecombinationIndexGroundTruthRecoveryTest] recombinant: ri=" + result.ri()
                + " pValue=" + result.pValue() + " significant=" + result.significant());

        assertThat(result.significant()).isTrue();
        assertThat(result.pValue()).isLessThan(0.05);
        assertThat(result.ri()).isGreaterThan(0.3);
    }

    @Test
    void doesNotFlagAClonalNonMosaicControlAsSignificant() {
        Random rng = new Random(SEED);
        char[] ancestral = randomNonDiagnosticBackground(rng);

        char[] pureA1 = ancestral.clone();
        char[] pureA2 = ancestral.clone();
        char[] pureB1 = ancestral.clone();
        char[] pureB2 = ancestral.clone();
        char[] clonalControl1 = ancestral.clone(); // clade-A allele throughout -- no mosaic structure anywhere
        char[] clonalControl2 = ancestral.clone();

        for (int pos = 0; pos < LENGTH; pos += DIAGNOSTIC_STEP) {
            char alleleA = ancestral[pos];
            char alleleB = alternateBase(alleleA);
            pureA1[pos] = alleleA;
            pureA2[pos] = alleleA;
            pureB1[pos] = alleleB;
            pureB2[pos] = alleleB;
            clonalControl1[pos] = alleleA;
            clonalControl2[pos] = alleleA;
        }

        SequenceAlignment alignment = SequenceAlignment.of(List.of(
                seq("pureA1", pureA1), seq("pureA2", pureA2),
                seq("pureB1", pureB1), seq("pureB2", pureB2),
                seq("clonalControl1", clonalControl1), seq("clonalControl2", clonalControl2)));

        RiResult result = new RecombinationIndexCalculator().computeFromAlignment(alignment);
        System.out.println("[RecombinationIndexGroundTruthRecoveryTest] clonal control: ri=" + result.ri()
                + " pValue=" + result.pValue() + " significant=" + result.significant());

        assertThat(result.significant()).isFalse();
        assertThat(result.ri()).isLessThan(0.05);
    }

    private static NucleotideSequence seq(String id, char[] chars) {
        return new NucleotideSequence(id, new String(chars));
    }

    private static char[] randomNonDiagnosticBackground(Random rng) {
        char[] seq = new char[LENGTH];
        for (int i = 0; i < LENGTH; i++) seq[i] = BASES[rng.nextInt(4)];
        return seq;
    }

    /** A fixed, different base -- just cycles to the next letter in ACGT, no randomness needed since this only needs to differ from {@code base}. */
    private static char alternateBase(char base) {
        return BASES[("ACGT".indexOf(base) + 1) % 4];
    }
}
