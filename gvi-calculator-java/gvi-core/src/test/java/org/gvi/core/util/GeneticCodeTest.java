package org.gvi.core.util;

import org.gvi.core.exception.GviComputationException;
import org.junit.jupiter.api.Test;

import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.assertThatThrownBy;

/**
 * {@link GeneticCode} is a hand-typed 64-entry table that every codon-level result rests on: dN/dS
 * counts synonymous versus non-synonymous sites through it, and CAI weights codons by it. A single
 * transposed letter would shift both silently -- no exception, no warning, just a wrong number -- and
 * the class was at 0% coverage.
 * <p>
 * The central test therefore does <em>not</em> retype the table, which would only reproduce whatever
 * typo it was checking for. It checks against NCBI translation table 1 as NCBI publishes it: a single
 * string in a different ordering (base1 outer, base2 middle, base3 inner, each cycling TCAG) from the
 * source's own layout. A transcription slip in one will not match the other.
 */
class GeneticCodeTest {

    /**
     * NCBI translation table 1, verbatim from the published {@code AAs} line. The codon for position
     * i is Base1[i] Base2[i] Base3[i] with the bases cycling as documented there.
     */
    private static final String NCBI_AAS =
            "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";

    private static final char[] TCAG = {'T', 'C', 'A', 'G'};

    /** The 64 codons in NCBI's own order, so the comparison is against a differently-ordered source. */
    private static String ncbiCodonAt(int i) {
        return "" + TCAG[i / 16] + TCAG[(i / 4) % 4] + TCAG[i % 4];
    }

    @Test
    void everyCodonAgreesWithNcbiTranslationTableOne() {
        assertThat(NCBI_AAS).hasSize(64);

        for (int i = 0; i < 64; i++) {
            String codon = ncbiCodonAt(i);
            assertThat(GeneticCode.translate(codon))
                    .as("codon %s (NCBI position %d)", codon, i)
                    .isEqualTo(NCBI_AAS.charAt(i));
        }
    }

    @Test
    void theTableHasAllSixtyFourCodonsAndNothingElse() {
        assertThat(GeneticCode.CODON_TABLE).hasSize(64);

        Set<String> expected = new HashSet<>();
        for (int i = 0; i < 64; i++) expected.add(ncbiCodonAt(i));
        assertThat(GeneticCode.CODON_TABLE.keySet()).isEqualTo(expected);
    }

    /**
     * Three stops and 61 sense codons is the defining property of the standard code. dN/dS site
     * counting treats stops specially, so an extra or missing stop changes the denominator.
     */
    @Test
    void thereAreExactlyThreeStopCodons() {
        List<String> stops = GeneticCode.CODON_TABLE.entrySet().stream()
                .filter(e -> e.getValue() == GeneticCode.STOP)
                .map(Map.Entry::getKey)
                .sorted()
                .toList();

        assertThat(stops).containsExactly("TAA", "TAG", "TGA");
        assertThat(GeneticCode.CODON_TABLE).hasSize(61 + stops.size());
    }

    /** Met and Trp are the only single-codon amino acids; this pins the degeneracy structure. */
    @Test
    void methionineAndTryptophanAreTheOnlyUniquelyEncodedAminoAcids() {
        List<Character> singletons = GeneticCode.SYNONYMOUS_CODONS.entrySet().stream()
                .filter(e -> e.getValue().size() == 1)
                .map(Map.Entry::getKey)
                .sorted()
                .toList();

        assertThat(singletons).containsExactly('M', 'W');
        assertThat(GeneticCode.SYNONYMOUS_CODONS.get('M')).containsExactly("ATG");
        assertThat(GeneticCode.SYNONYMOUS_CODONS.get('W')).containsExactly("TGG");
    }

    @Test
    void theCodeEncodesTwentyAminoAcidsPlusStop() {
        assertThat(GeneticCode.SYNONYMOUS_CODONS).hasSize(21);
        assertThat(GeneticCode.SYNONYMOUS_CODONS).doesNotContainKey('X');
    }

    /**
     * CAI divides a codon's usage by the maximum within its synonymous family, so the families must
     * partition the 64 codons exactly: a codon in two families, or in none, corrupts that ratio.
     */
    @Test
    void theSynonymousFamiliesPartitionTheTableExactly() {
        List<String> all = GeneticCode.SYNONYMOUS_CODONS.values().stream().flatMap(List::stream).toList();

        assertThat(all).as("no codon may appear in two families").doesNotHaveDuplicates();
        assertThat(all).hasSize(64);
        assertThat(new HashSet<>(all)).isEqualTo(GeneticCode.CODON_TABLE.keySet());

        // and every family really encodes the amino acid it is filed under
        GeneticCode.SYNONYMOUS_CODONS.forEach((aa, codons) ->
                assertThat(codons).allSatisfy(c ->
                        assertThat(GeneticCode.CODON_TABLE.get(c)).as("family %s contains %s", aa, c).isEqualTo(aa)));
    }

    /** The six-fold families are the ones a naive table most often gets wrong: they span two blocks. */
    @Test
    void theSixFoldDegenerateFamiliesAreComplete() {
        assertThat(GeneticCode.SYNONYMOUS_CODONS.get('L'))
                .containsExactlyInAnyOrder("TTA", "TTG", "CTT", "CTC", "CTA", "CTG");
        assertThat(GeneticCode.SYNONYMOUS_CODONS.get('S'))
                .containsExactlyInAnyOrder("TCT", "TCC", "TCA", "TCG", "AGT", "AGC");
        assertThat(GeneticCode.SYNONYMOUS_CODONS.get('R'))
                .containsExactlyInAnyOrder("CGT", "CGC", "CGA", "CGG", "AGA", "AGG");
    }

    // ── input handling ───────────────────────────────────────────────────────
    @Test
    void rnaAndLowercaseInputTranslateAsDna() {
        assertThat(GeneticCode.translate("aug")).isEqualTo('M');
        assertThat(GeneticCode.translate("AUG")).isEqualTo('M');
        assertThat(GeneticCode.translate("uaa")).isEqualTo(GeneticCode.STOP);
        assertThat(GeneticCode.translate("ATG")).isEqualTo('M');
    }

    /**
     * Gaps and ambiguity codes are ordinary in a real alignment. Translating one must fail loudly
     * rather than return a plausible amino acid -- the callers skip such codons, and they can only do
     * that if the code refuses to guess.
     */
    @Test
    void gapsAndAmbiguityCodesAreRefusedRatherThanGuessed() {
        for (String bad : new String[]{"AT-", "---", "ATN", "NNN", "ATR", "AT"}) {
            assertThatThrownBy(() -> GeneticCode.translate(bad))
                    .as("translating %s", bad)
                    .isInstanceOf(GviComputationException.class)
                    .hasMessageContaining(bad);
        }
    }

    @Test
    void isStandardCodonAgreesWithWhatTranslateAccepts() {
        assertThat(GeneticCode.isStandardCodon("ATG")).isTrue();
        assertThat(GeneticCode.isStandardCodon("aug")).isTrue();
        assertThat(GeneticCode.isStandardCodon("AT-")).isFalse();
        assertThat(GeneticCode.isStandardCodon("NNN")).isFalse();

        for (int i = 0; i < 64; i++) {
            assertThat(GeneticCode.isStandardCodon(ncbiCodonAt(i))).isTrue();
        }
    }

    @Test
    void isStopIdentifiesTheThreeStopsAndNothingElse() {
        assertThat(GeneticCode.isStop("TAA")).isTrue();
        assertThat(GeneticCode.isStop("TAG")).isTrue();
        assertThat(GeneticCode.isStop("TGA")).isTrue();
        assertThat(GeneticCode.isStop("UAA")).as("RNA spelling").isTrue();
        assertThat(GeneticCode.isStop("ATG")).isFalse();
        assertThat(GeneticCode.isStop("TGG")).as("TGG is Trp, not a stop -- the classic off-by-one").isFalse();
    }

    /** The tables are public and static; a caller mutating them would corrupt every later run. */
    @Test
    void thePublishedTablesAreImmutable() {
        assertThatThrownBy(() -> GeneticCode.CODON_TABLE.put("ATG", 'X'))
                .isInstanceOf(UnsupportedOperationException.class);
        assertThatThrownBy(() -> GeneticCode.SYNONYMOUS_CODONS.get('M').add("XXX"))
                .isInstanceOf(UnsupportedOperationException.class);
    }
}
