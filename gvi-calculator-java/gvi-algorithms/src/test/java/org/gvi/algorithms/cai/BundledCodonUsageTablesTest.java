package org.gvi.algorithms.cai;

import org.gvi.core.exception.GviInputException;
import org.gvi.core.model.CodonUsageTable;
import org.junit.jupiter.api.Test;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.assertThatThrownBy;

/**
 * Confirms every bundled species resource is present in the jar, parses
 * cleanly, and produces sane weights -- these are real data computed from
 * NCBI RefSeq CDS FASTA files (see {@link BundledCodonUsageTables}'s class
 * javadoc for full provenance: exact assemblies, CDS/codon counts), not
 * synthetic fixtures, so this test is really just checking the
 * bundling/loading mechanism, not the data itself.
 */
class BundledCodonUsageTablesTest {

    @Test
    void loadsEveryAdvertisedSpeciesWithoutError() {
        for (String species : BundledCodonUsageTables.availableSpecies()) {
            CodonUsageTable table = BundledCodonUsageTables.load(species);
            // ATG (Met) is always the sole codon for its amino acid, so its weight is always exactly 1.0 --
            // a cheap sanity check that the table loaded and computed weights at all.
            assertThat(table.weightOf("ATG")).isEqualTo(1.0);
        }
    }

    @Test
    void isCaseInsensitiveAndTrimsWhitespace() {
        assertThat(BundledCodonUsageTables.load(" Human ")).isNotNull();
        assertThat(BundledCodonUsageTables.load("HUMAN")).isNotNull();
    }

    @Test
    void rejectsAnUnknownSpeciesWithAClearMessage() {
        assertThatThrownBy(() -> BundledCodonUsageTables.load("dodo"))
                .isInstanceOf(GviInputException.class)
                .hasMessageContaining("dodo")
                .hasMessageContaining("human");
    }

    @Test
    void advertisesAllBundledSpecies() {
        assertThat(BundledCodonUsageTables.availableSpecies())
                .contains("human", "mouse", "pig", "wild_boar", "cattle", "buffalo", "sheep", "goat", "horse",
                        "ecoli", "aedes_aegypti");
    }

    @Test
    void wildBoarIsAnAliasForPigSinceTheyAreTheSameSpecies() {
        assertThat(BundledCodonUsageTables.load("wild_boar").weightOf("CTG"))
                .isEqualTo(BundledCodonUsageTables.load("pig").weightOf("CTG"));
    }
}
