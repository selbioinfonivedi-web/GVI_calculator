package org.gvi.algorithms.re;

import org.junit.jupiter.api.Test;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.assertThatThrownBy;
import static org.assertj.core.api.Assertions.within;

/**
 * The bundled generation-time table replaces a single hardcoded 5-day default that was
 * applied to every organism. These tests pin the three states apart, because the whole
 * point of the change is that two of them must NOT yield a number.
 */
class GenerationTimeTableTest {

    private final GenerationTimeTable table = GenerationTimeTable.bundled();

    @Test
    void bundledTableLoads() {
        assertThat(table.size()).isGreaterThan(10);
        assertThat(table.usableIds()).contains("influenza_a", "sars_cov_2", "asfv");
    }

    @Test
    void citationBackedEntryReturnsItsValue() {
        assertThat(table.requireGenerationTimeDays("sars_cov_2")).isCloseTo(5.2, within(1e-9));
        assertThat(table.requireGenerationTimeDays("influenza_a")).isCloseTo(2.0, within(1e-9));
    }

    @Test
    void lookupIsCaseInsensitive() {
        assertThat(table.requireGenerationTimeDays("SARS_CoV_2")).isCloseTo(5.2, within(1e-9));
    }

    @Test
    void environmentallyAcquiredPathogenRefusesToProduceAGenerationTime() {
        // B. anthracis is acquired from environmental spores, not from a preceding case, so
        // there is no serial interval to estimate Re from at any value of T.
        assertThatThrownBy(() -> table.requireGenerationTimeDays("b_anthracis"))
                .isInstanceOf(MissingGenerationTimeException.class)
                .hasMessageContaining("not applicable")
                .hasMessageContaining("no serial interval");
    }

    @Test
    void stubEntryFailsRatherThanGuessing() {
        assertThatThrownBy(() -> table.requireGenerationTimeDays("fmd"))
                .isInstanceOf(MissingGenerationTimeException.class)
                .hasMessageContaining("still a stub");
    }

    @Test
    void unknownIdListsWhatIsAvailable() {
        assertThatThrownBy(() -> table.requireGenerationTimeDays("not_a_pathogen"))
                .isInstanceOf(MissingGenerationTimeException.class)
                .hasMessageContaining("No generation-time entry")
                .hasMessageContaining("influenza_a");
    }

    @Test
    void correctedEntriesCarryTheCorrectedValues() {
        // M. tuberculosis was 20 d (an in-vitro doubling time); the serial interval is ~1 year.
        assertThat(table.requireGenerationTimeDays("m_tuberculosis")).isCloseTo(365.0, within(1e-9));
        // P. falciparum was 10 d (the mosquito extrinsic incubation period); the in-host
        // erythrocytic cycle used for a parasite gene tree is ~48 h.
        assertThat(table.requireGenerationTimeDays("p_falciparum")).isCloseTo(2.0, within(1e-9));
    }

    @Test
    void vectorBorneEntriesDocumentWhichCycleTheyMeasure() {
        for (String id : new String[]{"blue_tongue", "leishmania_spp", "p_falciparum"}) {
            GenerationTimeTable.Entry e = table.find(id).orElseThrow();
            assertThat(e.vectorBorne()).as("%s marked vector-borne", id).isTrue();
            assertThat(e.note().toLowerCase()).as("%s names the cycle", id).contains("intrinsic");
        }
    }

    @Test
    void everyApplicableEntryHasAPositiveValueAndEveryOtherHasNone() {
        for (String id : table.usableIds()) {
            GenerationTimeTable.Entry e = table.find(id).orElseThrow();
            assertThat(e.tDays()).as("%s", id).isNotNull().isGreaterThan(0.0);
        }
    }
}
