package org.gvi.core.util;

import org.junit.jupiter.api.Test;

import static org.assertj.core.api.Assertions.assertThat;

class FastaHeaderLocationHostParserTest {

    @Test
    void extractsLocationAndHostFromTheFourFieldConvention() {
        var result = FastaHeaderLocationHostParser.extract("JF416958.1|India|HUMAN|1957");
        assertThat(result).isPresent();
        assertThat(result.get().location()).isEqualTo("India");
        assertThat(result.get().host()).isEqualTo("HUMAN");
    }

    @Test
    void extractsFromTheUserReportedSheepDataHeader() {
        var result = FastaHeaderLocationHostParser.extract("KP821387.1|France|domestic_sheep|2001");
        assertThat(result).isPresent();
        assertThat(result.get().location()).isEqualTo("France");
        assertThat(result.get().host()).isEqualTo("domestic_sheep");
    }

    @Test
    void worksWithAFullIsoDateInTheLastField() {
        var result = FastaHeaderLocationHostParser.extract("ACC123|Kenya|bat|2019-04-02");
        assertThat(result).isPresent();
        assertThat(result.get().location()).isEqualTo("Kenya");
        assertThat(result.get().host()).isEqualTo("bat");
    }

    @Test
    void returnsEmptyForGisaidStyleTwoFieldHeader() {
        // "id|date|other" shapes must NOT be misread as "accession|location|host|date" --
        // conservative by design, no guessing at field meaning.
        assertThat(FastaHeaderLocationHostParser.extract("hCoV-19/USA/CA-1/2020|2020-03-15")).isEmpty();
    }

    @Test
    void returnsEmptyForThreeFieldNextstrainStyleHeader() {
        assertThat(FastaHeaderLocationHostParser.extract("sample|2021-06-01|USA")).isEmpty();
    }

    @Test
    void returnsEmptyWhenTheTrailingFieldIsNotADate() {
        assertThat(FastaHeaderLocationHostParser.extract("ACC1|India|HUMAN|not-a-date")).isEmpty();
    }

    @Test
    void returnsEmptyWhenAMiddleFieldIsBlank() {
        assertThat(FastaHeaderLocationHostParser.extract("ACC1||HUMAN|1957")).isEmpty();
        assertThat(FastaHeaderLocationHostParser.extract("ACC1|India||1957")).isEmpty();
    }

    @Test
    void returnsEmptyForNullOrPlainHeader() {
        assertThat(FastaHeaderLocationHostParser.extract(null)).isEmpty();
        assertThat(FastaHeaderLocationHostParser.extract("just_a_plain_sample_id")).isEmpty();
    }
}
