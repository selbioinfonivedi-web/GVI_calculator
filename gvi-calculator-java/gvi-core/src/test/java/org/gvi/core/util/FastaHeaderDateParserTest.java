package org.gvi.core.util;

import org.junit.jupiter.api.Test;

import java.time.LocalDate;

import static org.assertj.core.api.Assertions.assertThat;

class FastaHeaderDateParserTest {

    @Test
    void extractsIsoDateFromGisaidStyleHeader() {
        var result = FastaHeaderDateParser.extract("hCoV-19/USA/CA-1/2020|2020-03-15|EPI_ISL_12345");
        assertThat(result).isPresent();
        assertThat(result.get().date()).isEqualTo(LocalDate.of(2020, 3, 15));
        assertThat(result.get().yearOnly()).isFalse();
    }

    @Test
    void extractsIsoDateFromUnderscoreSeparatedHeader() {
        var result = FastaHeaderDateParser.extract("sample1_2021-06-01_India");
        assertThat(result).isPresent();
        assertThat(result.get().date()).isEqualTo(LocalDate.of(2021, 6, 1));
    }

    @Test
    void fallsBackToYearOnlyWhenNoFullDatePresent() {
        var result = FastaHeaderDateParser.extract("sample_collected_2021_somewhere");
        assertThat(result).isPresent();
        // Midpoint of the year, not January 1st -- an unbiased single-date estimate when only the year is
        // known (see TemporalUtil#midpointOfYear); resolving to Jan 1st would systematically bias every
        // downstream root-to-tip regression toward the start of the year.
        assertThat(result.get().date()).isEqualTo(org.gvi.core.util.TemporalUtil.midpointOfYear(2021));
        assertThat(result.get().yearOnly()).isTrue();
    }

    @Test
    void returnsEmptyWhenNoDateInHeader() {
        assertThat(FastaHeaderDateParser.extract("just_a_plain_sample_id")).isEmpty();
    }

    @Test
    void skipsInvalidCalendarDateAndKeepsScanning() {
        // "2024-99-99" isn't a real date (month 99 invalid) -- should be skipped, not crash,
        // and the parser should keep looking rather than giving up entirely.
        var result = FastaHeaderDateParser.extract("build-2024-99-99_realdate-2021-05-20");
        assertThat(result).isPresent();
        assertThat(result.get().date()).isEqualTo(LocalDate.of(2021, 5, 20));
    }

    @Test
    void returnsEmptyForNullHeader() {
        assertThat(FastaHeaderDateParser.extract(null)).isEmpty();
    }
}
