package org.gvi.ui.chart;

import org.junit.jupiter.api.Test;

import static org.assertj.core.api.Assertions.assertThat;

class ContributionBarsTest {

    @Test
    void barWidthIsProportionalToContributionRelativeToTheMax() {
        assertThat(ContributionBars.barWidthPx(0.5, 1.0, 260)).isEqualTo(130.0);
        assertThat(ContributionBars.barWidthPx(1.0, 1.0, 260)).isEqualTo(260.0);
        assertThat(ContributionBars.barWidthPx(0.0, 1.0, 260)).isEqualTo(2.0);
    }

    @Test
    void negativeContributionsUseTheirMagnitude() {
        assertThat(ContributionBars.barWidthPx(-0.5, 1.0, 260))
                .isEqualTo(ContributionBars.barWidthPx(0.5, 1.0, 260));
    }

    @Test
    void aZeroOrNegativeMaxFallsBackToAMinimumVisibleWidthRatherThanDividingByZero() {
        assertThat(ContributionBars.barWidthPx(0.5, 0.0, 260)).isEqualTo(2.0);
    }
}
