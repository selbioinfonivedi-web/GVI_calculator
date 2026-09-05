package org.gvi.ui;

import org.gvi.cli.PipelineResult;
import org.gvi.composite.IndexKey;
import org.gvi.composite.SensitivityResult;
import org.junit.jupiter.api.Test;

import java.util.List;
import java.util.Map;

import static org.assertj.core.api.Assertions.assertThat;

class SensitivityRowTest {

    @Test
    void sortsRowsBySpreadDescending() {
        List<SensitivityResult> sensitivity = List.of(
                new SensitivityResult(IndexKey.RE, 0.5, 0.48, 0.52),   // spread 0.04
                new SensitivityResult(IndexKey.MU, 0.5, 0.30, 0.70));  // spread 0.40

        PipelineResult result = new PipelineResult(Map.of(), Map.of(), null, Map.of(), Map.of(),
                Map.of(), Map.of(), Map.of(), Map.of(), sensitivity, List.of(), List.of(), null, List.of());

        List<SensitivityRow> rows = SensitivityRow.fromResult(result);
        assertThat(rows).hasSize(2);
        assertThat(rows.get(0).getIndex()).isEqualTo("mu");
        assertThat(rows.get(0).getSpread()).isEqualTo("0.4000");
        assertThat(rows.get(1).getIndex()).isEqualTo("Re");
    }

    @Test
    void emptyWhenNoSensitivityWasComputed() {
        PipelineResult result = new PipelineResult(Map.of(), Map.of(), null, Map.of(), Map.of(),
                Map.of(), Map.of(), Map.of(), Map.of(), List.of(), List.of(), List.of(), null, List.of());

        assertThat(SensitivityRow.fromResult(result)).isEmpty();
    }
}
