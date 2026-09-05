package org.gvi.ui;

import org.gvi.cli.PipelineResult;
import org.gvi.composite.GviResult;
import org.gvi.composite.IndexKey;
import org.gvi.core.spi.IndexResult;
import org.gvi.core.spi.SimpleIndexResult;
import org.junit.jupiter.api.Test;

import java.util.EnumMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import static org.assertj.core.api.Assertions.assertThat;

class SequenceRowTest {

    @Test
    void mapsPopulationAndPerSequenceIndicesIntoOneRowPerSequence() {
        Map<IndexKey, IndexResult> population = new EnumMap<>(IndexKey.class);
        population.put(IndexKey.PI, new SimpleIndexResult("pi", 0.0123, "test", List.of()));
        population.put(IndexKey.MU, new SimpleIndexResult("mu", 0.00045, "test", List.of()));

        Map<String, Map<IndexKey, IndexResult>> perSequence = new LinkedHashMap<>();
        Map<IndexKey, IndexResult> seq1 = new EnumMap<>(IndexKey.class);
        seq1.put(IndexKey.GD, new SimpleIndexResult("GD", 0.0067, "test", List.of()));
        seq1.put(IndexKey.MB, new SimpleIndexResult("MB", 12.0, "test", List.of()));
        perSequence.put("seq1", seq1);

        Map<String, GviResult> gviPerSequence = new LinkedHashMap<>();
        gviPerSequence.put("seq1", new GviResult(0.4321, List.of(), List.of(), 1.0, List.of()));

        PipelineResult result = new PipelineResult(population, Map.of(), null, perSequence, gviPerSequence,
                Map.of(), Map.of(), Map.of(), Map.of(), List.of(), List.of(), List.of(), null, List.of());

        List<SequenceRow> rows = SequenceRow.fromResult(result);
        assertThat(rows).hasSize(1);
        SequenceRow row = rows.get(0);

        assertThat(row.getSequenceId()).isEqualTo("seq1");
        assertThat(row.getGvi()).isEqualTo("0.4321");
        assertThat(row.getPi()).isEqualTo(String.format("%.4g", 0.0123));
        assertThat(row.getMu()).isEqualTo(String.format("%.4g", 0.00045));
        assertThat(row.getGd()).isEqualTo(String.format("%.4g", 0.0067));
        assertThat(row.getMb()).isEqualTo(String.format("%.4g", 12.0));
        // indices never computed for this sequence/dataset must show as a placeholder, not crash or show blank/null
        assertThat(row.getRe()).isEqualTo("-");
        assertThat(row.getDnds()).isEqualTo("-");
        assertThat(row.getCai()).isEqualTo("-");
        assertThat(row.getGc()).isEqualTo("-");
        assertThat(row.getRi()).isEqualTo("-");
    }

    @Test
    void showsPlaceholderWhenGviWasNotComputedForASequence() {
        Map<String, Map<IndexKey, IndexResult>> perSequence = new LinkedHashMap<>();
        perSequence.put("seq1", new EnumMap<>(IndexKey.class));

        PipelineResult result = new PipelineResult(Map.of(), Map.of(), null, perSequence, Map.of(),
                Map.of(), Map.of(), Map.of(), Map.of(), List.of(), List.of(), List.of(), null, List.of());

        List<SequenceRow> rows = SequenceRow.fromResult(result);
        assertThat(rows.get(0).getGvi()).isEqualTo("-");
    }
}
