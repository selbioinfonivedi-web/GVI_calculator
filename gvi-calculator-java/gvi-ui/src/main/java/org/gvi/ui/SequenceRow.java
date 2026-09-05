package org.gvi.ui;

import org.gvi.cli.PipelineResult;
import org.gvi.composite.GviResult;
import org.gvi.composite.IndexKey;
import org.gvi.core.spi.IndexResult;

import java.util.EnumMap;
import java.util.Map;

/**
 * One row of the results table: a query sequence's GVI and its component
 * index values (population-wide indices like pi/mu/Re are repeated on
 * every row for easy comparison against per-sequence ones like GD/MB).
 * Plain JavaBean-style getters so JavaFX's {@code PropertyValueFactory}
 * can bind table columns to them by name without extra glue code.
 */
public final class SequenceRow {

    private final String sequenceId;
    private final Double gvi;
    private final Map<IndexKey, IndexResult> values;

    SequenceRow(String sequenceId, Double gvi, Map<IndexKey, IndexResult> values) {
        this.sequenceId = sequenceId;
        this.gvi = gvi;
        this.values = values;
    }

    public static java.util.List<SequenceRow> fromResult(PipelineResult result) {
        java.util.List<SequenceRow> rows = new java.util.ArrayList<>();
        for (var entry : result.perSequenceIndices().entrySet()) {
            String id = entry.getKey();
            Map<IndexKey, IndexResult> combined = new EnumMap<>(entry.getValue());
            combined.putAll(result.populationIndices());
            GviResult gviResult = result.gviPerSequence().get(id);
            rows.add(new SequenceRow(id, gviResult == null ? null : gviResult.gvi(), combined));
        }
        return rows;
    }

    public String getSequenceId() {
        return sequenceId;
    }

    public String getGvi() {
        return gvi == null ? "-" : String.format("%.4f", gvi);
    }

    public String getMu() {
        return format(IndexKey.MU);
    }

    public String getRe() {
        return format(IndexKey.RE);
    }

    public String getPi() {
        return format(IndexKey.PI);
    }

    public String getMb() {
        return format(IndexKey.MB);
    }

    public String getDnds() {
        return format(IndexKey.DNDS);
    }

    public String getGd() {
        return format(IndexKey.GD);
    }

    public String getCai() {
        return format(IndexKey.CAI);
    }

    public String getGc() {
        return format(IndexKey.GC);
    }

    public String getRi() {
        return format(IndexKey.RI);
    }

    private String format(IndexKey key) {
        IndexResult r = values.get(key);
        return r == null ? "-" : String.format("%.4g", r.primaryValue());
    }
}
