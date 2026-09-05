package org.gvi.ui.pipeline;

import org.gvi.composite.IndexKey;

import java.util.List;

/**
 * One row in {@link PipelineStagesView}. {@code indexKeys} lists which composite-GVI indices this
 * stage's final status is derived from -- empty for structural stages (Sequence Input, GVI Scoring,
 * Final Report) that aren't tied to a single index.
 */
public record PipelineStage(String id, String label, List<IndexKey> indexKeys) {
}
