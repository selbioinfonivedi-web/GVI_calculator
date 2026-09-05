package org.gvi.cli;

import org.gvi.algorithms.cai.CaiResult;
import org.gvi.algorithms.dnds.DnDsResult;
import org.gvi.algorithms.dnds.SlidingWindowDnDsResult;
import org.gvi.algorithms.dnds.ml.MlDnDsResult;
import org.gvi.algorithms.ri.ReassortmentResult;
import org.gvi.composite.GviResult;
import org.gvi.composite.IndexKey;
import org.gvi.composite.SensitivityResult;
import org.gvi.core.spi.IndexResult;

import java.util.List;
import java.util.Map;

/**
 * Everything one pipeline run produced.
 * <p>
 * The headline result is {@code datasetGvi}: ONE composite GVI score for
 * the whole uploaded file, built from the indices that are already
 * dataset-wide by construction (pi, RI, mu, Re, in {@code populationIndices})
 * plus the per-sequence indices (GD/MB/dN-dS/CAI/GC) averaged across every
 * sequence in the file (in {@code datasetIndices}) -- see
 * {@code GviPipeline.aggregateAcrossSequences}.
 * <p>
 * {@code perSequenceIndices}/{@code gviPerSequence} are kept as supporting
 * detail (which specific sequence is driving the dataset result), not the
 * headline answer. {@code perSequenceIndices}' CAI/DNDS entries are
 * themselves an aggregate across all annotated genes (mean CAI, max dN/dS
 * -- see GviPipeline); the full per-gene breakdown and, for dN/dS, the
 * sliding-window site-resolution detail are carried separately here for
 * reporting.
 * <p>
 * {@code sensitivity} (Section 5.9/8's sensitivity analysis) shows, for
 * the whole-file GVI, how much each contributing index's weight would need
 * to move the composite score if perturbed +/-20% -- empty if the
 * whole-file GVI itself couldn't be computed.
 * <p>
 * {@code mlDnDsPerGene} (opt-in via {@code --ml-dnds}) is a SEPARATE,
 * dataset-wide (not per-sequence) maximum-likelihood dN/dS estimate per
 * gene -- {@code org.gvi.algorithms.dnds.ml.MlCodonDnDsEstimator}'s GY94/
 * {@code codeml} M0-equivalent fit, jointly using the whole tree rather
 * than pairwise-to-reference counting. It does not replace or feed into
 * {@code dnDsPerGene}/the composite GVI's dN/dS component (still the
 * already-validated Nei-Gojobori method); it is additional, cross-checking
 * detail, empty for genes that exceeded its taxa/codon caps.
 * <p>
 * {@code reassortment} (opt-in via 2+ {@code --segment label=path} options) is a completely separate
 * analysis from everything else here: it compares independently-aligned genome SEGMENTS (e.g.
 * Bluetongue's 10, influenza's 8) for the same isolates, one {@link ReassortmentResult} per segment
 * pair, and is not folded into the composite GVI (which is spec'd around indices from one alignment).
 * Empty unless segments were supplied.
 */
public record PipelineResult(
        Map<IndexKey, IndexResult> populationIndices,
        Map<IndexKey, IndexResult> datasetIndices,
        GviResult datasetGvi,
        Map<String, Map<IndexKey, IndexResult>> perSequenceIndices,
        Map<String, GviResult> gviPerSequence,
        Map<String, List<CaiResult>> caiPerGene,
        Map<String, List<DnDsResult>> dnDsPerGene,
        Map<String, List<SlidingWindowDnDsResult>> dnDsSlidingWindows,
        Map<String, MlDnDsResult> mlDnDsPerGene,
        List<SensitivityResult> sensitivity,
        List<String> skipped,
        List<String> warnings,
        DatasetSummary datasetSummary,
        List<ReassortmentResult> reassortment
) {
}
