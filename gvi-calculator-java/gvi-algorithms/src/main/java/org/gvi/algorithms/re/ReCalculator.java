package org.gvi.algorithms.re;

import org.gvi.core.exception.GviComputationException;
import org.gvi.core.exception.GviInputException;
import org.gvi.core.model.IncidencePoint;
import org.gvi.core.model.SequenceAlignment;

import java.util.ArrayList;
import java.util.List;

/**
 * Facade for Index 2 (Section 5.2): picks the incidence-based Cori
 * estimator when a case time series is supplied (preferred -- it is the
 * direct, better-validated method), and falls back to the phylodynamic
 * lineage-accumulation approximation when only dated sequences are
 * available. Which estimator ran is recorded on the result.
 */
public final class ReCalculator {

    private final CoriReEstimator cori;
    private final PhylodynamicReEstimator phylo;
    private final SerialInterval serialInterval;
    private final double defaultGenerationTimeDays;

    public ReCalculator() {
        this(new CoriReEstimator(), new PhylodynamicReEstimator(), SerialInterval.defaultProfile(), 5.0);
    }

    public ReCalculator(CoriReEstimator cori, PhylodynamicReEstimator phylo, SerialInterval serialInterval, double defaultGenerationTimeDays) {
        this.cori = cori;
        this.phylo = phylo;
        this.serialInterval = serialInterval;
        this.defaultGenerationTimeDays = defaultGenerationTimeDays;
    }

    public ReResult compute(List<IncidencePoint> incidence, SequenceAlignment alignment) {
        return compute(incidence, alignment, defaultGenerationTimeDays);
    }

    public ReResult compute(List<IncidencePoint> incidence, SequenceAlignment alignment, double generationTimeDays) {
        boolean hasIncidence = incidence != null && !incidence.isEmpty();
        boolean hasAlignment = alignment != null;

        if (!hasIncidence && !hasAlignment) {
            throw new GviInputException("Re requires either a case-incidence time series or a dated sequence alignment");
        }

        if (hasIncidence) {
            try {
                return cori.compute(incidence, serialInterval);
            } catch (GviComputationException e) {
                if (!hasAlignment) throw e;
                ReResult fallback = phylo.compute(alignment, generationTimeDays);
                List<String> diagnostics = new ArrayList<>(fallback.diagnostics());
                diagnostics.add(0, "Incidence-based estimation failed (" + e.getMessage() + "); used phylodynamic fallback instead");
                return new ReResult(fallback.method(), fallback.estimateDate(), fallback.reMean(), fallback.interval(),
                        fallback.series(), fallback.doublingTimeDays(), fallback.category(), diagnostics);
            }
        }
        return phylo.compute(alignment, generationTimeDays);
    }
}
