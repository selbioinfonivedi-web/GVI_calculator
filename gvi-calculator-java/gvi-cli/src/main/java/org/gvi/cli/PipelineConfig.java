package org.gvi.cli;

import org.gvi.algorithms.gd.GdMethod;
import org.gvi.algorithms.mu.GenomeType;
import org.gvi.core.model.OrganismClass;

import java.nio.file.Path;
import java.util.Map;
import java.util.Set;

/** Everything {@link GviPipeline} needs for one run, gathered from CLI options. */
public record PipelineConfig(
        Path fastaPath,
        Path metadataPath,
        Path gffPath,
        Path codonUsagePath,
        String codonUsageSpecies,
        Path incidencePath,
        Double referenceGcPercent,
        String referenceId,
        Set<String> indicesRequested, // lowercase keys: mu,re,pi,mb,dnds,gd,cai,gc,ri or "all"
        OrganismClass organismClass,
        String pathogenId, // key into the bundled generation-time table; null when --generation-time-days is given directly
        boolean generationTimeExplicit, // true when the user actually passed --generation-time-days
        double generationTimeDays,
        GdMethod gdMethod,
        double beta0,
        double betaScaleFactor,
        boolean highAccuracyMu,
        Double gammaAlpha,
        String substitutionModel,
        boolean bootstrapSupport,
        boolean leastSquaresDating,
        boolean relaxedClockMu,
        boolean bdskyRe,
        boolean mlDnds,
        Path weightsPath,
        GenomeType genomeType,
        Map<String, Path> segmentFastas, // label -> path; 2+ entries opts into reassortment detection, entirely separate from fastaPath
        Path gffOut, // if set and no --gff was supplied (or it failed to parse), the native ORF-scan prediction is written here as a real GFF3
        boolean trimToCovered, // restrict the alignment to columns every sequence covers before computing anything
        int minTrimmedColumns, // refuse to trim below this many columns
        boolean autoMode, // derive everything derivable from the alignment alone, and say what was derived
        Path derivedOutputDir // where --auto writes the files it derived, so they can be inspected and reused
) {
    /**
     * Back-compatible constructor for callers predating organism-class gating and the
     * per-pathogen generation-time table. Defaults to {@link OrganismClass#UNSPECIFIED}
     * (no gating) and treats the supplied generation time as explicit, which preserves
     * the prior behaviour exactly.
     */
    public PipelineConfig(Path fastaPath, Path metadataPath, Path gffPath, Path codonUsagePath,
                          String codonUsageSpecies, Path incidencePath, Double referenceGcPercent,
                          String referenceId, Set<String> indicesRequested, double generationTimeDays,
                          GdMethod gdMethod, double beta0, double betaScaleFactor, boolean highAccuracyMu,
                          Double gammaAlpha, String substitutionModel, boolean bootstrapSupport,
                          boolean leastSquaresDating, boolean bdskyRe, boolean mlDnds, Path weightsPath,
                          GenomeType genomeType, Map<String, Path> segmentFastas, Path gffOut) {
        this(fastaPath, metadataPath, gffPath, codonUsagePath, codonUsageSpecies, incidencePath,
                referenceGcPercent, referenceId, indicesRequested, OrganismClass.UNSPECIFIED, null, true,
                generationTimeDays, gdMethod, beta0, betaScaleFactor, highAccuracyMu, gammaAlpha,
                substitutionModel, bootstrapSupport, leastSquaresDating, false, bdskyRe, mlDnds, weightsPath,
                genomeType, segmentFastas, gffOut, false, 100, false, null);
    }

    public boolean wants(String indexKey) {
        return indicesRequested.contains("all") || indicesRequested.contains(indexKey);
    }

    /** Null-safe accessor -- older callers may leave the field unset. */
    public OrganismClass organismClassOrUnspecified() {
        return organismClass == null ? OrganismClass.UNSPECIFIED : organismClass;
    }
}
