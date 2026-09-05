package org.gvi.web;

import com.fasterxml.jackson.annotation.JsonIgnoreProperties;

import java.util.List;

/**
 * One analysis request, posted as JSON.
 * <p>
 * Files arrive as raw text rather than as a multipart upload: the browser reads them with
 * {@code FileReader} and posts their contents. That avoids hand-rolling a multipart parser against the
 * JDK's bare {@code HttpServer} (which has no such support), and these inputs are single-locus alignments
 * of a few dozen sequences, not references -- {@link AnalysisService#MAX_INPUT_CHARS} caps them.
 * <p>
 * {@code organismClass} and {@code genomeType} are nullable but should not be. Neither can be inferred
 * from sequence, and both change what the numbers mean: organism class gates whether host-relative CAI
 * is biologically meaningful at all, and genome type sets which ceiling mu is normalized against (a DNA
 * rate judged on the RNA scale contributes essentially nothing regardless of how fast it is for a DNA
 * genome). The UI asks for both and warns when they are left unset.
 */
@JsonIgnoreProperties(ignoreUnknown = true)
public class AnalyzeRequest {
    public String fasta;
    public String metadata;
    public String gff;
    public String codonUsage;
    public String codonUsageSpecies;
    public String incidence;

    public String organismClass;
    public String genomeType;
    public String pathogenId;
    /**
     * Mean generation time in days. Null means "resolve it from {@link #pathogenId} against the
     * bundled table, or fall back to the default".
     * <p>
     * This is exposed because it is the single most consequential number a caller can get wrong for
     * Re. The birth-death process is scaled by {@code delta = 365.25 / generationTimeDays}, so a
     * generation time wrong by a factor of k moves Re by roughly the same factor, and Re carries the
     * largest weight in the composite.
     */
    public Double generationTimeDays;
    public Double referenceGc;
    public String referenceId;
    public String gdMethod;

    public boolean trimToCovered;
    public boolean auto;
    public boolean perSequence;
    public boolean bdskyRe = true;

    public List<String> indices;
}
