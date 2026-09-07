package org.gvi.cli;

import org.gvi.algorithms.cai.CaiReferenceTable;
import org.gvi.algorithms.cai.CaiResult;
import org.gvi.algorithms.cai.CodonAdaptationCalculator;
import org.gvi.algorithms.cai.GcContentCalculator;
import org.gvi.algorithms.cai.GcResult;
import org.gvi.algorithms.cai.GcReferenceTable;
import org.gvi.algorithms.dnds.DnDsCalculator;
import org.gvi.algorithms.dnds.DnDsGeneReferenceTable;
import org.gvi.algorithms.dnds.DnDsResult;
import org.gvi.algorithms.dnds.SlidingWindowDnDsResult;
import org.gvi.algorithms.gd.GdReferenceTable;
import org.gvi.algorithms.gd.GdResult;
import org.gvi.algorithms.gd.GeneticDistanceCalculator;
import org.gvi.algorithms.mb.MbReferenceTable;
import org.gvi.algorithms.mb.MbResult;
import org.gvi.algorithms.mb.MutationBurdenCalculator;
import org.gvi.algorithms.mu.EvolutionaryRateCalculator;
import org.gvi.algorithms.mu.MuResult;
import org.gvi.algorithms.phylo.BootstrapSupportCalculator;
import org.gvi.algorithms.pi.NucleotideDiversityCalculator;
import org.gvi.algorithms.pi.PiResult;
import org.gvi.algorithms.re.ReCalculator;
import org.gvi.algorithms.re.ReResult;
import org.gvi.algorithms.ri.ReassortmentIndexCalculator;
import org.gvi.algorithms.ri.ReassortmentResult;
import org.gvi.algorithms.ri.RecombinationIndexCalculator;
import org.gvi.algorithms.ri.RiResult;
import org.gvi.composite.CompositeGviEngine;
import org.gvi.composite.CompositeWeights;
import org.gvi.composite.GviResult;
import org.gvi.composite.IndexKey;
import org.gvi.composite.NormalizationRange;
import org.gvi.composite.SensitivityResult;
import org.gvi.core.exception.GviException;
import org.gvi.core.io.*;
import org.gvi.core.model.*;
import org.gvi.core.spi.IndexResult;
import org.gvi.core.spi.SimpleIndexResult;

import java.util.*;

/**
 * Wires the 8 index modules and the composite engine together against
 * whatever inputs the user actually supplied (Section 6: each index
 * independently computable AND combinable). Every optional computation is
 * individually try/caught -- one missing input or one pathological
 * sequence never aborts the whole run, it just gets recorded in
 * {@link PipelineResult#skipped()}.
 */
public final class GviPipeline {

    public PipelineResult run(PipelineConfig config) {
        List<String> skipped = new ArrayList<>();
        List<String> warnings = new ArrayList<>();

        FastaReader.Result fastaResult = FastaReader.read(config.fastaPath());
        collectWarnings(warnings, "fasta", fastaResult.report());

        List<NucleotideSequence> sequences = fastaResult.sequences();
        if (config.metadataPath() != null) {
            MetadataCsvReader.Result metaResult = MetadataCsvReader.read(config.metadataPath());
            collectWarnings(warnings, "metadata", metaResult.report());
            MetadataJoiner.Result joined = MetadataJoiner.join(sequences, metaResult.rows());
            collectWarnings(warnings, "metadata-join", joined.report());
            sequences = joined.sequences();
        }

        SequenceAlignment alignment = SequenceAlignment.of(sequences, config.referenceId());

        // Two cheap checks on the alignment itself, before anything is computed from it. Both come
        // from corpus datasets that produced a confident-looking GVI from inputs that could not
        // support one -- the diagnosis existed, but only scattered across per-index exclusions after
        // the whole run. Stating it once up front, in the alignment's own terms, is a better place
        // to learn it. Warnings, not gates: a deliberately divergent panel is a legitimate input.
        org.gvi.core.util.AlignmentPreflight.Report preflight =
                org.gvi.core.util.AlignmentPreflight.check(alignment);
        for (String finding : preflight.findings()) {
            warnings.add("Input: " + finding);
        }

        // Optional trim to the shared coverage window. Must happen before any index runs, since the
        // whole point is that the discarded columns are absent data rather than observed difference.
        org.gvi.core.util.AlignmentTrimmer.Result trim = null;
        // --auto trims only when the alignment would otherwise fail the coverage gate. Trimming a
        // well-covered alignment would discard real data for nothing, so the decision is made from
        // the same gap-fraction threshold the quality gate uses rather than applied unconditionally.
        boolean autoTrim = config.autoMode()
                && QualityGate.gapFraction(alignment) >= QualityGate.GAP_FRACTION_LIMIT;
        if (config.trimToCovered() || autoTrim) {
            trim = org.gvi.core.util.AlignmentTrimmer.trimToFullyCovered(alignment, config.minTrimmedColumns());
            warnings.add((autoTrim && !config.trimToCovered() ? "--auto (coverage gate would fail): " : "--trim-to-covered: ")
                    + trim.description());
            if (trim.trimmed()) alignment = trim.alignment(); else trim = null;
        }

        NucleotideSequence reference = alignment.getReference();

        Map<IndexKey, IndexResult> population = new EnumMap<>(IndexKey.class);
        Map<String, Map<IndexKey, IndexResult>> perSequence = new LinkedHashMap<>();
        for (NucleotideSequence q : alignment.getQueries()) {
            perSequence.put(q.getId(), new EnumMap<>(IndexKey.class));
        }

        computePi(config, alignment, population, skipped);
        computeRi(config, alignment, population, skipped);
        List<ReassortmentResult> reassortment = computeReassortment(config, warnings, skipped);
        computeMu(config, alignment, population, skipped, warnings);
        computeRe(config, alignment, population, skipped, warnings);
        computeGd(config, alignment, reference, perSequence, skipped);
        computeMb(config, alignment, reference, perSequence, skipped);

        GeneLoadResult geneLoad = loadGenes(config, alignment, warnings);
        List<GeneAnnotation> genes = remapGenesAfterTrim(geneLoad.genes(), trim, warnings);
        writeGffOut(config, geneLoad, reference, warnings);
        Map<String, List<CaiResult>> caiPerGene = new LinkedHashMap<>();
        Map<String, List<DnDsResult>> dnDsPerGene = new LinkedHashMap<>();
        Map<String, List<SlidingWindowDnDsResult>> dnDsSlidingWindows = new LinkedHashMap<>();
        Map<String, DnDsResult> winningDnDsPerSequence = new LinkedHashMap<>();
        computeCai(config, alignment, genes, perSequence, caiPerGene, skipped, warnings);
        computeGc(config, alignment, reference, perSequence, skipped, warnings);
        computeDnDs(config, alignment, reference, genes, perSequence, dnDsPerGene, dnDsSlidingWindows, winningDnDsPerSequence, skipped);
        Map<String, org.gvi.algorithms.dnds.ml.MlDnDsResult> mlDnDsPerGene = computeMlDnDs(config, alignment, genes, warnings, skipped);

        CompositeGviEngine engine = new CompositeGviEngine();
        CompositeWeights weights = loadWeights(config, warnings);

        Aggregated aggregated = aggregateAcrossSequences(perSequence, winningDnDsPerSequence, config.organismClassOrUnspecified());
        Map<IndexKey, IndexResult> datasetIndices = aggregated.indices();
        Map<IndexKey, IndexResult> datasetCombined = new EnumMap<>(datasetIndices);
        datasetCombined.putAll(population);

        // Quality gating runs BEFORE the composite, so an index whose own diagnostics say its
        // inputs are invalid is dropped from the weighted score instead of contributing to it at
        // full weight. The engine's existing missing-index weight renormalization does the rest.
        // Excluded indices are still computed and still reported in their own sections.
        QualityGate gate = QualityGate.evaluate(alignment, population, datasetIndices, aggregated.dnDsSelection(),
                config.genomeType());
        for (QualityGate.Exclusion exclusion : gate.exclusions()) {
            skipped.add(exclusion.key().label() + " (excluded from composite GVI): " + exclusion.reason());
        }
        Map<IndexKey, IndexResult> compositeInputs = gate.filter(datasetCombined);
        Map<IndexKey, NormalizationRange> ranges = normalizationRanges(config);

        GviResult datasetGvi = null;
        List<SensitivityResult> sensitivity = List.of();
        try {
            datasetGvi = engine.compute(compositeInputs, weights, ranges);
            sensitivity = engine.sensitivityAnalysis(compositeInputs, weights, ranges, SENSITIVITY_PERTURBATION_FRACTION);
            addCoverageAdvisory(datasetGvi, warnings);
        } catch (GviException e) {
            warnings.add("Whole-file GVI: " + e.getMessage());
        }

        Map<String, GviResult> gviPerSequence = new LinkedHashMap<>();
        for (var entry : perSequence.entrySet()) {
            Map<IndexKey, IndexResult> combined = new EnumMap<>(entry.getValue());
            combined.putAll(population);
            try {
                gviPerSequence.put(entry.getKey(), engine.compute(gate.filter(combined), weights, ranges));
            } catch (GviException e) {
                warnings.add("GVI for '" + entry.getKey() + "': " + e.getMessage());
            }
        }

        DatasetSummary datasetSummary = buildDatasetSummary(alignment);
        addDataQualityAdvisories(alignment, population, datasetIndices, warnings, geneLoad.autoPredicted());
        addLineageAdvisory(alignment, population, gate, warnings);
        writeDerivedInputs(config, alignment, genes, geneLoad.autoPredicted(), warnings);

        return new PipelineResult(population, datasetIndices, datasetGvi, perSequence, gviPerSequence,
                caiPerGene, dnDsPerGene, dnDsSlidingWindows, mlDnDsPerGene, sensitivity, skipped, warnings,
                datasetSummary, reassortment);
    }

    /**
     * Normalization ranges for this run. Only mu is genome-type dependent: a DNA genome's
     * substitution rate measured against the RNA-scaled ceiling normalizes to ~0 no matter how
     * fast it is for a DNA genome, which is the same split MuReferenceTable already applies to
     * the category labels.
     */
    /**
     * Warns when too little of the weighting scheme survived gating for the headline GVI to be compared
     * against another dataset's. See {@link GviResult#MIN_COMPARABLE_COVERAGE} for why a low-coverage score
     * is actively misleading rather than merely imprecise: the engine renormalizes the surviving weights to
     * sum to 1, so a GVI built from two indices still lands in [0,1] and still looks like a GVI.
     */
    private void addCoverageAdvisory(GviResult gvi, List<String> warnings) {
        if (gvi.comparable()) return;
        warnings.add(String.format(java.util.Locale.ROOT,
                "NOT COMPARABLE: this GVI of %.4f was computed from only %s (threshold for cross-dataset comparison "
                        + "is %.0f%%). The surviving indices were renormalized to sum to 1, so the score still falls in "
                        + "[0,1] and still looks like a normal GVI -- but it summarizes a different, much smaller set of "
                        + "evidence than a fully-populated score does. Do NOT rank this value against other pathogens' "
                        + "GVIs. A LOW value here is especially easy to misread: it can mean the contributing indices "
                        + "happened to be small, not that virulence is low. Resolve the exclusions listed under "
                        + "'Skipped' and re-run before drawing any comparative conclusion.",
                gvi.gvi(), gvi.coverageSummary(), GviResult.MIN_COMPARABLE_COVERAGE * 100));
    }

    private Map<IndexKey, NormalizationRange> normalizationRanges(PipelineConfig config) {
        double muMax = config.genomeType() == org.gvi.algorithms.mu.GenomeType.DNA
                ? NormalizationRange.MU_MAX_DNA
                : NormalizationRange.MU_MAX_RNA;
        return NormalizationRange.defaults(muMax);
    }

    /**
     * When mu was rejected for want of temporal signal, the documented remedy is to split the
     * alignment by lineage and estimate a rate within each -- but that advice is unactionable
     * without knowing where the split falls. This reports the partition, so the user has something
     * concrete to act on instead of a general instruction.
     */
    private void addLineageAdvisory(SequenceAlignment alignment, Map<IndexKey, IndexResult> population,
                                    QualityGate gate, List<String> warnings) {
        if (!gate.isExcluded(IndexKey.MU)) return;
        if (!(population.get(IndexKey.MU) instanceof MuResult)) return;
        var partition = org.gvi.core.util.LineagePartitioner.detect(alignment);
        if (partition.partitioned()) {
            warnings.add("mu remedy -- " + partition.description());
        }
    }

    /**
     * Rewrites gene coordinates onto the trimmed alignment's numbering.
     * <p>
     * Trimming renumbers every column, so a GFF3 supplied against the original alignment no longer
     * describes it -- a CDS at [1, 1146] against a 456-column trimmed alignment simply fails to
     * load, costing exactly the dN/dS and CAI that trimming was supposed to rescue. A gene whose
     * span is entirely discarded is dropped with an explanation rather than silently mangled.
     */
    private List<GeneAnnotation> remapGenesAfterTrim(List<GeneAnnotation> genes,
                                                     org.gvi.core.util.AlignmentTrimmer.Result trim,
                                                     List<String> warnings) {
        if (trim == null || !trim.trimmed() || genes.isEmpty()) return genes;

        List<GeneAnnotation> remapped = new ArrayList<>(genes.size());
        List<String> dropped = new ArrayList<>();
        for (GeneAnnotation g : genes) {
            var mapped = trim.remap(g.start(), g.end());
            if (mapped.isEmpty()) {
                dropped.add(g.geneName());
                continue;
            }
            long[] se = mapped.get();
            remapped.add(new GeneAnnotation(g.geneName(), se[0], se[1], g.strand()));
        }

        warnings.add(String.format(
                "--trim-to-covered: remapped %d gene annotation(s) onto the trimmed coordinates%s. Gene spans supplied "
                        + "against the untrimmed alignment would otherwise fall outside it and be skipped, silently "
                        + "costing the dN/dS and CAI the trim was meant to recover.",
                remapped.size(),
                dropped.isEmpty() ? "" : "; dropped " + dropped.size() + " gene(s) whose span was entirely discarded ("
                        + String.join(", ", dropped) + ")"));
        return remapped;
    }

    /**
     * Writes out what {@code --auto} derived from the alignment alone, so the user can see it,
     * correct it, and feed it back in on a later run rather than re-deriving every time.
     * <p>
     * Derivation that stays invisible is hard to trust: a collection date parsed out of a FASTA
     * header is a real inference about the data, and the analyst should be able to check it against
     * the accession rather than take it on faith.
     */
    private void writeDerivedInputs(PipelineConfig config, SequenceAlignment alignment,
                                    List<GeneAnnotation> genes, boolean genesAutoPredicted,
                                    List<String> warnings) {
        if (!config.autoMode() || config.derivedOutputDir() == null) return;

        java.nio.file.Path dir = config.derivedOutputDir();
        String stem = config.fastaPath().getFileName().toString().replaceAll("\\.[^.]*$", "");
        List<String> written = new ArrayList<>();

        // Metadata parsed out of the FASTA headers.
        int dated = 0;
        StringBuilder csv = new StringBuilder("sequence_id,collection_date,location,host\n");
        for (NucleotideSequence s : alignment.getSequences()) {
            if (s.getCollectionDate().isPresent()) dated++;
            csv.append(s.getId()).append(',')
               .append(s.getCollectionDate().map(Object::toString).orElse("")).append(',')
               .append(s.getLocation().orElse("")).append(',')
               .append(s.getHost().orElse("")).append('\n');
        }
        java.nio.file.Path metaPath = dir.resolve(stem + ".derived-metadata.csv");
        try {
            java.nio.file.Files.writeString(metaPath, csv.toString());
            written.add(metaPath.getFileName() + " (" + dated + "/" + alignment.size() + " sequences carried a parseable date)");
        } catch (java.io.IOException e) {
            warnings.add("--auto: could not write derived metadata to " + metaPath + ": " + e.getMessage());
        }

        // Gene boundaries, but only when they were predicted rather than supplied.
        if (genesAutoPredicted && !genes.isEmpty()) {
            StringBuilder gff = new StringBuilder("##gff-version 3\n");
            String seqId = alignment.getReference().getId();
            for (GeneAnnotation g : genes) {
                gff.append(seqId).append("\tgvi-auto\tCDS\t").append(g.start()).append('\t').append(g.end())
                   .append("\t.\t").append(g.strand()).append("\t0\tID=").append(g.geneName())
                   .append(";Name=").append(g.geneName()).append('\n');
            }
            java.nio.file.Path gffPath = dir.resolve(stem + ".predicted-genes.gff3");
            try {
                java.nio.file.Files.writeString(gffPath, gff.toString());
                written.add(gffPath.getFileName() + " (" + genes.size() + " predicted ORF(s))");
            } catch (java.io.IOException e) {
                warnings.add("--auto: could not write predicted genes to " + gffPath + ": " + e.getMessage());
            }
        }

        if (!written.isEmpty()) {
            warnings.add("--auto: derived inputs written for inspection -- " + String.join("; ", written)
                    + ". Supply these with --metadata/--gff on a later run to skip re-deriving them, or correct them first.");
        }
    }

    /** Gap fraction above which an alignment is unusually gap-heavy for a normal viral MSA. */
    private static final double HIGH_GAP_FRACTION_THRESHOLD = 0.15;
    /** Mutation-burden fraction of alignment length above which divergence looks implausibly high. */
    private static final double HIGH_MB_FRACTION_THRESHOLD = 0.15;
    /** dN/dS above which the whole-file estimate is worth a second look. */
    private static final double NOTABLE_DNDS_THRESHOLD = 1.5;
    /** A sequence pair's real-base (non-gap) spans overlapping less than this fraction of the shorter one's own span counts as "poorly overlapping" -- a signature of two different primer windows padded to the same alignment length rather than a real shared amplicon. */
    private static final double MIN_PAIRWISE_SPAN_OVERLAP = 0.5;
    /** Fraction of all sequence pairs that must be poorly-overlapping before it's flagged as a dataset-wide pattern rather than one or two odd sequences. */
    private static final double POOR_OVERLAP_PAIR_FRACTION_THRESHOLD = 0.15;

    /**
     * Holistic, general-purpose sanity checks over the FINISHED result -- distinct from the per-step
     * try/catch warnings collected while computing individual indices. These don't require knowing
     * anything about the specific virus/species (this tool is offline and has no reference database of
     * "expected" values per organism) -- they flag internally-inconsistent patterns that any correctly
     * prepared viral alignment should not produce, and point at the most common real-world cause:
     * multiple genetically distinct lineages pooled into one alignment, or the input not actually being
     * aligned despite equal sequence lengths.
     */
    private void addDataQualityAdvisories(SequenceAlignment alignment, Map<IndexKey, IndexResult> population,
                                           Map<IndexKey, IndexResult> datasetIndices, List<String> warnings,
                                           boolean genesAutoPredicted) {
        IndexResult muResult = population.get(IndexKey.MU);
        if (muResult instanceof org.gvi.algorithms.mu.MuResult mu && (mu.muSubPerSiteYear() < 0 || !mu.reliable())) {
            warnings.add(String.format(java.util.Locale.ROOT,
                    "Data quality advisory -- mu: the evolutionary rate came out %s (R²=%.2f, reliable threshold 0.30). "
                            + "A single molecular-clock rate is only meaningful when every sequence shares one continuously "
                            + "evolving lineage. Common causes: collection dates don't span enough real time, or the "
                            + "alignment pools multiple genetically distinct lineages/strains/serotypes that were "
                            + "introduced independently rather than evolving from each other. If your sequences span more "
                            + "than one lineage, split them into separate per-lineage alignments and re-run each separately.",
                    mu.muSubPerSiteYear() < 0 ? "negative" : "unreliable", mu.rSquared()));
        }

        long gapChars = 0;
        long totalChars = 0;
        for (NucleotideSequence seq : alignment.getSequences()) {
            String s = seq.getSequence();
            totalChars += s.length();
            for (int i = 0; i < s.length(); i++) {
                char c = s.charAt(i);
                if (c == '-' || c == '.' || c == '?') gapChars++;
            }
        }
        double gapFraction = totalChars > 0 ? (double) gapChars / totalChars : 0;
        if (gapFraction >= HIGH_GAP_FRACTION_THRESHOLD) {
            warnings.add(String.format(java.util.Locale.ROOT,
                    "Data quality advisory -- alignment: %.0f%% of the alignment is gap characters (-/./?), which is "
                            + "unusually high for a viral multi-FASTA. This tool requires an ALREADY-ALIGNED input; a high "
                            + "gap fraction usually means the sequences were padded/concatenated to equal length rather "
                            + "than genuinely aligned. Re-align with a real multiple-sequence-alignment tool before "
                            + "re-running, e.g. \"mafft --auto input.fasta > aligned.fasta\" or \"muscle -align "
                            + "input.fasta -output aligned.fasta\".",
                    gapFraction * 100));
        }

        checkAmpliconWindowOverlap(alignment, warnings);

        IndexResult mbResult = datasetIndices.get(IndexKey.MB);
        if (mbResult != null && alignment.length() > 0) {
            double mbFraction = mbResult.primaryValue() / alignment.length();
            if (mbFraction >= HIGH_MB_FRACTION_THRESHOLD) {
                warnings.add(String.format(java.util.Locale.ROOT,
                        "Data quality advisory -- mutation burden: the average %.1f variant events per sequence is %.0f%% "
                                + "of the %d bp alignment length, implausibly high for sequences assumed to be closely "
                                + "related. Likely causes: a misalignment (see the alignment-quality advisory above if "
                                + "present), the wrong reference sequence selected (check Reference sequence id), or the "
                                + "sequences genuinely aren't as closely related as assumed.",
                        mbResult.primaryValue(), mbFraction * 100, alignment.length()));
            }
        }

        IndexResult dndsResult = datasetIndices.get(IndexKey.DNDS);
        if (dndsResult != null && genesAutoPredicted) {
            // Fires regardless of the resulting value -- unlike the notable-value check below, an auto-predicted
            // ORF means the premise itself (this marker is protein-coding) was never confirmed, and a
            // plausible-looking dN/dS from a chance ORF in a non-coding marker (e.g. rRNA) is just as wrong as an
            // extreme one; a value that happens to look reasonable was previously never flagged at all.
            warnings.add(String.format(java.util.Locale.ROOT,
                    "Data quality advisory -- dN/dS/CAI premise: no --gff was supplied, so gene boundaries for this "
                            + "dN/dS (%.3f) and any CAI result were auto-predicted by a 6-frame ORF scan, not confirmed as "
                            + "real coding sequence. That scan finds SOME ATG...stop run by chance in most alignments, "
                            + "including genuinely non-coding markers (rRNA, ITS, intergenic regions) -- this caveat applies "
                            + "even when the number itself looks biologically plausible, not just when it looks extreme. If "
                            + "this marker is not protein-coding, exclude dnds/cai via --indices instead of trusting this "
                            + "value; if it is, supply --gff for real gene coordinates.",
                    dndsResult.primaryValue()));
        } else if (dndsResult != null && dndsResult.primaryValue() >= NOTABLE_DNDS_THRESHOLD) {
            warnings.add(String.format(java.util.Locale.ROOT,
                    "Data quality advisory -- dN/dS: the whole-file estimate (%.3f) is high enough to suggest either "
                            + "strong positive selection or, more commonly, a misalignment/frameshift artifact (a shifted "
                            + "reading frame turns synonymous change into apparent nonsynonymous change). If this gene is "
                            + "expected to be conserved, double-check the coding frame and gene boundaries (--gff) before "
                            + "trusting this value.",
                    dndsResult.primaryValue()));
        }
    }

    /**
     * Catches a real-world failure mode the gap-fraction check above (deliberately) doesn't: sequences
     * pulled from different studies that used different, non-overlapping primer pairs for "the same gene."
     * MAFFT still aligns them -- large internal gap blocks pad each sequence out to the shared alignment
     * length -- and the result can stay under the general gap-fraction threshold even though specific pairs
     * of sequences share almost no real (non-gap) sequence in common. Any pairwise diversity/distance/dN/dS
     * computed between such a pair is measuring amplicon-boundary artifact, not biology (confirmed this is
     * a real, recurring pattern in independent real-world testing -- see the project's validation notes).
     * Flags the dataset as a whole when a substantial fraction of all pairs are affected, rather than
     * treating one or two short/partial sequences as automatically suspicious.
     */
    private void checkAmpliconWindowOverlap(SequenceAlignment alignment, List<String> warnings) {
        List<NucleotideSequence> seqs = alignment.getSequences();
        int n = seqs.size();
        if (n < 3) return; // need enough pairs for "a substantial fraction" to mean anything

        int[] firstBase = new int[n];
        int[] lastBase = new int[n];
        for (int i = 0; i < n; i++) {
            String s = seqs.get(i).getSequence();
            int first = -1, last = -1;
            for (int c = 0; c < s.length(); c++) {
                char ch = s.charAt(c);
                if (ch != '-' && ch != '.' && ch != '?') {
                    if (first == -1) first = c;
                    last = c;
                }
            }
            firstBase[i] = first;
            lastBase[i] = last;
        }

        long poorPairs = 0, totalPairs = 0;
        for (int i = 0; i < n; i++) {
            if (firstBase[i] == -1) continue; // all-gap sequence -- not this check's concern
            for (int j = i + 1; j < n; j++) {
                if (firstBase[j] == -1) continue;
                totalPairs++;
                int overlap = Math.max(0, Math.min(lastBase[i], lastBase[j]) - Math.max(firstBase[i], firstBase[j]) + 1);
                int minSpan = Math.min(lastBase[i] - firstBase[i] + 1, lastBase[j] - firstBase[j] + 1);
                if (minSpan > 0 && (double) overlap / minSpan < MIN_PAIRWISE_SPAN_OVERLAP) {
                    poorPairs++;
                }
            }
        }
        if (totalPairs == 0) return;
        double poorFraction = (double) poorPairs / totalPairs;
        if (poorFraction < POOR_OVERLAP_PAIR_FRACTION_THRESHOLD) return;

        int refIdx = seqs.indexOf(alignment.getReference());
        List<String> examples = new ArrayList<>();
        if (refIdx >= 0 && firstBase[refIdx] != -1) {
            for (int i = 0; i < n && examples.size() < 3; i++) {
                if (i == refIdx || firstBase[i] == -1) continue;
                int overlap = Math.max(0, Math.min(lastBase[refIdx], lastBase[i]) - Math.max(firstBase[refIdx], firstBase[i]) + 1);
                int minSpan = Math.min(lastBase[refIdx] - firstBase[refIdx] + 1, lastBase[i] - firstBase[i] + 1);
                if (minSpan > 0 && (double) overlap / minSpan < MIN_PAIRWISE_SPAN_OVERLAP) {
                    examples.add(seqs.get(i).getId());
                }
            }
        }

        warnings.add(String.format(java.util.Locale.ROOT,
                "Data quality advisory -- sequence span overlap: %.0f%% of sequence pairs share less than %.0f%% overlap "
                        + "in their own real-base (non-gap) span, even though the overall alignment gap fraction may look "
                        + "normal. This is the signature of pooling sequences amplified with different, non-overlapping "
                        + "primer pairs under one gene name -- the alignment tool pads each one out to a common length with "
                        + "gaps rather than genuinely overlapping them. Any diversity/distance/dN/dS computed between a "
                        + "poorly-overlapping pair reflects amplicon-boundary artifact, not real biology.%s Consider "
                        + "restricting the dataset to sequences that share one consistent amplicon window before re-running.",
                poorFraction * 100, MIN_PAIRWISE_SPAN_OVERLAP * 100,
                examples.isEmpty() ? "" : " Example low-overlap-with-reference sequence(s): " + String.join(", ", examples) + "."));
    }

    private DatasetSummary buildDatasetSummary(SequenceAlignment alignment) {
        java.time.LocalDate earliest = null;
        java.time.LocalDate latest = null;
        for (NucleotideSequence seq : alignment.getSequences()) {
            var date = seq.getCollectionDate();
            if (date.isEmpty()) continue;
            if (earliest == null || date.get().isBefore(earliest)) earliest = date.get();
            if (latest == null || date.get().isAfter(latest)) latest = date.get();
        }
        return new DatasetSummary(alignment.size(), alignment.length(), alignment.getReference().getId(),
                earliest, latest);
    }

    /** Section 5.9/8's default perturbation for sensitivity analysis: how much each index's weight is nudged +/-. */
    private static final double SENSITIVITY_PERTURBATION_FRACTION = 0.2;

    /**
     * Aggregates the per-sequence indices (GD, MB, dN/dS, CAI, GC -- each
     * computed per query sequence against the reference) into ONE
     * dataset-wide value per index. GD, MB, CAI, and GC are combined via
     * the mean across all sequences that have them -- a legitimate,
     * standard population-level statistic for each (mean divergence from
     * reference, mean mutation burden, mean codon adaptation, mean GC
     * deviation -- not an approximation of something more "real").
     * dN/dS is handled separately in {@link #run}: it is POOLED (raw
     * Nei-Gojobori counts summed, then one ratio computed from the totals
     * -- see {@link DnDsCalculator#pool}), not averaged, because averaging
     * individual omega ratios is statistically biased. Combined with the
     * indices that are already dataset-wide by construction (pi, RI, mu,
     * Re -- computed once from the whole alignment/tree, never averaged
     * either), this is what the whole-file GVI is built from.
     */
    private Aggregated aggregateAcrossSequences(Map<String, Map<IndexKey, IndexResult>> perSequence,
                                                                  Map<String, DnDsResult> winningDnDsPerSequence,
                                                                  OrganismClass organismClass) {
        Map<IndexKey, IndexResult> dataset = new EnumMap<>(IndexKey.class);
        org.gvi.algorithms.dnds.NeiGojoboriSelectionTest.Result dnDsSelection = null;
        for (IndexKey key : new IndexKey[]{IndexKey.GD, IndexKey.MB, IndexKey.CAI, IndexKey.GC}) {
            List<Double> values = new ArrayList<>();
            for (Map<IndexKey, IndexResult> seqIndices : perSequence.values()) {
                IndexResult r = seqIndices.get(key);
                if (r != null) values.add(r.primaryValue());
            }
            if (values.isEmpty()) continue;
            double mean = values.stream().mapToDouble(Double::doubleValue).average().orElse(0.0);
            String category = retargetHgtText(categorize(key, mean), organismClass)
                    + " -- mean across " + values.size() + " sequence(s)";
            List<String> diagnostics = List.of("Mean of " + values.size() + " per-sequence value(s); range ["
                    + String.format("%.4g", java.util.Collections.min(values)) + ", "
                    + String.format("%.4g", java.util.Collections.max(values)) + "]");
            dataset.put(key, new SimpleIndexResult(key.label(), mean, category, diagnostics));
        }

        if (!winningDnDsPerSequence.isEmpty()) {
            DnDsResult pooled = new DnDsCalculator().pool("dataset", "pooled", new ArrayList<>(winningDnDsPerSequence.values()));
            List<String> diagnostics = new ArrayList<>(pooled.diagnostics());
            diagnostics.add(0, "Pooled from the winning (max-omega) gene of each of " + winningDnDsPerSequence.size() + " sequence(s)");

            // Attach the Nei-Gojobori Z-test to the pooled counts. The test was already implemented
            // and reported, but nothing consumed its verdict -- so an omega of 2.5 that is not
            // statistically distinguishable from 1 was still presented as positive selection AND
            // still carried full weight in the composite. Keeping the DnDsResult (rather than
            // flattening straight to a SimpleIndexResult) preserves the raw counts the test needs.
            var selection = new org.gvi.algorithms.dnds.NeiGojoboriSelectionTest().test(pooled);
            diagnostics.add(String.format(java.util.Locale.ROOT,
                    "Nei-Gojobori Z-test on the pooled counts: Z=%.3f, two-tailed p=%.4g, verdict=%s. %s",
                    selection.zScore(), selection.pValueTwoTailed(), selection.verdict(), selection.note()));

            // Kept as a SimpleIndexResult so the serialized shape of dataset_indices is unchanged;
            // the significance verdict travels alongside the map rather than inside it.
            dataset.put(IndexKey.DNDS, new SimpleIndexResult("dN/dS", pooled.omega(), pooled.category(), diagnostics));
            dnDsSelection = selection;
        }
        return new Aggregated(dataset, dnDsSelection);
    }

    /**
     * Dataset-level indices plus the pooled dN/dS significance verdict. The verdict is carried
     * separately rather than embedded in the {@link IndexResult} because the report serializes
     * that map directly, and widening the value type would change the published JSON shape.
     */
    record Aggregated(Map<IndexKey, IndexResult> indices,
                      org.gvi.algorithms.dnds.NeiGojoboriSelectionTest.Result dnDsSelection) {
    }

    private String categorize(IndexKey key, double meanValue) {
        return switch (key) {
            case GD -> {
                var band = GdReferenceTable.classify(meanValue);
                yield band.classification() + " (" + band.timeSinceDivergence() + ", confidence: " + band.confidence() + ")";
            }
            case MB -> {
                var band = MbReferenceTable.classify((int) Math.round(meanValue));
                yield band.pathogenExample() + " (" + band.timeline() + "); " + band.controlImplication();
            }
            case CAI -> {
                var band = org.gvi.algorithms.cai.CaiReferenceTable.classify(meanValue);
                yield band.status() + " (" + band.interpretation() + "; " + band.evolutionaryAge() + ")";
            }
            case GC -> {
                var band = GcReferenceTable.classify(meanValue);
                yield band.status() + " (" + band.interpretation() + "; " + band.evolutionaryAge() + ")";
            }
            default -> "N/A";
        };
    }

    /**
     * Strips HGT wording from a GC category for organisms that cannot acquire sequence that way.
     * Applied at BOTH the per-sequence and dataset-aggregate level -- the aggregate rebuilds its
     * own category string from the reference table, so rewriting only the per-sequence results
     * would leave the headline dataset category still asserting HGT.
     */
    private static String retargetHgtText(String category, OrganismClass organismClass) {
        if (category == null || organismClass == OrganismClass.UNSPECIFIED || organismClass.hgtCompetent()) {
            return category;
        }
        return category
                .replace("HGT suspected", "compositional shift (not HGT for this organism class)")
                .replace("No HGT signal", "No compositional anomaly")
                .replace("Minor recombination / reassortment", "Minor compositional drift");
    }

    private void computePi(PipelineConfig config, SequenceAlignment alignment, Map<IndexKey, IndexResult> population, List<String> skipped) {
        if (!config.wants("pi")) return;
        try {
            PiResult r = new NucleotideDiversityCalculator().compute(alignment);
            population.put(IndexKey.PI, r);
        } catch (GviException e) {
            skipped.add("pi: " + e.getMessage());
        }
    }

    private void computeRi(PipelineConfig config, SequenceAlignment alignment, Map<IndexKey, IndexResult> population, List<String> skipped) {
        if (!config.wants("ri")) return;
        try {
            RiResult r = new RecombinationIndexCalculator().computeFromAlignment(alignment);
            population.put(IndexKey.RI, r);
        } catch (GviException e) {
            skipped.add("RI: " + e.getMessage());
        }
    }

    /**
     * Reassortment, unlike within-locus recombination (PHI test above), is a distinct biological process
     * (segments of a segmented-genome pathogen swapping between co-infecting strains) that requires
     * comparing multiple INDEPENDENTLY-aligned segment files for the same isolates -- opt-in via 2+
     * {@code --segment label=path} options, entirely separate from {@code --fasta} (the primary
     * single-locus alignment every other index runs against). Not folded into the composite GVI (which is
     * spec'd around single-alignment indices); reported as its own section. Every segment pair is
     * independently try/caught so one bad segment file doesn't lose the others.
     */
    private List<ReassortmentResult> computeReassortment(PipelineConfig config, List<String> warnings, List<String> skipped) {
        Map<String, java.nio.file.Path> segmentPaths = config.segmentFastas();
        if (segmentPaths == null || segmentPaths.size() < 2) return List.of();

        Map<String, SequenceAlignment> segments = new LinkedHashMap<>();
        for (var entry : segmentPaths.entrySet()) {
            try {
                FastaReader.Result r = FastaReader.read(entry.getValue());
                collectWarnings(warnings, "segment '" + entry.getKey() + "'", r.report());
                segments.put(entry.getKey(), SequenceAlignment.of(r.sequences()));
            } catch (GviException e) {
                skipped.add("Reassortment segment '" + entry.getKey() + "' (" + entry.getValue() + "): " + e.getMessage());
            }
        }

        List<String> labels = new ArrayList<>(segments.keySet());
        List<ReassortmentResult> results = new ArrayList<>();
        ReassortmentIndexCalculator calc = new ReassortmentIndexCalculator();
        for (int i = 0; i < labels.size(); i++) {
            for (int j = i + 1; j < labels.size(); j++) {
                String labelA = labels.get(i);
                String labelB = labels.get(j);
                try {
                    results.add(calc.compare(labelA, segments.get(labelA), labelB, segments.get(labelB)));
                } catch (GviException e) {
                    skipped.add("Reassortment '" + labelA + "' vs '" + labelB + "': " + e.getMessage());
                }
            }
        }
        return results;
    }

    private void computeMu(PipelineConfig config, SequenceAlignment alignment, Map<IndexKey, IndexResult> population,
                           List<String> skipped, List<String> warnings) {
        if (!config.wants("mu")) return;
        EvolutionaryRateCalculator calc = new EvolutionaryRateCalculator();

        if (config.highAccuracyMu()) {
            try {
                // Joint ML branch-length + substitution-model-parameter optimization (optionally across all 6
                // named models via --substitution-model auto): slow but a real accuracy upgrade over plain
                // JC69+NJ. Opt-in via --high-accuracy-mu; falls through to the faster estimators below on failure
                // (e.g. dataset exceeds the ML size cap) rather than failing the whole run.
                MuResult r = calc.computeGtrGammaAware(alignment, config.substitutionModel(), config.gammaAlpha(), config.genomeType());
                population.put(IndexKey.MU, r);
                return;
            } catch (GviException highAccuracyFailure) {
                warnings.add("mu: --high-accuracy-mu unavailable (" + highAccuracyFailure.getMessage() + "); falling back to the standard estimator");
            }
        }

        if (config.leastSquaresDating()) {
            try {
                // Joint least-squares fit of mu AND every internal node's date directly against every edge
                // (the same objective LSD2/IQ-TREE --date solves) -- opt-in via --lsd-mu; needs every sequence
                // dated, falls through to the standard estimator below on failure rather than failing the run.
                MuResult r = calc.computeLeastSquaresDating(alignment, config.genomeType());
                population.put(IndexKey.MU, r);
                return;
            } catch (GviException lsdFailure) {
                warnings.add("mu: --lsd-mu unavailable (" + lsdFailure.getMessage() + "); falling back to the standard estimator");
            }
        }

        try {
            // Prefer the tree-aware estimator (accounts for shared ancestry between sequences);
            // fall back to pairwise-to-reference if the tree can't be built (e.g. too many taxa, too few sequences).
            MuResult r = calc.computeTreeAware(alignment, config.bootstrapSupport(), BootstrapSupportCalculator.DEFAULT_REPLICATES, config.genomeType());
            population.put(IndexKey.MU, r);
        } catch (GviException treeFailure) {
            try {
                MuResult r = calc.compute(alignment, config.genomeType());
                population.put(IndexKey.MU, r);
                warnings.add("mu: tree-aware estimator unavailable (" + treeFailure.getMessage()
                        + "); used the simpler pairwise-to-reference estimator instead");
            } catch (GviException e) {
                skipped.add("mu: " + e.getMessage() + (config.metadataPath() == null ? " (no --metadata supplied)" : ""));
            }
        }
    }

    /** {@link org.gvi.cli.GviCli}'s unmodified default for {@code --generation-time-days} -- used only to detect whether a user has actually overridden it. */
    private static final double DEFAULT_GENERATION_TIME_DAYS = 5.0;

    private void computeRe(PipelineConfig config, SequenceAlignment alignment, Map<IndexKey, IndexResult> population,
                            List<String> skipped, List<String> warnings) {
        if (!config.wants("re")) return;
        try {
            List<IncidencePoint> incidence = List.of();
            if (config.incidencePath() != null) {
                IncidenceCsvReader.Result r = IncidenceCsvReader.read(config.incidencePath());
                collectWarnings(warnings, "incidence", r.report());
                incidence = r.points();
            }

            // Resolve the generation time HERE rather than in a front end.
            //
            // It used to be resolved in GviCli before the config was built, so a pathogen whose
            // table entry is still a stub aborted the whole run with exit 2 -- no GVI at all, not
            // even the eight indices that need no generation time. That is the wrong disposal: a
            // missing generation time is absent data, and every other gate in this pipeline
            // excludes its own index and lets the rest score. MissingGenerationTimeException is a
            // GviException, so the catch below now turns it into "Re: <reason>" in the skipped
            // list, exactly like any other index that could not be computed.
            double generationTimeDays = org.gvi.algorithms.re.GenerationTimeTable.resolve(
                    config.generationTimeExplicit(), config.generationTimeDays(),
                    config.pathogenId(), config.generationTimeDays(),
                    org.gvi.algorithms.re.GenerationTimeTable.bundled());

            // Warn only when the generation time was neither typed explicitly nor resolved from the
            // per-pathogen table -- i.e. the run really is inheriting the historical 5-day default.
            boolean generationTimeUnsourced = !config.generationTimeExplicit()
                    && (config.pathogenId() == null || config.pathogenId().isBlank());
            if (incidence.isEmpty() && generationTimeUnsourced
                    && config.generationTimeDays() == DEFAULT_GENERATION_TIME_DAYS) {
                // Without real case counts, Re is estimated from a growth rate fitted to the tree and then
                // converted via the Euler-Lotka renewal equation, which needs a real serial-interval/generation-time
                // estimate for this specific pathogen -- 5 days is a fast-viral-epidemic default, not a safe
                // assumption. A short default mechanically compresses any growth-rate difference toward Re=1
                // (Re = 1 + r*generationTime/365 in the linear regime), so an unexamined default can make a
                // genuinely fast- or slow-spreading pathogen both read as "stable" for the same underlying reason.
                warnings.add("Re: no --incidence data and --generation-time-days was left at its default (5 days, "
                        + "tuned for a fast viral epidemic) -- Re is highly sensitive to this value and a wrong "
                        + "generation time can silently compress a real growth-rate difference toward Re=1, reading "
                        + "as false stability. Set --generation-time-days to a value appropriate for this pathogen's "
                        + "actual transmission cycle (e.g. much longer for a vector-borne parasite with a multi-week "
                        + "tick/vector cycle than for a directly-transmitted respiratory virus) before trusting this Re.");
            }

            // Native birth-death-sampling ML fit (BDSKY-equivalent), tried whenever no incidence series
            // is available (Cori et al. remains preferred whenever real case counts exist).
            //
            // This is now the DEFAULT rather than opt-in. The LTT/Euler-Lotka fallback below fits a
            // regression slope to a lineages-through-time curve, and on this project's corpus it
            // returned Re in 1.0007-1.0262 for every pathogen -- a spread too narrow to distinguish
            // anything, and largely insensitive to the generation time it is supposedly driven by.
            // BDSKY fits a real generative population model to the branching times instead: on the
            // same alignments it separates FMD (1.0917) from Bluetongue (1.0180), a spread roughly
            // eighty times wider. It is capped at 60 taxa and costs more, so the LTT fallback remains
            // for larger trees and for any dataset where the fit fails.
            if (incidence.isEmpty()) {
                try {
                    ReResult r = new org.gvi.algorithms.re.bdsky.BdskyReEstimator().compute(alignment, generationTimeDays);
                    population.put(IndexKey.RE, r);
                    return;
                } catch (GviException bdskyFailure) {
                    warnings.add("Re: the birth-death (BDSKY) fit was unavailable for this dataset ("
                            + bdskyFailure.getMessage() + "), so Re falls back to the lineages-through-time "
                            + "regression, which is markedly less able to distinguish growth rates -- treat the "
                            + "resulting Re as indicative only.");
                }
            }

            ReResult r = new ReCalculator().compute(incidence, alignment, generationTimeDays);
            population.put(IndexKey.RE, r);
        } catch (GviException e) {
            skipped.add("Re: " + e.getMessage());
        }
    }

    private void computeGd(PipelineConfig config, SequenceAlignment alignment, NucleotideSequence reference,
                            Map<String, Map<IndexKey, IndexResult>> perSequence, List<String> skipped) {
        if (!config.wants("gd")) return;
        GeneticDistanceCalculator calc = new GeneticDistanceCalculator();
        for (NucleotideSequence q : alignment.getQueries()) {
            try {
                GdResult r = calc.compute(reference, q, config.gdMethod());
                perSequence.get(q.getId()).put(IndexKey.GD, r);
            } catch (GviException e) {
                skipped.add("GD for '" + q.getId() + "': " + e.getMessage());
            }
        }
    }

    private void computeMb(PipelineConfig config, SequenceAlignment alignment, NucleotideSequence reference,
                            Map<String, Map<IndexKey, IndexResult>> perSequence, List<String> skipped) {
        if (!config.wants("mb")) return;
        MutationBurdenCalculator calc = new MutationBurdenCalculator();
        for (NucleotideSequence q : alignment.getQueries()) {
            try {
                MbResult r = calc.computeFromAlignment(reference, q);
                perSequence.get(q.getId()).put(IndexKey.MB, r);
            } catch (GviException e) {
                skipped.add("MB for '" + q.getId() + "': " + e.getMessage());
            }
        }
    }

    /**
     * Whether {@link #loadGenes} resolved gene coordinates from a real, user-supplied GFF3
     * ({@code autoPredicted=false}) or had to guess them with a 6-frame ORF scan
     * ({@code autoPredicted=true}) -- {@link #addDataQualityAdvisories} uses this to decide whether
     * dN/dS and CAI need an unconditional "this marker's coding status was never confirmed" caveat,
     * not just a caveat when the resulting value happens to look extreme.
     */
    private record GeneLoadResult(List<GeneAnnotation> genes, boolean autoPredicted) {
    }

    /** Persists the native ORF scan's predicted coordinates as a real GFF3 (--gff-out) -- a no-op unless genes were actually auto-predicted for this run. */
    private void writeGffOut(PipelineConfig config, GeneLoadResult geneLoad, NucleotideSequence reference, List<String> warnings) {
        if (config.gffOut() == null) return;
        if (!geneLoad.autoPredicted()) {
            warnings.add("--gff-out: no file written -- a real --gff was supplied and used successfully, so there's no auto-prediction to persist");
            return;
        }
        GffWriter.write(config.gffOut(), reference.getId(), geneLoad.genes());
        warnings.add("--gff-out: wrote " + geneLoad.genes().size() + " auto-predicted gene(s) to " + config.gffOut()
                + " -- inspect/correct before trusting it as a real annotation, or pass it back in via --gff on a later run");
    }

    private GeneLoadResult loadGenes(PipelineConfig config, SequenceAlignment alignment, List<String> warnings) {
        if (config.gffPath() != null) {
            GffReader.Result r;
            try {
                r = GffReader.read(config.gffPath());
            } catch (GviException e) {
                warnings.add("GFF3 (" + config.gffPath() + "): " + e.getMessage() + " -- falling back to native ORF prediction");
                return new GeneLoadResult(predictOrfs(alignment, warnings), true);
            }
            collectWarnings(warnings, "gff", r.report());
            if (!r.genes().isEmpty()) return new GeneLoadResult(r.genes(), false);
        }
        return new GeneLoadResult(predictOrfs(alignment, warnings), true);
    }

    /** Real ORFs in a compact viral genome are dramatically longer than the chance ORFs a 6-frame scan finds everywhere in a whole genome (an ~11kb sequence has dozens of >=30-codon runs by pure chance) -- keeping only ORFs within this fraction of the single longest one found is what actually separates real signal from that noise; confirmed empirically (real KFDV polyprotein: 3416 codons; next-longest, noise: 179 codons -- a 19x gap). */
    private static final double ORF_RELATIVE_LENGTH_THRESHOLD = 0.3;

    /**
     * No {@code --gff} supplied: rather than blindly treating the whole
     * sequence as one ORF, run a real (if simple) native ORF prediction --
     * {@link org.gvi.algorithms.orf.OrfFinder}, a standard 6-frame
     * ATG-to-stop scan -- on the reference sequence, using its coordinates
     * for every sequence in the alignment (consistent with how a supplied
     * GFF3's coordinates are already used dataset-wide). Falls back to the
     * single-whole-sequence-ORF treatment only if prediction finds nothing.
     * Filters down to ORFs of a comparable size to the longest one found
     * (see {@link #ORF_RELATIVE_LENGTH_THRESHOLD}) -- without this, dozens
     * of short chance ORFs from a whole-genome scan pollute dN/dS's
     * max-omega aggregation with noise (confirmed against real KFDV data:
     * an unfiltered scan drove the composite dN/dS to an implausible 45,
     * driven entirely by a handful of essentially meaningless codons in a
     * spurious 30-40 codon "ORF").
     */
    private List<GeneAnnotation> predictOrfs(SequenceAlignment alignment, List<String> warnings) {
        List<GeneAnnotation> allOrfs = new org.gvi.algorithms.orf.OrfFinder().findOrfs(alignment.getReference().getSequence());
        if (allOrfs.isEmpty()) {
            warnings.add("No --gff supplied and native ORF prediction found no open reading frame >= "
                    + org.gvi.algorithms.orf.OrfFinder.DEFAULT_MIN_ORF_CODONS
                    + " codons on the reference sequence; falling back to treating the whole sequence as a single ORF");
            return List.of(new GeneAnnotation("whole_genome", 1, alignment.length(), '+'));
        }

        long longest = allOrfs.stream().mapToLong(GeneAnnotation::length).max().orElseThrow();
        List<GeneAnnotation> predicted = allOrfs.stream()
                .filter(g -> g.length() >= longest * ORF_RELATIVE_LENGTH_THRESHOLD)
                .toList();

        warnings.add("No --gff supplied; predicted " + predicted.size() + " open reading frame(s) natively from the "
                + "reference sequence (ATG...stop 6-frame scan, >= " + org.gvi.algorithms.orf.OrfFinder.DEFAULT_MIN_ORF_CODONS
                + " codons, kept if within " + (int) (ORF_RELATIVE_LENGTH_THRESHOLD * 100) + "% of the longest one found -- "
                + allOrfs.size() + " candidate(s) before that filter) instead of treating the whole sequence as one ORF -- "
                + "a real but simpler technique than a trained gene-finder (no coding-potential model, no splice handling); "
                + "supply --gff for definitive gene coordinates if available.");
        return predicted;
    }

    /**
     * Computes CAI against every annotated gene (Section 5.8), not just the
     * first -- required to correctly handle genomes with multiple or
     * overlapping ORFs (common in compact viral genomes), which a
     * single-gene assumption would silently mis-annotate. The composite-facing
     * value is the mean CAI across genes (host translational adaptation is a
     * genome-wide property, not concentrated in one hotspot gene the way
     * selection pressure can be -- see computeDnDs for the contrasting
     * max-based aggregation there). Every gene's individual result is kept
     * in {@code caiPerGene} for full transparency in the report.
     */
    private void computeCai(PipelineConfig config, SequenceAlignment alignment, List<GeneAnnotation> genes,
                             Map<String, Map<IndexKey, IndexResult>> perSequence, Map<String, List<CaiResult>> caiPerGene,
                             List<String> skipped, List<String> warnings) {
        if (!config.wants("cai")) return;

        // A host-relative CAI presupposes the pathogen translates its proteins on HOST ribosomes,
        // using the host's tRNA pool -- true for viruses, false for bacteria and eukaryotic
        // parasites, which carry their own. Scoring those against a host table (or against E. coli
        // as a generic microbial stand-in) measures similarity to an unrelated third organism, not
        // host adaptation. Note the objection is to the reference SET, not to CAI: Sharp & Li (1987)
        // define CAI against highly expressed genes of the same organism, and scoring a bacterium
        // against its own ribosomal-protein genes would be perfectly valid.
        if (config.organismClassOrUnspecified().hasOwnTranslationMachinery()) {
            skipped.add("CAI: not computed for organism-class '"
                    + config.organismClassOrUnspecified().name().toLowerCase(java.util.Locale.ROOT)
                    + "'. Host-relative CAI measures adaptation to the host's tRNA pool, which only applies to organisms "
                    + "that translate on host ribosomes (viruses). This organism has its own ribosomes and tRNA genes, so "
                    + "a host codon table would measure similarity to an unrelated organism rather than host adaptation. "
                    + "A valid alternative is CAI against this organism's OWN highly-expressed genes (ribosomal proteins, "
                    + "elongation factors), which is how Sharp & Li originally defined it.");
            return;
        }

        if (config.organismClassOrUnspecified() == org.gvi.core.model.OrganismClass.UNSPECIFIED) {
            warnings.add("CAI: --organism-class was not supplied, so CAI is computed on the UNVERIFIED assumption that "
                    + "this organism translates its proteins on host ribosomes. That holds for viruses; for a bacterium "
                    + "or a eukaryotic parasite a host codon table measures similarity to an unrelated organism rather "
                    + "than host adaptation, and CAI would be excluded. Pass --organism-class to have this checked.");
        }

        CodonUsageTable table;
        if (config.codonUsagePath() != null) {
            table = CodonUsageTableReader.read(config.codonUsagePath(), config.codonUsagePath().getFileName().toString());
        } else if (config.codonUsageSpecies() != null && !config.codonUsageSpecies().isBlank()) {
            try {
                table = org.gvi.algorithms.cai.BundledCodonUsageTables.load(config.codonUsageSpecies());
                warnings.add("CAI: no --codon-usage supplied; using the bundled '" + config.codonUsageSpecies()
                        + "' reference table (real data from the Kazusa Codon Usage Database, not this dataset's own "
                        + "host) -- supply --codon-usage for a genuine dataset-specific host reference instead.");
            } catch (org.gvi.core.exception.GviInputException e) {
                // A misspelled species is the caller getting the flag wrong, not data being unavailable.
                // Skipping it here silently drops CAI from the composite and renormalizes the remaining
                // weights, so the run still prints a confident-looking GVI computed from fewer indices --
                // the failure mode is a plausible wrong number, not an obvious one. Abort instead.
                throw new org.gvi.core.exception.GviInputException(e.getMessage()
                        + " -- fix the --codon-usage-species value; the run is aborted rather than scored "
                        + "without CAI, because dropping an index silently changes the composite.");
            } catch (GviException e) {
                skipped.add("CAI: " + e.getMessage());
                return;
            }
        } else {
            skipped.add("CAI: no --codon-usage table or --codon-usage-species supplied");
            return;
        }
        CodonAdaptationCalculator calc = new CodonAdaptationCalculator();
        for (NucleotideSequence q : alignment.getQueries()) {
            List<CaiResult> perGene = new ArrayList<>();
            List<String> geneFailures = new ArrayList<>();
            for (GeneAnnotation gene : genes) {
                try {
                    String cds = org.gvi.algorithms.common.CodonUtil.extractCds(q.getSequence(), gene.start(), gene.end(), gene.strand());
                    perGene.add(calc.compute(q.getId(), gene.geneName(), cds, table));
                } catch (GviException e) {
                    geneFailures.add(gene.geneName() + ": " + e.getMessage());
                }
            }
            if (!geneFailures.isEmpty()) {
                skipped.add("CAI for '" + q.getId() + "', gene(s) skipped: " + String.join("; ", geneFailures));
            }
            if (perGene.isEmpty()) continue;
            caiPerGene.put(q.getId(), perGene);

            double meanCai = perGene.stream().mapToDouble(CaiResult::cai).average().orElse(0.0);
            var band = CaiReferenceTable.classify(meanCai);
            String category = band.status() + " (" + band.interpretation() + "; " + band.evolutionaryAge()
                    + ") -- mean across " + perGene.size() + " gene(s)";
            List<String> diagnostics = new ArrayList<>();
            for (CaiResult g : perGene) diagnostics.add(g.geneName() + ": CAI=" + String.format("%.4f", g.cai()));
            perSequence.get(q.getId()).put(IndexKey.CAI, new SimpleIndexResult("CAI", meanCai, category, diagnostics));
        }
    }

    /**
     * When {@code --reference-gc} isn't supplied, falls back to the
     * alignment's own reference sequence's GC% as the comparison baseline
     * -- derivable directly from the sequence data, unlike an external
     * host-genome target -- rather than skipping the index entirely. This
     * changes what "deviation" means (from THIS dataset's reference, not
     * from an assumed external target) and is flagged clearly in a warning
     * so it's never mistaken for the same thing.
     */
    private void computeGc(PipelineConfig config, SequenceAlignment alignment, NucleotideSequence reference,
                            Map<String, Map<IndexKey, IndexResult>> perSequence, List<String> skipped, List<String> warnings) {
        if (!config.wants("gc")) return;
        GcContentCalculator calc = new GcContentCalculator();
        Double referenceGcPercent = config.referenceGcPercent();
        if (referenceGcPercent == null) {
            try {
                referenceGcPercent = calc.compute(reference, 0.0).observedGcPercent();
                warnings.add(String.format("GC_Deviation: no --reference-gc supplied; using the alignment's own reference "
                                + "sequence's GC%% (%.2f%%, from '%s') as the baseline instead -- deviation is measured from "
                                + "THIS dataset's reference, not an external host-genome target; supply --reference-gc "
                                + "explicitly for that.",
                        referenceGcPercent, reference.getId()));
            } catch (GviException e) {
                skipped.add("GC_Deviation: no --reference-gc supplied and the reference sequence's own GC% could not be computed ("
                        + e.getMessage() + ")");
                return;
            }
        }
        boolean hgtLanguageApplies = config.organismClassOrUnspecified().hgtCompetent();
        boolean organismKnown = config.organismClassOrUnspecified() != org.gvi.core.model.OrganismClass.UNSPECIFIED;
        if (organismKnown && !hgtLanguageApplies) {
            warnings.add("GC_Deviation: organism-class is '"
                    + config.organismClassOrUnspecified().name().toLowerCase(java.util.Locale.ROOT)
                    + "', so a large GC deviation is NOT reported as horizontal gene transfer. Classical HGT "
                    + "(conjugation, transformation, transduction) is a documented mechanism in bacteria, but is not how "
                    + "viruses or eukaryotic parasites typically acquire compositional shifts -- for those, mutational "
                    + "pressure (e.g. APOBEC/ADAR editing) or host-driven composition bias is the usual explanation. "
                    + "The deviation value itself is unchanged; only its interpretation is.");
        }

        for (NucleotideSequence q : alignment.getQueries()) {
            try {
                GcResult r = calc.compute(q, referenceGcPercent, reference);
                perSequence.get(q.getId()).put(IndexKey.GC, retargetHgtWording(r, hgtLanguageApplies, organismKnown));
            } catch (GviException e) {
                skipped.add("GC_Deviation for '" + q.getId() + "': " + e.getMessage());
            }
        }
    }

    /**
     * Replaces the reference table's HGT wording when the organism cannot plausibly acquire
     * sequence by classical horizontal gene transfer. The numeric deviation and its band are
     * untouched -- only the causal claim attached to them changes, because "HGT suspected" on a
     * virus asserts a mechanism that does not apply.
     */
    private IndexResult retargetHgtWording(GcResult r, boolean hgtApplies, boolean organismKnown) {
        if (hgtApplies || !organismKnown || r.category() == null) return r;
        String rewritten = r.category()
                .replace("HGT suspected", "compositional shift (not HGT -- see organism-class note)")
                .replace("No HGT signal", "No compositional anomaly")
                .replace("Minor recombination / reassortment", "Minor compositional drift");
        if (rewritten.equals(r.category())) return r;
        List<String> diagnostics = new ArrayList<>(r.diagnostics());
        diagnostics.add("HGT wording suppressed for this organism class; the deviation value and band are unchanged.");
        return new SimpleIndexResult("GC_Deviation", r.primaryValue(), rewritten, diagnostics);
    }

    private static final int SLIDING_WINDOW_CODONS = 30;
    private static final int SLIDING_WINDOW_STEP_CODONS = 10;

    /**
     * Computes dN/dS against every annotated gene (Section 5.5), not just
     * the first (see computeCai for why this matters for overlapping/
     * multi-ORF genomes). The composite-facing value is the MAXIMUM omega
     * across genes: per Section 5.5's beta-integration language ("dN/dS > 1
     * in transmissibility-linked genes... directly increases predicted Re
     * and beta"), the single most concerning gene should drive the
     * composite signal rather than being diluted by an average across
     * genes under routine purifying selection. Every gene's full result
     * (in {@code dnDsPerGene}) and a site-resolved sliding-window scan (in
     * {@code dnDsSlidingWindows}) are kept for the detailed report.
     */
    private void computeDnDs(PipelineConfig config, SequenceAlignment alignment, NucleotideSequence reference,
                              List<GeneAnnotation> genes, Map<String, Map<IndexKey, IndexResult>> perSequence,
                              Map<String, List<DnDsResult>> dnDsPerGene, Map<String, List<SlidingWindowDnDsResult>> dnDsSlidingWindows,
                              Map<String, DnDsResult> winningDnDsPerSequence, List<String> skipped) {
        if (!config.wants("dnds")) return;
        DnDsCalculator calc = new DnDsCalculator();
        for (NucleotideSequence q : alignment.getQueries()) {
            List<DnDsResult> perGene = new ArrayList<>();
            List<SlidingWindowDnDsResult> windows = new ArrayList<>();
            List<String> geneFailures = new ArrayList<>();
            for (GeneAnnotation gene : genes) {
                try {
                    String refCds = org.gvi.algorithms.common.CodonUtil.extractCds(reference.getSequence(), gene.start(), gene.end(), gene.strand());
                    String qryCds = org.gvi.algorithms.common.CodonUtil.extractCds(q.getSequence(), gene.start(), gene.end(), gene.strand());
                    perGene.add(calc.compute(q.getId(), gene.geneName(), refCds, qryCds));
                    windows.addAll(calc.slidingWindow(q.getId(), gene.geneName(), refCds, qryCds, SLIDING_WINDOW_CODONS, SLIDING_WINDOW_STEP_CODONS));
                } catch (GviException e) {
                    geneFailures.add(gene.geneName() + ": " + e.getMessage());
                }
            }
            if (!geneFailures.isEmpty()) {
                skipped.add("dN/dS for '" + q.getId() + "', gene(s) skipped: " + String.join("; ", geneFailures));
            }
            if (perGene.isEmpty()) continue;
            dnDsPerGene.put(q.getId(), perGene);
            if (!windows.isEmpty()) dnDsSlidingWindows.put(q.getId(), windows);

            DnDsResult worst = perGene.stream().max(Comparator.comparingDouble(DnDsResult::omega)).orElseThrow();
            winningDnDsPerSequence.put(q.getId(), worst);
            List<String> diagnostics = new ArrayList<>();
            diagnostics.add("Composite value is the MAXIMUM omega across " + perGene.size() + " gene(s), from '" + worst.geneName() + "'");
            for (DnDsResult g : perGene) diagnostics.add(g.geneName() + ": omega=" + String.format("%.4f", g.omega()));
            perSequence.get(q.getId()).put(IndexKey.DNDS, new SimpleIndexResult("dN/dS", worst.omega(), worst.category(), diagnostics));
        }
    }

    /**
     * Opt-in ({@code --ml-dnds}): a SEPARATE, dataset-wide (not per-sequence) maximum-likelihood dN/dS
     * per gene via {@link org.gvi.algorithms.dnds.ml.MlCodonDnDsEstimator} (GY94/{@code codeml}
     * M0-equivalent, joint-tree ML fit) -- additional cross-checking detail alongside the default
     * Nei-Gojobori dN/dS, not a replacement for it (see {@link PipelineResult}'s javadoc). Tried per
     * gene independently so one oversized/undersized gene doesn't block the others; skipped genes get
     * a clear reason recorded rather than silently vanishing.
     */
    private Map<String, org.gvi.algorithms.dnds.ml.MlDnDsResult> computeMlDnDs(PipelineConfig config, SequenceAlignment alignment,
                                                                                 List<GeneAnnotation> genes, List<String> warnings, List<String> skipped) {
        Map<String, org.gvi.algorithms.dnds.ml.MlDnDsResult> results = new LinkedHashMap<>();
        if (!config.wants("dnds") || !config.mlDnds()) return results;

        org.gvi.algorithms.dnds.ml.MlCodonDnDsEstimator estimator = new org.gvi.algorithms.dnds.ml.MlCodonDnDsEstimator();
        for (GeneAnnotation gene : genes) {
            try {
                org.gvi.algorithms.dnds.ml.MlDnDsResult r = estimator.compute(alignment, gene);
                results.put(gene.geneName(), r);
                warnings.add(String.format("ML dN/dS (--ml-dnds) for gene '%s': omega=%.4f kappa=%.4f (GY94/codeml "
                                + "M0-equivalent, %d taxa, %d codons, %d cycles) -- cross-check only, not used in the composite GVI",
                        gene.geneName(), r.omega(), r.kappa(), r.taxonCount(), r.codonsUsed(), r.cyclesUsed()));
            } catch (GviException e) {
                skipped.add("ML dN/dS (--ml-dnds) for gene '" + gene.geneName() + "': " + e.getMessage());
            }
        }
        return results;
    }

    /**
     * Opt-in ({@code --weights}): loads a user-supplied composite-weights JSON instead of the spec's
     * default midpoints -- the read/merge/renormalize/unknown-label handling all lives in
     * {@link CompositeWeightsReader}; this just decides whether to call it and reports the outcome.
     */
    private CompositeWeights loadWeights(PipelineConfig config, List<String> warnings) {
        if (config.weightsPath() == null) return CompositeWeights.defaults();
        try {
            CompositeWeightsReader.Result r = CompositeWeightsReader.read(config.weightsPath());
            for (String w : r.warnings()) warnings.add("Custom weights (" + config.weightsPath() + "): " + w);
            warnings.add("Using custom composite weights from " + config.weightsPath()
                    + " (any index not mentioned there kept its spec-default weight; the full set was renormalized to sum to 1.0)");
            return r.weights();
        } catch (GviException e) {
            warnings.add("Custom weights (" + config.weightsPath() + ") could not be loaded (" + e.getMessage()
                    + "); using spec defaults instead");
            return CompositeWeights.defaults();
        }
    }

    private void collectWarnings(List<String> warnings, String source, ParseReport report) {
        for (ParseWarning w : report.getWarnings()) {
            warnings.add(source + ": " + w);
        }
    }
}
