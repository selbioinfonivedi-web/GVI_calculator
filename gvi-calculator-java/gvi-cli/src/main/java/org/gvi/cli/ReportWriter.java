package org.gvi.cli;

import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.SerializationFeature;
import com.fasterxml.jackson.datatype.jsr310.JavaTimeModule;
import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVPrinter;
import org.gvi.composite.GviComponent;
import org.gvi.composite.GviResult;
import org.gvi.composite.IndexKey;
import org.gvi.composite.SensitivityResult;
import org.gvi.core.spi.IndexResult;

import java.io.IOException;
import java.io.PrintStream;
import java.io.Writer;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Comparator;
import java.util.EnumMap;
import java.util.List;
import java.util.Map;

/** Formats a {@link PipelineResult} as human-readable text, JSON, or a per-sequence CSV summary. */
public final class ReportWriter {

    private final ObjectMapper mapper;

    public ReportWriter() {
        this.mapper = new ObjectMapper()
                .registerModule(new JavaTimeModule())
                .enable(SerializationFeature.INDENT_OUTPUT)
                .disable(SerializationFeature.WRITE_DATES_AS_TIMESTAMPS);
    }

    /** Convenience overload: whole-file result only, no per-sequence section -- see {@link #writeText(PipelineResult, PrintStream, boolean)}. */
    public void writeText(PipelineResult result, PrintStream out) {
        writeText(result, out, false);
    }

    /**
     * @param includePerSequence whether to also print the per-sequence detail section. The headline
     *                           result is always the whole-file one regardless of this flag -- this only
     *                           controls whether the (potentially very long, one block per sequence)
     *                           supporting-evidence section is included too.
     */
    public void writeText(PipelineResult result, PrintStream out, boolean includePerSequence) {
        out.println("=== Genomic Virulence Index (GVI) Report ===");
        out.println();
        out.println("############################################");
        out.println("##  WHOLE-FILE RESULT (all sequences combined)");
        out.println("############################################");
        if (result.datasetGvi() != null) {
            out.printf("%nGVI = %.4f   [%s]%n", result.datasetGvi().gvi(), result.datasetGvi().coverageSummary());
            // A low-coverage GVI is numerically indistinguishable from a well-supported one, so the caveat has
            // to sit next to the number itself -- not only in the warnings block far below it.
            if (!result.datasetGvi().comparable()) {
                out.println("  *** NOT COMPARABLE across datasets -- too little of the weighting scheme survived quality");
                out.println("      gating. See the warning below before ranking this value against anything. ***");
            }
            out.println();
            for (GviComponent c : result.datasetGvi().components()) {
                out.printf("  %-14s raw=%-12.6g norm=%.4f weight=%.4f contribution=%.4f%n",
                        c.key().label(), c.rawValue(), c.normalizedValue(), c.effectiveWeight(), c.contribution());
            }
            if (!result.datasetGvi().excludedIndices().isEmpty()) {
                out.println("  excluded: " + result.datasetGvi().excludedIndices());
            }
        } else {
            out.println("\n(Whole-file GVI could not be computed -- see warnings below.)");
        }
        out.println();
        out.println("-- Component index detail (dataset-wide) --");
        for (var entry : result.populationIndices().entrySet()) {
            printIndex(out, "  ", entry.getKey(), entry.getValue());
        }
        for (var entry : result.datasetIndices().entrySet()) {
            printIndex(out, "  ", entry.getKey(), entry.getValue());
        }

        if (!result.sensitivity().isEmpty()) {
            out.println();
            out.println("-- Sensitivity analysis (each index's weight perturbed +/-20%, most sensitive first) --");
            List<SensitivityResult> sorted = result.sensitivity().stream()
                    .sorted(Comparator.comparingDouble(SensitivityResult::spread).reversed())
                    .toList();
            for (SensitivityResult s : sorted) {
                out.printf("  %-14s GVI range [%.4f, %.4f]  (spread=%.4f, base=%.4f)%n",
                        s.key().label(), s.gviAtLowWeight(), s.gviAtHighWeight(), s.spread(), s.baseGvi());
            }
        }

        if (includePerSequence) {
            out.println();
            out.println("-- Per-sequence detail (supporting evidence, not the headline result) --");
            for (var seqEntry : result.perSequenceIndices().entrySet()) {
                String id = seqEntry.getKey();
                out.println();
                out.println("[" + id + "]");
                for (var entry : seqEntry.getValue().entrySet()) {
                    printIndex(out, "    ", entry.getKey(), entry.getValue());
                }
                GviResult gvi = result.gviPerSequence().get(id);
                if (gvi != null) {
                    out.printf("    GVI (this sequence only) = %.4f%n", gvi.gvi());
                }
            }
        } else if (!result.perSequenceIndices().isEmpty()) {
            out.println();
            out.printf("-- Per-sequence detail omitted (%d sequence(s) -- pass --per-sequence to include it) --%n",
                    result.perSequenceIndices().size());
        }

        if (!result.mlDnDsPerGene().isEmpty()) {
            out.println();
            out.println("-- ML dN/dS (--ml-dnds, GY94/codeml M0-equivalent, cross-check only -- not used in the composite GVI) --");
            for (var entry : result.mlDnDsPerGene().entrySet()) {
                var r = entry.getValue();
                out.printf("  %-20s omega=%.4f kappa=%.4f logL=%.2f (%d taxa, %d codons, %d cycles)%n",
                        entry.getKey(), r.omega(), r.kappa(), r.logLikelihood(), r.taxonCount(), r.codonsUsed(), r.cyclesUsed());
            }
        }

        if (!result.reassortment().isEmpty()) {
            out.println();
            out.println("-- Reassortment (segment concordance, --segment; separate from the composite GVI) --");
            for (var r : result.reassortment()) {
                out.printf("  %-10s vs %-10s  mantel_r=%-8.4f p=%-8.4f reassortment_index=%.4f (%d taxa)%n",
                        r.segmentALabel(), r.segmentBLabel(), r.mantelR(), r.pValue(), r.reassortmentIndex(), r.taxaCompared());
                out.println("      " + r.category());
            }
        }

        if (!result.skipped().isEmpty()) {
            out.println();
            out.println("-- Skipped --");
            result.skipped().forEach(s -> out.println("  - " + s));
        }
        if (!result.warnings().isEmpty()) {
            out.println();
            out.println("-- Warnings --");
            result.warnings().forEach(w -> out.println("  - " + w));
        }
    }

    private void printIndex(PrintStream out, String indent, IndexKey key, IndexResult r) {
        out.printf("%s%-14s value=%-12.6g category=%s%n", indent, key.label(), r.primaryValue(), r.category());
        for (String d : r.diagnostics()) {
            out.println(indent + "    * " + d);
        }
    }

    /** Convenience overload: whole-file result only, no per-sequence-keyed fields -- see {@link #writeJson(PipelineResult, Path, boolean)}. */
    public void writeJson(PipelineResult result, Path path) {
        writeJson(result, path, false);
    }

    /**
     * @param includePerSequence whether to also include the per-sequence-keyed fields
     *                           ({@code per_sequence_indices}, {@code gvi_per_sequence},
     *                           {@code cai_per_gene}, {@code dnds_per_gene}, {@code dnds_sliding_windows}
     *                           -- each keyed by sequence id). {@code ml_dnds_per_gene} is dataset-wide
     *                           (keyed by gene, not sequence) and is always included regardless.
     */
    public void writeJson(PipelineResult result, Path path, boolean includePerSequence) {
        try {
            mapper.writeValue(path.toFile(), toJsonModel(result, includePerSequence));
        } catch (IOException e) {
            throw new org.gvi.core.exception.GviInputException("Could not write JSON report to " + path + ": " + e.getMessage(), e);
        }
    }

    private Map<String, Object> toJsonModel(PipelineResult result, boolean includePerSequence) {
        Map<String, Object> root = new java.util.LinkedHashMap<>();
        root.put("dataset_gvi", result.datasetGvi());
        root.put("dataset_indices", result.datasetIndices());
        root.put("population_indices", result.populationIndices());
        if (includePerSequence) {
            root.put("per_sequence_indices", result.perSequenceIndices());
            root.put("gvi_per_sequence", result.gviPerSequence());
            root.put("cai_per_gene", result.caiPerGene());
            root.put("dnds_per_gene", result.dnDsPerGene());
            root.put("dnds_sliding_windows", result.dnDsSlidingWindows());
        }
        root.put("ml_dnds_per_gene", result.mlDnDsPerGene());
        root.put("reassortment", result.reassortment());
        root.put("sensitivity", result.sensitivity());
        root.put("skipped", result.skipped());
        root.put("warnings", result.warnings());
        return root;
    }

    /** One row: the whole-file result. Per-sequence detail is in the JSON report; the CSV is the single-answer format. */
    public void writeCsv(PipelineResult result, Path path) {
        try (Writer w = Files.newBufferedWriter(path, StandardCharsets.UTF_8);
             CSVPrinter printer = new CSVPrinter(w, CSVFormat.DEFAULT)) {
            // comparable/coverage_pct/indices_scored sit immediately after GVI on purpose. A CSV is what
            // gets pasted into a spreadsheet and sorted, and a GVI built from two indices is renormalized
            // into [0,1] exactly like one built from nine -- so without these columns the file invites a
            // ranking it cannot support.
            printer.printRecord("GVI", "comparable", "coverage_pct", "indices_scored",
                    "mu", "Re", "pi", "MB", "dN/dS", "GD", "CAI", "GC_Deviation", "RI", "sequences_in_file");
            Map<IndexKey, IndexResult> combined = new EnumMap<>(result.datasetIndices());
            combined.putAll(result.populationIndices());
            GviResult dataset = result.datasetGvi();
            printer.printRecord(
                    dataset != null ? dataset.gvi() : "",
                    dataset == null ? "" : (dataset.comparable() ? "yes" : "NO"),
                    dataset == null ? "" : String.format(java.util.Locale.ROOT, "%.1f", dataset.effectiveWeightSum() * 100),
                    dataset == null ? "" : dataset.components().size(),
                    valueOrBlank(combined, IndexKey.MU),
                    valueOrBlank(combined, IndexKey.RE),
                    valueOrBlank(combined, IndexKey.PI),
                    valueOrBlank(combined, IndexKey.MB),
                    valueOrBlank(combined, IndexKey.DNDS),
                    valueOrBlank(combined, IndexKey.GD),
                    valueOrBlank(combined, IndexKey.CAI),
                    valueOrBlank(combined, IndexKey.GC),
                    valueOrBlank(combined, IndexKey.RI),
                    result.perSequenceIndices().size()
            );
        } catch (IOException e) {
            throw new org.gvi.core.exception.GviInputException("Could not write CSV report to " + path + ": " + e.getMessage(), e);
        }
    }

    /** Optional: the previous per-sequence CSV format, for anyone who wants the per-sample breakdown as a spreadsheet. */
    public void writePerSequenceCsv(PipelineResult result, Path path) {
        try (Writer w = Files.newBufferedWriter(path, StandardCharsets.UTF_8);
             CSVPrinter printer = new CSVPrinter(w, CSVFormat.DEFAULT)) {
            printer.printRecord("sequence_id", "GVI", "mu", "Re", "pi", "MB", "dN/dS", "GD", "CAI", "GC_Deviation", "RI");
            for (var entry : result.perSequenceIndices().entrySet()) {
                String id = entry.getKey();
                Map<IndexKey, IndexResult> combined = new EnumMap<>(entry.getValue());
                combined.putAll(result.populationIndices());
                GviResult gvi = result.gviPerSequence().get(id);
                printer.printRecord(
                        id,
                        gvi != null ? gvi.gvi() : "",
                        valueOrBlank(combined, IndexKey.MU),
                        valueOrBlank(combined, IndexKey.RE),
                        valueOrBlank(combined, IndexKey.PI),
                        valueOrBlank(combined, IndexKey.MB),
                        valueOrBlank(combined, IndexKey.DNDS),
                        valueOrBlank(combined, IndexKey.GD),
                        valueOrBlank(combined, IndexKey.CAI),
                        valueOrBlank(combined, IndexKey.GC),
                        valueOrBlank(combined, IndexKey.RI)
                );
            }
        } catch (IOException e) {
            throw new org.gvi.core.exception.GviInputException("Could not write CSV report to " + path + ": " + e.getMessage(), e);
        }
    }

    private Object valueOrBlank(Map<IndexKey, IndexResult> map, IndexKey key) {
        IndexResult r = map.get(key);
        return r == null ? "" : r.primaryValue();
    }
}
