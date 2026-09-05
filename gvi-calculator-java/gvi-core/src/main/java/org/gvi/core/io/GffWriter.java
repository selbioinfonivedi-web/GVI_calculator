package org.gvi.core.io;

import org.gvi.core.exception.GviInputException;
import org.gvi.core.model.GeneAnnotation;

import java.io.IOException;
import java.io.Writer;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

/**
 * Writes {@link GeneAnnotation}s as a real GFF3 file -- the counterpart to {@link GffReader}, and what
 * {@code --gff-out} uses to persist the coordinates {@code GviPipeline}'s native ORF scan predicted when
 * no {@code --gff} was supplied. Previously that prediction only ever existed in memory for one run;
 * writing it out lets a user inspect it, hand-correct it, or feed it back in via {@code --gff} on a
 * later run instead of re-guessing every time.
 */
public final class GffWriter {

    private static final String SOURCE = "gvi-calculator-orf-scan";

    private GffWriter() {
    }

    public static void write(Path path, String seqId, List<GeneAnnotation> genes) {
        try (Writer w = Files.newBufferedWriter(path, StandardCharsets.UTF_8)) {
            w.write("##gff-version 3\n");
            w.write("# Auto-predicted by gvi-calculator's native 6-frame ORF scan (no --gff was supplied for this "
                    + "run) -- not a curated annotation. Coding status was never confirmed; treat gene boundaries "
                    + "here as a starting point, not ground truth.\n");
            for (GeneAnnotation gene : genes) {
                w.write(String.join("\t",
                        seqId, SOURCE, "CDS",
                        Long.toString(gene.start()), Long.toString(gene.end()),
                        ".", String.valueOf(gene.strand()), "0",
                        "ID=" + gene.geneName() + ";Name=" + gene.geneName()));
                w.write("\n");
            }
        } catch (IOException e) {
            throw new GviInputException("Could not write GFF3 to " + path + ": " + e.getMessage(), e);
        }
    }
}
