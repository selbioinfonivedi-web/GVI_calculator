package org.gvi.algorithms.cai;

import org.gvi.algorithms.common.CodonUtil;
import org.gvi.core.exception.GviComputationException;
import org.gvi.core.model.CodonUsageTable;
import org.gvi.core.util.GeneticCode;

import java.util.ArrayList;
import java.util.List;

/**
 * CAI half of Index 8 (Section 5.8), Sharp &amp; Li 1987:
 * CAI = exp( (1/L) * sum(ln w(codon_i)) ), the geometric mean of per-codon
 * relative adaptiveness weights, computed in log-space for numerical
 * stability. Stop codons and single-codon amino acids (Met, Trp) are
 * excluded from L per convention (their w is always 1 / undefined and
 * carries no discriminative information).
 */
public final class CodonAdaptationCalculator {

    public CaiResult compute(String sequenceId, String geneName, String cdsNucleotides, CodonUsageTable table) {
        CodonUtil.SplitResult split = CodonUtil.splitIntoCodons(cdsNucleotides);
        List<String> codons = split.codons();

        double sumLogW = 0.0;
        int used = 0, excludedStop = 0, excludedSingleAa = 0, excludedInvalid = 0;

        for (String codon : codons) {
            if (!GeneticCode.isStandardCodon(codon)) {
                excludedInvalid++;
                continue;
            }
            char aa = GeneticCode.translate(codon);
            if (aa == GeneticCode.STOP) {
                excludedStop++;
                continue;
            }
            if (GeneticCode.SYNONYMOUS_CODONS.get(aa).size() == 1) {
                excludedSingleAa++;
                continue;
            }
            double w = table.weightOf(codon);
            sumLogW += Math.log(w);
            used++;
        }

        if (used == 0) {
            throw new GviComputationException("CAI for '" + sequenceId + "'/'" + geneName
                    + "': no codons usable after excluding stops/invalid/single-codon amino acids");
        }

        double cai = Math.exp(sumLogW / used);
        var band = CaiReferenceTable.classify(cai);
        String category = band.status() + " (" + band.interpretation() + "; " + band.evolutionaryAge() + ")";

        List<String> diagnostics = new ArrayList<>();
        if (split.droppedTrailingBases() > 0) {
            diagnostics.add(split.droppedTrailingBases() + " trailing base(s) dropped (CDS length not a multiple of 3)");
        }
        if (excludedInvalid > 0) {
            diagnostics.add(excludedInvalid + " codon(s) contained gaps/ambiguity codes and were excluded");
        }
        diagnostics.add(used + " codons used (excluded " + excludedStop + " stop, " + excludedSingleAa + " single-codon-AA, " + excludedInvalid + " invalid)");

        return new CaiResult(sequenceId, geneName, cai, used, excludedStop, excludedSingleAa, excludedInvalid, category, diagnostics);
    }
}
