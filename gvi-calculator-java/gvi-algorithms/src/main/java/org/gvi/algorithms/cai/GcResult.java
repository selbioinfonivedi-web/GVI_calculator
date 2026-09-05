package org.gvi.algorithms.cai;

import org.gvi.core.spi.IndexResult;

import java.util.List;

/**
 * Result of the GC Content Deviation half of Index 8 (Section 5.8).
 * {@code mutationalSignature} is null unless the actual reference sequence
 * was supplied (see {@link GcContentCalculator#compute(NucleotideSequence,
 * double, NucleotideSequence)}) -- it tests whether the compositional shift
 * is explained by known host RNA-editing pressure (APOBEC3/ADAR) rather
 * than the index's default "HGT suspected" reading.
 */
public record GcResult(String sequenceId, double observedGcPercent, double referenceGcPercent,
                        double deviationPercent, String category, List<String> diagnostics,
                        MutationalSignatureResult mutationalSignature) implements IndexResult {

    @Override
    public String indexName() {
        return "GC_Deviation";
    }

    @Override
    public double primaryValue() {
        return deviationPercent;
    }
}
