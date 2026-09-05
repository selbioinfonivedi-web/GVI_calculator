package org.gvi.algorithms.cai;

import java.util.List;

/**
 * Whether the C->T and A->G substitutions between a reference and query
 * genome are statistically enriched at the known APOBEC3/ADAR dinucleotide
 * hotspot context, versus the same substitution elsewhere in the genome --
 * see {@link MutationalSignatureAnalyzer}. Gives GC_Deviation's compositional
 * shift a mechanistic explanation (host RNA-editing pressure) as an
 * alternative to the index's default "HGT suspected" interpretation.
 */
public record MutationalSignatureResult(SignatureTest apobecLikeC2T, SignatureTest adarLikeA2G, List<String> diagnostics) {

    public record SignatureTest(String label, int hotspotContextMutated, int hotspotContextTotal,
                                 int otherContextMutated, int otherContextTotal,
                                 double oddsRatio, double pValue, boolean significant) {

        static SignatureTest inconclusive(String label, int hotspotMutated, int hotspotTotal, int otherMutated, int otherTotal) {
            return new SignatureTest(label, hotspotMutated, hotspotTotal, otherMutated, otherTotal, Double.NaN, Double.NaN, false);
        }
    }

    public boolean anySignatureDetected() {
        return apobecLikeC2T.significant() || adarLikeA2G.significant();
    }
}
