package org.gvi.algorithms.cai;

import org.gvi.core.model.NucleotideSequence;
import org.junit.jupiter.api.Test;

import static org.assertj.core.api.Assertions.assertThat;

class MutationalSignatureAnalyzerTest {

    private final MutationalSignatureAnalyzer analyzer = new MutationalSignatureAnalyzer();

    /**
     * 15 C's preceded by T (hotspot context) + 15 C's preceded by G (other
     * context); mutate 12/15 hotspot C's to T and 0/15 other-context C's --
     * a textbook APOBEC3 TCW-style pattern that a real signature scan must
     * flag as significant.
     */
    @Test
    void detectsApobecLikeSignatureWhenCToTIsEnrichedAtHotspotContext() {
        StringBuilder refB = new StringBuilder();
        StringBuilder qryB = new StringBuilder();
        for (int i = 0; i < 15; i++) {
            refB.append("TC");
            qryB.append(i < 12 ? "TT" : "TC"); // 12 of 15 hotspot C's mutate to T
        }
        for (int i = 0; i < 15; i++) {
            refB.append("GC");
            qryB.append("GC"); // 0 of 15 other-context C's mutate
        }
        NucleotideSequence reference = new NucleotideSequence("ref", refB.toString());
        NucleotideSequence query = new NucleotideSequence("qry", qryB.toString());

        MutationalSignatureResult result = analyzer.analyze(reference, query);

        assertThat(result.apobecLikeC2T().significant()).isTrue();
        assertThat(result.apobecLikeC2T().pValue()).isLessThan(0.05);
        assertThat(result.apobecLikeC2T().oddsRatio()).isGreaterThan(1.0);
        assertThat(result.apobecLikeC2T().hotspotContextMutated()).isEqualTo(12);
        assertThat(result.apobecLikeC2T().otherContextMutated()).isEqualTo(0);
        assertThat(result.anySignatureDetected()).isTrue();

        // no 'A' bases anywhere in this fixture -- the ADAR test has nothing to work with and must say so, not fabricate a verdict
        assertThat(result.adarLikeA2G().pValue()).isNaN();
        assertThat(result.adarLikeA2G().significant()).isFalse();
    }

    @Test
    void doesNotFlagSignatureWhenMutationRateIsUniformAcrossContexts() {
        StringBuilder refB = new StringBuilder();
        StringBuilder qryB = new StringBuilder();
        for (int i = 0; i < 15; i++) {
            refB.append("TC");
            qryB.append(i < 2 ? "TT" : "TC"); // 2 of 15 mutate, same rate as below
        }
        for (int i = 0; i < 15; i++) {
            refB.append("GC");
            qryB.append(i < 2 ? "GT" : "GC"); // 2 of 15 mutate -- same rate, no context enrichment
        }
        NucleotideSequence reference = new NucleotideSequence("ref", refB.toString());
        NucleotideSequence query = new NucleotideSequence("qry", qryB.toString());

        MutationalSignatureResult result = analyzer.analyze(reference, query);

        assertThat(result.apobecLikeC2T().significant()).isFalse();
        assertThat(result.anySignatureDetected()).isFalse();
    }

    @Test
    void reportsInconclusiveWhenTooFewSitesInEitherContext() {
        // only 3 hotspot-context C's -- well under MIN_SITES_PER_CONTEXT
        NucleotideSequence reference = new NucleotideSequence("ref", "TCTCTC" + "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGC");
        NucleotideSequence query = new NucleotideSequence("qry", "TTTTTT" + "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGC");

        MutationalSignatureResult result = analyzer.analyze(reference, query);

        assertThat(result.apobecLikeC2T().pValue()).isNaN();
        assertThat(result.apobecLikeC2T().significant()).isFalse();
    }
}
