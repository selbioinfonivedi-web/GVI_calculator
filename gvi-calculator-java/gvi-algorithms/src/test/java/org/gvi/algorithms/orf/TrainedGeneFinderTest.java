package org.gvi.algorithms.orf;

import org.gvi.core.model.GeneAnnotation;
import org.junit.jupiter.api.Test;

import java.util.List;

import static org.assertj.core.api.Assertions.assertThat;

/**
 * The whole point of a trained gene finder over a length-only heuristic: it must (1) keep a real
 * gene even when that gene is far shorter than the longest ORF in the genome, and (2) reject a
 * spurious ORF of comparable length whose composition doesn't match the genome's own trained codon
 * usage. Neither property is exercised by testing length filtering alone -- this is the same kind of
 * "does the new thing actually do more than the old thing" check the project's own history already
 * needed once for the length-only heuristic (its javadoc cites a real KFDV case where an unfiltered
 * scan drove a composite index to an implausible value, entirely from noise ORFs the length filter
 * let through).
 */
class TrainedGeneFinderTest {

    // 12 sense codons, deliberately disjoint from NOISE_CODONS, cycled to build every "real" gene's body.
    private static final String[] REAL_CODONS = {
            "GCC", "GAA", "CTG", "AAG", "TTC", "CAG", "ATC", "GGC", "ACC", "CCG", "TGG", "CGC"
    };
    // 12 different sense codons, never used by a real gene, cycled to build the one spurious ORF.
    private static final String[] NOISE_CODONS = {
            "AAT", "TTT", "CCA", "GGA", "ACA", "TCA", "AGA", "GTA", "CAT", "TAT", "GAT", "CTA"
    };

    private record Placed(String genome, long start, long end) {
    }

    private static String gene(String[] codons, int bodyCodons) {
        StringBuilder sb = new StringBuilder("ATG");
        for (int i = 0; i < bodyCodons; i++) sb.append(codons[i % codons.length]);
        sb.append("TAA");
        return sb.toString();
    }

    @Test
    void keepsAShortRealGeneAndRejectsALongerSpuriousOneOfMismatchedComposition() {
        // Frame-0-safe filler: "AAT" never spells ATG at a codon boundary and (being one of the
        // NOISE_CODONS, not a stop) never closes a segment either -- it exists only to keep the
        // genome's overall base composition from being dominated by either codon set, and an
        // explicit "TAA" always separates it from whatever comes next.
        String filler = "AAT".repeat(60) + "TAA";

        StringBuilder genome = new StringBuilder(filler);
        long[] trainingCoords = new long[6]; // start1, end1, start2, end2, start3, end3
        for (int i = 0; i < 3; i++) {
            String g = gene(REAL_CODONS, 100); // 102 codons -- clears MIN_TRAINING_ORF_CODONS (90)
            trainingCoords[i * 2] = genome.length() + 1L;
            genome.append(g);
            trainingCoords[i * 2 + 1] = genome.length();
            genome.append(filler);
        }

        String shortRealGene = gene(REAL_CODONS, 35); // 37 codons -- clears OrfFinder's own 30-codon candidate floor, but well under the 90-codon training floor
        long shortRealStart = genome.length() + 1L;
        genome.append(shortRealGene);
        long shortRealEnd = genome.length();
        genome.append(filler);

        String spuriousOrf = gene(NOISE_CODONS, 40); // 42 codons -- longer than the real short gene, wrong composition
        long spuriousStart = genome.length() + 1L;
        genome.append(spuriousOrf);
        long spuriousEnd = genome.length();
        genome.append(filler);

        String sequence = genome.toString();
        List<TrainedGeneFinder.ScoredGene> scored = new TrainedGeneFinder().scoreAll(sequence);

        double shortRealScore = scoreOf(scored, shortRealStart, shortRealEnd);
        double spuriousScore = scoreOf(scored, spuriousStart, spuriousEnd);
        System.out.println("[TrainedGeneFinderTest] short real gene (37 codons) log-odds=" + shortRealScore
                + "; spurious ORF (42 codons, wrong composition) log-odds=" + spuriousScore);

        assertThat(shortRealScore).as("a real gene's own codon usage must score positive against the trained model")
                .isGreaterThan(0.0);
        assertThat(spuriousScore).as("length alone (42 > 22 codons) must not save a composition that never appears in training")
                .isLessThan(0.0);

        List<GeneAnnotation> kept = new TrainedGeneFinder().findGenes(sequence);
        assertThat(kept).as("the short real gene must be kept despite being far shorter than the training genes")
                .anyMatch(g -> g.start() == shortRealStart && g.end() == shortRealEnd);
        assertThat(kept).as("the longer spurious ORF must be rejected despite being longer than the real gene it's compared against")
                .noneMatch(g -> g.start() == spuriousStart && g.end() == spuriousEnd);
        for (int i = 0; i < 3; i++) {
            long s = trainingCoords[i * 2], e = trainingCoords[i * 2 + 1];
            assertThat(kept).as("training gene %d must itself be kept", i + 1).anyMatch(g -> g.start() == s && g.end() == e);
        }
    }

    private double scoreOf(List<TrainedGeneFinder.ScoredGene> scored, long start, long end) {
        return scored.stream()
                .filter(s -> s.gene().start() == start && s.gene().end() == end && s.gene().strand() == '+')
                .mapToDouble(TrainedGeneFinder.ScoredGene::codingLogOdds)
                .findFirst()
                .orElseThrow(() -> new AssertionError("no candidate found at [" + start + "," + end + "]"));
    }
}
