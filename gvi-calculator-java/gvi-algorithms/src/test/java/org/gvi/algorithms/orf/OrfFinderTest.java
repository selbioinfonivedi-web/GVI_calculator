package org.gvi.algorithms.orf;

import org.gvi.core.model.GeneAnnotation;
import org.junit.jupiter.api.Test;

import java.util.List;

import static org.assertj.core.api.Assertions.assertThat;

class OrfFinderTest {

    @Test
    void findsAKnownForwardOrfAtTheExpectedCoordinates() {
        String prefix = "A".repeat(10);           // junk: codon AAA (Lysine), no start/stop signal in any frame
        String orf = "ATG" + "GCC".repeat(32) + "TAA"; // start + 32 codons + stop = 34 codons, >= default 30-codon minimum
        String suffix = "A".repeat(10);
        String sequence = prefix + orf + suffix;

        List<GeneAnnotation> found = new OrfFinder().findOrfs(sequence);

        List<GeneAnnotation> forwardOrfs = found.stream().filter(g -> g.strand() == '+').toList();
        assertThat(forwardOrfs).hasSize(1);
        GeneAnnotation result = forwardOrfs.get(0);
        assertThat(result.start()).isEqualTo(prefix.length() + 1L); // 1-based
        assertThat(result.end()).isEqualTo(prefix.length() + orf.length());
        assertThat(result.length()).isEqualTo(orf.length());
    }

    @Test
    void rejectsAnOrfShorterThanTheMinimumLength() {
        // start + 5 codons + stop = 7 codons, well under the 30-codon default minimum
        String shortOrf = "ATG" + "GCC".repeat(5) + "TAA";
        String sequence = "A".repeat(10) + shortOrf + "A".repeat(10);

        List<GeneAnnotation> found = new OrfFinder().findOrfs(sequence);

        assertThat(found).noneMatch(g -> g.strand() == '+' && g.start() == 11);
    }

    @Test
    void aLowerMinimumLengthRecoversTheSameShortOrf() {
        String shortOrf = "ATG" + "GCC".repeat(5) + "TAA"; // 7 codons
        String sequence = "A".repeat(10) + shortOrf + "A".repeat(10);

        List<GeneAnnotation> found = new OrfFinder(5).findOrfs(sequence);

        assertThat(found).anyMatch(g -> g.strand() == '+' && g.start() == 11 && g.length() == shortOrf.length());
    }

    @Test
    void takesTheFirstAtgInAStopToStopSegmentNotALaterOne() {
        // Two in-frame ATGs before the stop -- the ORF should start at the FIRST one (the longer possible ORF),
        // the standard "longest ORF per stop-to-stop segment" convention.
        String orf = "ATG" + "ATG" + "GCC".repeat(31) + "TAA"; // first ATG, then a second in-frame ATG, then 31 more codons + stop
        String sequence = "A".repeat(9) + orf + "A".repeat(9); // 9 (not 10) to keep the whole thing in frame 0 cleanly

        List<GeneAnnotation> found = new OrfFinder().findOrfs(sequence);

        List<GeneAnnotation> forwardOrfs = found.stream().filter(g -> g.strand() == '+').toList();
        assertThat(forwardOrfs).hasSize(1);
        assertThat(forwardOrfs.get(0).start()).isEqualTo(10L); // the FIRST Atg's position, not the second
    }
}
