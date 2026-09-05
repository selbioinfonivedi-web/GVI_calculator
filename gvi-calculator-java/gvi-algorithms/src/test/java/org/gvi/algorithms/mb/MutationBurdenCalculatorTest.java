package org.gvi.algorithms.mb;

import org.gvi.core.model.NucleotideSequence;
import org.gvi.core.model.VariantRecord;
import org.junit.jupiter.api.Test;

import java.util.List;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.within;

class MutationBurdenCalculatorTest {

    private final MutationBurdenCalculator calc = new MutationBurdenCalculator(10);

    @Test
    void vcfPathFiltersLowCoverageVariants() {
        List<VariantRecord> variants = List.of(
                new VariantRecord("chr1", 100, "A", "G", 25, "s1"), // passes
                new VariantRecord("chr1", 200, "C", "T", 5, "s1"),  // below threshold
                new VariantRecord("chr1", 300, "AT", "A", null, "s1") // unknown depth, included
        );
        MbResult r = calc.computeFromVariants("s1", variants);
        assertThat(r.mutationBurden()).isEqualTo(2);
        assertThat(r.excludedLowCoverage()).isEqualTo(1);
        assertThat(r.excludedAmbiguousOrUnknownCoverage()).isEqualTo(1);
        assertThat(r.source()).isEqualTo(MbResult.MbSource.VCF_WITH_COVERAGE_FILTER);
    }

    @Test
    void alignmentPathCountsSnpsIndividuallyAndCollapsesIndelRuns() {
        // ref: ACGT---ACGT   query: ACGTGGGACGT  -> insertion of GGG = 1 indel event, no SNPs
        NucleotideSequence ref = new NucleotideSequence("ref", "ACGT---ACGT");
        NucleotideSequence qry = new NucleotideSequence("qry", "ACGTGGGACGT");
        MbResult r = calc.computeFromAlignment(ref, qry);
        assertThat(r.mutationBurden()).isEqualTo(1);
        assertThat(r.source()).isEqualTo(MbResult.MbSource.ALIGNMENT_NO_COVERAGE_FILTER);
    }

    @Test
    void alignmentPathCountsMultipleSnpsSeparately() {
        NucleotideSequence ref = new NucleotideSequence("ref", "ACGTACGT");
        NucleotideSequence qry = new NucleotideSequence("qry", "AGGTACGA"); // positions 1 and 7 differ
        MbResult r = calc.computeFromAlignment(ref, qry);
        assertThat(r.mutationBurden()).isEqualTo(2);
    }

    @Test
    void alignmentPathCollapsesDeletionRunIntoOneEvent() {
        // ref has extra bases that query lacks -> deletion run of 3
        NucleotideSequence ref = new NucleotideSequence("ref", "ACGTGGGACGT");
        NucleotideSequence qry = new NucleotideSequence("qry", "ACGT---ACGT");
        MbResult r = calc.computeFromAlignment(ref, qry);
        assertThat(r.mutationBurden()).isEqualTo(1);
    }

    @Test
    void identicalSequencesHaveZeroBurden() {
        NucleotideSequence ref = new NucleotideSequence("ref", "ACGTACGT");
        NucleotideSequence qry = new NucleotideSequence("qry", "ACGTACGT");
        MbResult r = calc.computeFromAlignment(ref, qry);
        assertThat(r.mutationBurden()).isZero();
        assertThat(r.category()).contains("Within-household");
    }

    @Test
    void outbreakAgeWeeksIsDimensionallyCorrect() {
        // mu = 1e-3 sub/site/yr, genome = 30000 nt -> 30 substitutions/year genome-wide
        // MB = 15 -> years = 0.5 -> weeks = 0.5 * 52.1775 = 26.08875
        double weeks = MutationBurdenCalculator.estimateOutbreakAgeWeeks(15, 1e-3, 30000);
        assertThat(weeks).isCloseTo(26.08875, within(1e-6));
    }
}
