package org.gvi.algorithms.phylo.model;

import org.gvi.core.model.NucleotideSequence;
import org.gvi.core.model.SequenceAlignment;
import org.junit.jupiter.api.Test;

import java.util.List;
import java.util.Set;

import static org.assertj.core.api.Assertions.assertThat;

class SubstitutionModelSelectorTest {

    private final SubstitutionModelSelector selector = new SubstitutionModelSelector();

    /**
     * IMPORTANT design point (found by an earlier, wrong version of this
     * test failing): base composition bias and transition/transversion
     * RATE bias are different things and must not be conflated. A dataset
     * built from only A and G bases (with C/T never appearing at all)
     * doesn't actually test kappa -- F81's frequency-weighted mechanism
     * alone already explains "A always becomes G" purely because G is
     * common and C/T are ~never observed, with no need for a rate-ratio
     * parameter. AIC correctly preferred F81 on that data; that was the
     * selector working right and the test's data being wrong.
     * <p>
     * This version keeps composition balanced (reference cycles
     * A,C,G,T equally) and introduces ONLY transitions (A<->G, C<->T) at
     * roughly symmetric positions, so a real ts/tv rate signal exists that
     * composition alone can't explain -- isolating the thing this test is
     * actually meant to check.
     */
    @Test
    void prefersATransitionTransversionAwareModelWhenTheDataHasStrongTsTvBias() {
        int repeats = 100;
        StringBuilder refBuilder = new StringBuilder();
        for (int i = 0; i < repeats; i++) refBuilder.append("ACGT");
        char[] base = refBuilder.toString().toCharArray(); // 400bp, exactly 25% each base
        NucleotideSequence ref = new NucleotideSequence("ref", new String(base));

        NucleotideSequence q1 = transitionMutate(base, "q1", 0, 15);
        NucleotideSequence q2 = transitionMutate(base, "q2", 15, 35);

        SequenceAlignment aln = SequenceAlignment.of(List.of(ref, q1, q2), "ref");

        SubstitutionModelSelector.SelectionResult result = selector.select(aln, NamedSubstitutionModels.all(), 0);

        assertThat(result.allFits()).hasSize(6);
        // fits must be sorted best (lowest AIC) first
        for (int i = 1; i < result.allFits().size(); i++) {
            assertThat(result.allFits().get(i).aic()).isGreaterThanOrEqualTo(result.allFits().get(i - 1).aic());
        }
        assertThat(result.best().modelName()).isEqualTo(result.allFits().get(0).modelName());

        Set<String> tsTvAwareModels = Set.of("K80", "HKY85", "TN93", "GTR");
        assertThat(tsTvAwareModels).as("best model (%s) should account for the strong ts/tv bias", result.best().modelName())
                .contains(result.best().modelName());
    }

    /** Mutates positions [fromIdx, toIdx) each to its OWN transition partner (A<->G, C<->T) -- never a transversion. */
    private static NucleotideSequence transitionMutate(char[] base, String id, int fromIdx, int toIdx) {
        char[] copy = base.clone();
        for (int i = fromIdx; i < toIdx; i++) {
            copy[i] = switch (base[i]) {
                case 'A' -> 'G';
                case 'G' -> 'A';
                case 'C' -> 'T';
                case 'T' -> 'C';
                default -> base[i];
            };
        }
        return new NucleotideSequence(id, new String(copy));
    }
}
