package org.gvi.algorithms.phylo.model;

import org.junit.jupiter.api.Test;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.within;

class NamedSubstitutionModelsTest {

    private static final double[] EQUAL_FREQ = {0.25, 0.25, 0.25, 0.25};

    /**
     * Every named model must reduce to JC69's known closed form when its
     * extra parameters are set to their "no effect" values (kappa=1 means
     * transitions and transversions are equally likely, i.e. no
     * transition/transversion bias; equal frequencies) -- this is the
     * defining nested-model relationship, and a strong correctness check
     * since GtrModel's JC69 behavior is itself independently verified
     * (GtrModelTest) against the textbook closed form.
     */
    @Test
    void k80WithKappaOneReducesToJukesCantor() {
        assertReducesToJc69(NamedSubstitutionModels.k80().buildModel(new double[]{1.0}, EQUAL_FREQ));
    }

    @Test
    void hky85WithKappaOneAndEqualFreqReducesToJukesCantor() {
        assertReducesToJc69(NamedSubstitutionModels.hky85().buildModel(new double[]{1.0}, EQUAL_FREQ));
    }

    @Test
    void tn93WithBothKappasOneAndEqualFreqReducesToJukesCantor() {
        assertReducesToJc69(NamedSubstitutionModels.tn93().buildModel(new double[]{1.0, 1.0}, EQUAL_FREQ));
    }

    @Test
    void gtrWithAllRatesOneAndEqualFreqReducesToJukesCantor() {
        assertReducesToJc69(NamedSubstitutionModels.gtr().buildModel(new double[]{1, 1, 1, 1, 1}, EQUAL_FREQ));
    }

    @Test
    void f81WithEqualFrequenciesReducesToJukesCantor() {
        assertReducesToJc69(NamedSubstitutionModels.f81().buildModel(new double[0], EQUAL_FREQ));
    }

    private void assertReducesToJc69(GtrModel model) {
        GtrModel jc = GtrModel.jukesCantor();
        for (double t : new double[]{0.05, 0.3, 1.0}) {
            double[][] actual = model.transitionProbabilities(t);
            double[][] expected = jc.transitionProbabilities(t);
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    assertThat(actual[i][j]).as("P[%d][%d] at t=%s", i, j, t).isCloseTo(expected[i][j], within(1e-9));
                }
            }
        }
    }

    @Test
    void parameterCountsMatchStandardConvention() {
        assertThat(NamedSubstitutionModels.jc69().totalSubstitutionParameterCount()).isEqualTo(0);
        assertThat(NamedSubstitutionModels.f81().totalSubstitutionParameterCount()).isEqualTo(3);
        assertThat(NamedSubstitutionModels.k80().totalSubstitutionParameterCount()).isEqualTo(1);
        assertThat(NamedSubstitutionModels.hky85().totalSubstitutionParameterCount()).isEqualTo(4);
        assertThat(NamedSubstitutionModels.tn93().totalSubstitutionParameterCount()).isEqualTo(5);
        assertThat(NamedSubstitutionModels.gtr().totalSubstitutionParameterCount()).isEqualTo(8);
    }

    @Test
    void byNameIsCaseInsensitiveAndRejectsUnknownModels() {
        assertThat(NamedSubstitutionModels.byName("gtr").modelName()).isEqualTo("GTR");
        assertThat(NamedSubstitutionModels.byName("HkY85").modelName()).isEqualTo("HKY85");
        org.assertj.core.api.Assertions.assertThatThrownBy(() -> NamedSubstitutionModels.byName("nonsense"))
                .isInstanceOf(IllegalArgumentException.class);
    }
}
