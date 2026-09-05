package org.gvi.algorithms.re;

import org.gvi.core.model.IncidencePoint;
import org.junit.jupiter.api.Test;

import java.time.LocalDate;
import java.util.ArrayList;
import java.util.List;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.within;

class CoriReEstimatorTest {

    @Test
    void recoversTheTrueReproductionNumberFromARenewalEquationSimulatedSeries() {
        // Simulate incidence deterministically from the SAME renewal equation the
        // estimator inverts: I(t) = R_true * sum_u I(t-u) w(u). If the estimator's
        // math is correct, its posterior mean should converge tightly to R_true.
        double rTrue = 1.5;
        SerialInterval si = new SerialInterval(5.0, 2.0, 20);
        LocalDate start = LocalDate.of(2024, 1, 1);

        double[] incidence = new double[80];
        for (int i = 0; i < si.tMax(); i++) incidence[i] = 100.0; // seed period
        for (int t = si.tMax(); t < incidence.length; t++) {
            double lambda = 0;
            for (int u = 1; u <= si.tMax(); u++) {
                lambda += incidence[t - u] * si.weight(u);
            }
            incidence[t] = rTrue * lambda;
        }

        List<IncidencePoint> points = new ArrayList<>();
        for (int i = 0; i < incidence.length; i++) {
            points.add(new IncidencePoint(start.plusDays(i), incidence[i]));
        }

        CoriReEstimator estimator = new CoriReEstimator();
        ReResult result = estimator.compute(points, si);

        assertThat(result.reMean()).isCloseTo(rTrue, within(0.02));
        assertThat(result.method()).isEqualTo(ReMethod.CORI_INCIDENCE);
        assertThat(result.interval()).isNotNull();
        assertThat(result.interval().lower()).isLessThan(result.reMean());
        assertThat(result.interval().upper()).isGreaterThan(result.reMean());
    }

    @Test
    void classifiesDecliningEpidemicAsControlled() {
        double rTrue = 0.6;
        SerialInterval si = new SerialInterval(5.0, 2.0, 20);
        LocalDate start = LocalDate.of(2024, 1, 1);

        double[] incidence = new double[80];
        for (int i = 0; i < si.tMax(); i++) incidence[i] = 1000.0;
        for (int t = si.tMax(); t < incidence.length; t++) {
            double lambda = 0;
            for (int u = 1; u <= si.tMax(); u++) lambda += incidence[t - u] * si.weight(u);
            incidence[t] = rTrue * lambda;
        }
        List<IncidencePoint> points = new ArrayList<>();
        for (int i = 0; i < incidence.length; i++) points.add(new IncidencePoint(start.plusDays(i), incidence[i]));

        ReResult result = new CoriReEstimator().compute(points, si);
        assertThat(result.reMean()).isLessThan(0.8);
        assertThat(result.category()).contains("Controlled");
    }
}
