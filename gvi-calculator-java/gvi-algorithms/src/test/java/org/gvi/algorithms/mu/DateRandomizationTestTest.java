package org.gvi.algorithms.mu;

import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import static org.assertj.core.api.Assertions.assertThat;

/**
 * The date-randomization test exists because R² cannot distinguish a real molecular clock from
 * lineage structure. These tests pin exactly that distinction.
 */
class DateRandomizationTestTest {

    private static List<DateRandomizationTest.TemporalPoint> points(double[][] raw) {
        List<DateRandomizationTest.TemporalPoint> out = new ArrayList<>();
        for (double[] r : raw) out.add(new DateRandomizationTest.TemporalPoint(r[0], r[1]));
        return out;
    }

    @Test
    void genuineClockSignalPasses() {
        // Distance accumulates steadily with date: a real clock at 0.01 sub/site/yr.
        double[][] raw = new double[12][];
        for (int i = 0; i < 12; i++) {
            double year = 2010 + i;
            raw[i] = new double[]{year, (year - 2010) * 0.01};
        }
        var r = new DateRandomizationTest().run(points(raw), 0.01);
        assertThat(r.passed()).isTrue();
        assertThat(r.explanation()).contains("passed");
    }

    @Test
    void lineageStructureWithNoRealTemporalSignalFails() {
        // Two tight clusters of distance, dates assigned independently of which cluster a point is
        // in. A regression through this yields a confident-looking slope that has nothing to do with
        // elapsed time -- the exact failure mode seen on the real Sheep & Goat Pox alignment.
        Random rng = new Random(7);
        double[][] raw = new double[16][];
        for (int i = 0; i < 16; i++) {
            double year = 2010 + rng.nextInt(10);
            double distance = (i % 2 == 0) ? 0.001 : 0.045; // lineage A vs lineage B
            raw[i] = new double[]{year, distance};
        }
        // Whatever slope the real assignment gives, shuffled dates reproduce the same kind of value.
        var pts = points(raw);
        double realSlope = slopeOf(pts);
        var r = new DateRandomizationTest().run(pts, realSlope);
        assertThat(r.passed()).isFalse();
        assertThat(r.explanation()).contains("FAILED").contains("independently introduced");
    }

    @Test
    void identicalDatesCannotSupportAClock() {
        double[][] raw = {{2020, 0.01}, {2020, 0.02}, {2020, 0.03}, {2020, 0.04}};
        var r = new DateRandomizationTest().run(points(raw), 0.01);
        assertThat(r.passed()).isFalse();
        assertThat(r.explanation()).contains("same collection date");
    }

    @Test
    void tooFewPointsIsReportedAsUnverifiedNotConfirmed() {
        var r = new DateRandomizationTest().run(points(new double[][]{{2020, 0.01}, {2021, 0.02}}), 0.01);
        assertThat(r.explanation()).contains("unverified");
        assertThat(r.replicatesRun()).isZero();
    }

    @Test
    void isDeterministic() {
        double[][] raw = new double[10][];
        for (int i = 0; i < 10; i++) raw[i] = new double[]{2010 + i, i * 0.01};
        var a = new DateRandomizationTest().run(points(raw), 0.01);
        var b = new DateRandomizationTest().run(points(raw), 0.01);
        assertThat(a.randomizedMin()).isEqualTo(b.randomizedMin());
        assertThat(a.randomizedMax()).isEqualTo(b.randomizedMax());
    }

    private static double slopeOf(List<DateRandomizationTest.TemporalPoint> pts) {
        double mx = pts.stream().mapToDouble(DateRandomizationTest.TemporalPoint::decimalYear).average().orElse(0);
        double my = pts.stream().mapToDouble(DateRandomizationTest.TemporalPoint::distanceToRoot).average().orElse(0);
        double num = 0, den = 0;
        for (var p : pts) {
            double dx = p.decimalYear() - mx;
            num += dx * (p.distanceToRoot() - my);
            den += dx * dx;
        }
        return den == 0 ? 0 : num / den;
    }
}
