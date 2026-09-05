package org.gvi.selftest;

import org.gvi.algorithms.cai.CodonAdaptationCalculator;
import org.gvi.algorithms.cai.GcContentCalculator;
import org.gvi.algorithms.dnds.DnDsCalculator;
import org.gvi.algorithms.gd.GdMethod;
import org.gvi.algorithms.gd.GdResult;
import org.gvi.algorithms.gd.GeneticDistanceCalculator;
import org.gvi.algorithms.mb.MbResult;
import org.gvi.algorithms.mb.MutationBurdenCalculator;
import org.gvi.algorithms.mu.EvolutionaryRateCalculator;
import org.gvi.algorithms.mu.MuResult;
import org.gvi.algorithms.phylo.NeighborJoining;
import org.gvi.algorithms.phylo.PhyloTree;
import org.gvi.algorithms.pi.NucleotideDiversityCalculator;
import org.gvi.algorithms.pi.PiResult;
import org.gvi.algorithms.re.CoriReEstimator;
import org.gvi.algorithms.re.ReResult;
import org.gvi.algorithms.re.SerialInterval;
import org.gvi.algorithms.ri.RecombinationIndexCalculator;
import org.gvi.algorithms.ri.RiResult;
import org.gvi.composite.CompositeGviEngine;
import org.gvi.composite.CompositeWeights;
import org.gvi.composite.GviResult;
import org.gvi.composite.IndexKey;
import org.gvi.composite.WeightCalibrator;
import org.gvi.core.io.FastaReader;
import org.gvi.core.model.CodonUsageTable;
import org.gvi.core.model.IncidencePoint;
import org.gvi.core.model.NucleotideSequence;
import org.gvi.core.model.SequenceAlignment;

import java.time.LocalDate;
import java.util.ArrayList;
import java.util.EnumMap;
import java.util.List;
import java.util.Map;

/**
 * Bundled offline self-diagnostic suite (Section 8.3 of the build spec):
 * each check runs one algorithm module against a small, hand-verifiable
 * dataset built into this class (no external files, no network -- works
 * air-gapped) and compares the result to a value that was worked out by
 * hand (see comments) against a numerical tolerance. Run via
 * {@code --self-test} on first launch or on demand, so a user can confirm
 * this specific installation is computing correctly before trusting it
 * with real data.
 */
public final class SelfTestRunner {

    private static final double TOL = 1e-6;

    public SelfTestReport runAll() {
        List<SelfTestCase> cases = new ArrayList<>();
        cases.add(run("Genetic Distance (Jukes-Cantor)", this::checkGd));
        cases.add(run("Nucleotide Diversity (pi)", this::checkPi));
        cases.add(run("Mutation Burden", this::checkMb));
        cases.add(run("Evolutionary Rate (mu)", this::checkMu));
        cases.add(run("dN/dS (Nei-Gojobori)", this::checkDnDs));
        cases.add(run("Codon Adaptation Index", this::checkCai));
        cases.add(run("GC Content Deviation", this::checkGc));
        cases.add(run("Effective Reproduction Number (Re)", this::checkRe));
        cases.add(run("Recombination Index (PHI test)", this::checkRi));
        cases.add(run("Composite GVI", this::checkComposite));
        cases.add(run("Neighbor-Joining tree construction", this::checkNeighborJoining));
        cases.add(run("Weight calibration", this::checkWeightCalibration));
        cases.add(run("Crash resilience (malformed FASTA input)", this::checkCrashResilience));
        return SelfTestReport.of(cases);
    }

    private interface Check {
        void run();
    }

    private SelfTestCase run(String name, Check check) {
        try {
            check.run();
            return new SelfTestCase(name, true, null);
        } catch (AssertionError e) {
            return new SelfTestCase(name, false, e.getMessage());
        } catch (Exception e) {
            return new SelfTestCase(name, false, "Unexpected exception: " + e);
        }
    }

    private void assertClose(String label, double expected, double actual, double tol) {
        if (Math.abs(expected - actual) > tol) {
            throw new AssertionError(label + ": expected " + expected + " but got " + actual + " (tolerance " + tol + ")");
        }
    }

    // p=0.1 over 20 sites -> d_JC = -0.75*ln(1-4/3*0.1)
    private void checkGd() {
        String ref = "A".repeat(20);
        StringBuilder q = new StringBuilder(ref);
        q.setCharAt(0, 'G');
        q.setCharAt(1, 'G');
        GdResult r = new GeneticDistanceCalculator().compute(
                new NucleotideSequence("ref", ref), new NucleotideSequence("q", q.toString()), GdMethod.JUKES_CANTOR);
        double expected = -0.75 * Math.log(1 - (4.0 / 3.0) * 0.1);
        assertClose("Jukes-Cantor distance", expected, r.distance(), TOL);
    }

    // 3 sequences, hand-computed pi = 0.4/3 (see gvi-algorithms NucleotideDiversityCalculatorTest)
    private void checkPi() {
        SequenceAlignment aln = SequenceAlignment.of(List.of(
                new NucleotideSequence("s1", "AAAAAAAAAA"),
                new NucleotideSequence("s2", "AAAAAAAAAG"),
                new NucleotideSequence("s3", "AAAAAAAAGG")
        ));
        PiResult r = new NucleotideDiversityCalculator(1L, 20, 100_000L).compute(aln);
        assertClose("nucleotide diversity", 0.4 / 3.0, r.pi(), TOL);
    }

    // insertion of 3bp collapses to exactly 1 mutation event
    private void checkMb() {
        MbResult r = new MutationBurdenCalculator(10).computeFromAlignment(
                new NucleotideSequence("ref", "ACGT---ACGT"), new NucleotideSequence("qry", "ACGTGGGACGT"));
        if (r.mutationBurden() != 1) {
            throw new AssertionError("mutation burden: expected 1, got " + r.mutationBurden());
        }
    }

    // exact 2-point regression slope: (d2-d1)/(t2-t1)
    private void checkMu() {
        String base = "A".repeat(1000);
        NucleotideSequence ref = new NucleotideSequence("ref", base).withMetadata(LocalDate.of(2020, 1, 1), null, null);
        StringBuilder mutant = new StringBuilder(base);
        for (int i = 0; i < 10; i++) mutant.setCharAt(i, 'G'); // p = 0.01
        NucleotideSequence q = new NucleotideSequence("q", mutant.toString()).withMetadata(LocalDate.of(2021, 1, 1), null, null);
        SequenceAlignment aln = SequenceAlignment.of(List.of(ref, q), "ref");
        MuResult r = new EvolutionaryRateCalculator().compute(aln);
        double expectedD = -0.75 * Math.log(1 - (4.0 / 3.0) * 0.01);
        double expectedMu = expectedD / 1.0; // ~1 year apart
        assertClose("mu (2-point slope)", expectedMu, r.muSubPerSiteYear(), 1e-3);
    }

    // TTT(Phe)->TTA(Leu), single-codon gene: hand-computed S=0.5, N=2.5, Sd=0, Nd=1
    // -> pS=0 (dS=0, triggers the dS==0 fallback where omega=dN), pN=0.4 -> dN=-0.75*ln(1-4/3*0.4)=0.571605...
    private void checkDnDs() {
        // TTT(Phe)->TTA(Leu), single-codon gene: hand-computed S=0.5, N=2.5, Sd=0, Nd=1.
        // Zero observed synonymous differences with S>0 triggers the continuity correction
        // (pS = 0.5/(S+1)) documented in DnDsCalculator.omegaFrom, rather than a literal dS=0.
        var r = new DnDsCalculator().compute("selftest", "geneX", "TTT", "TTA");
        double expectedPS = 0.5 / (0.5 + 1);
        double expectedDS = -0.75 * Math.log(1 - (4.0 / 3.0) * expectedPS);
        double expectedPN = 1.0 / 2.5;
        double expectedDN = -0.75 * Math.log(1 - (4.0 / 3.0) * expectedPN);
        assertClose("dN/dS: dS", expectedDS, r.dS(), 1e-9);
        assertClose("dN/dS: dN", expectedDN, r.dN(), 1e-9);
        assertClose("dN/dS: omega", expectedDN / expectedDS, r.omega(), 1e-9);
    }

    // all-optimal-codon gene must have CAI exactly 1.0
    private void checkCai() {
        CodonUsageTable table = CodonUsageTable.fromFrequencies("selftest", Map.of("TTT", 1.0, "TTC", 0.0));
        var r = new CodonAdaptationCalculator().compute("s1", "geneX", "TTT".repeat(10), table);
        assertClose("CAI", 1.0, r.cai(), TOL);
    }

    // 4 G/C of 10 bases -> 40% GC
    private void checkGc() {
        var r = new GcContentCalculator().compute(new NucleotideSequence("s1", "GCGCAAAAAA"), 38.0);
        assertClose("GC%", 40.0, r.observedGcPercent(), TOL);
    }

    // deterministic renewal-equation simulation: estimator must recover R_true=1.4 to within 0.05
    private void checkRe() {
        double rTrue = 1.4;
        SerialInterval si = new SerialInterval(5.0, 2.0, 20);
        LocalDate start = LocalDate.of(2024, 1, 1);
        double[] incidence = new double[80];
        for (int i = 0; i < si.tMax(); i++) incidence[i] = 100.0;
        for (int t = si.tMax(); t < incidence.length; t++) {
            double lambda = 0;
            for (int u = 1; u <= si.tMax(); u++) lambda += incidence[t - u] * si.weight(u);
            incidence[t] = rTrue * lambda;
        }
        List<IncidencePoint> points = new ArrayList<>();
        for (int i = 0; i < incidence.length; i++) points.add(new IncidencePoint(start.plusDays(i), incidence[i]));
        ReResult r = new CoriReEstimator().compute(points, si);
        assertClose("Re", rTrue, r.reMean(), 0.05);
    }

    // cross-cutting bipartitions (recombination signature) must be flagged significant
    private void checkRi() {
        int[] region1 = {0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95};
        int[] region2 = {100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 185, 190, 195};
        List<NucleotideSequence> seqs = new ArrayList<>();
        for (int i = 0; i < 24; i++) {
            char[] bases = new char[200];
            java.util.Arrays.fill(bases, 'A');
            char a1 = i < 12 ? 'A' : 'G';
            char a2 = i % 2 == 0 ? 'A' : 'G';
            for (int p : region1) bases[p] = a1;
            for (int p : region2) bases[p] = a2;
            seqs.add(new NucleotideSequence("ind" + i, new String(bases)));
        }
        RiResult r = new RecombinationIndexCalculator(100, 200, 42L).computeFromAlignment(SequenceAlignment.of(seqs));
        if (!r.significant()) {
            throw new AssertionError("RI: expected a statistically significant recombination signal, got p=" + r.pValue());
        }
    }

    // two equal-weight components at raw=5 in range [0,10] -> normalized 0.5 each -> GVI=0.5
    private void checkComposite() {
        Map<IndexKey, Double> weightMap = new EnumMap<>(IndexKey.class);
        weightMap.put(IndexKey.MU, 0.5);
        weightMap.put(IndexKey.RE, 0.5);
        CompositeWeights weights = CompositeWeights.of(weightMap);
        Map<IndexKey, org.gvi.composite.NormalizationRange> ranges = new EnumMap<>(IndexKey.class);
        ranges.put(IndexKey.MU, new org.gvi.composite.NormalizationRange(0, 10));
        ranges.put(IndexKey.RE, new org.gvi.composite.NormalizationRange(0, 10));
        Map<IndexKey, org.gvi.core.spi.IndexResult> available = new EnumMap<>(IndexKey.class);
        available.put(IndexKey.MU, fake("mu", 5.0));
        available.put(IndexKey.RE, fake("Re", 5.0));
        GviResult r = new CompositeGviEngine().compute(available, weights, ranges);
        assertClose("composite GVI", 0.5, r.gvi(), TOL);
    }

    private org.gvi.core.spi.IndexResult fake(String name, double value) {
        return new org.gvi.core.spi.SimpleIndexResult(name, value, "selftest", List.of());
    }

    // Neighbor-Joining is provably exact on an additive distance matrix (see NeighborJoiningTest for the derivation)
    private void checkNeighborJoining() {
        String[] labels = {"A", "B", "C", "D"};
        double[][] distances = {
                {0, 3, 5, 8},
                {3, 0, 6, 9},
                {5, 6, 0, 5},
                {8, 9, 5, 0}
        };
        NeighborJoining.Result result = NeighborJoining.build(labels, distances);
        PhyloTree tree = PhyloTree.rootAt(result, "A");
        assertClose("NJ root-to-tip(B)", 3.0, tree.rootToTip("B"), TOL);
        assertClose("NJ root-to-tip(C)", 5.0, tree.rootToTip("C"), TOL);
        assertClose("NJ root-to-tip(D)", 8.0, tree.rootToTip("D"), TOL);
    }

    // noise-free synthetic data with known weights [0.5, 0.5] must be recovered by the optimizer
    private void checkWeightCalibration() {
        List<IndexKey> keys = List.of(IndexKey.MU, IndexKey.RE);
        List<WeightCalibrator.Observation> observations = new ArrayList<>();
        double[][] samples = {{0.2, 0.8}, {0.9, 0.1}, {0.4, 0.6}, {0.7, 0.3}, {0.5, 0.5}, {0.1, 0.9}};
        for (double[] s : samples) {
            Map<IndexKey, Double> values = new EnumMap<>(IndexKey.class);
            values.put(IndexKey.MU, s[0]);
            values.put(IndexKey.RE, s[1]);
            observations.add(new WeightCalibrator.Observation(values, 0.5 * s[0] + 0.5 * s[1]));
        }
        WeightCalibrator.CalibrationResult result = new WeightCalibrator().calibrate(observations, keys);
        assertClose("calibrated weight(mu)", 0.5, result.weights().get(IndexKey.MU), 0.02);
        assertClose("calibrated weight(Re)", 0.5, result.weights().get(IndexKey.RE), 0.02);
    }

    // deliberately malformed input must be rejected cleanly, never crash the JVM
    private void checkCrashResilience() {
        String malformed = ">seq1\nACGT123XYZ!!\n>seq2 duplicate-of-nothing\n\n>seq2\nACGT\n";
        FastaReader.Result result = FastaReader.read(new java.io.StringReader(malformed), "selftest-malformed");
        if (result.report().getSkippedCount() == 0) {
            throw new AssertionError("Expected malformed records to be flagged as skipped, none were");
        }
    }
}
