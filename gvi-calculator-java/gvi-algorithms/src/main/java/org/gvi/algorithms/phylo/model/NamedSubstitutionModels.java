package org.gvi.algorithms.phylo.model;

import java.util.List;

/**
 * The standard nested hierarchy of nucleotide substitution models, each a
 * constrained special case of GTR (Section: see class javadoc on
 * {@link RateParameterization}):
 * <ul>
 *   <li><b>JC69</b> -- equal frequencies, equal rates. 0 free substitution params.</li>
 *   <li><b>F81</b> -- empirical frequencies, equal rates. 3 free params (frequencies).</li>
 *   <li><b>K80</b> -- equal frequencies, one transition/transversion ratio kappa. 1 free param.</li>
 *   <li><b>HKY85</b> -- empirical frequencies + kappa. 4 free params. Arguably the most
 *       commonly used "reasonable default" model in practice.</li>
 *   <li><b>TN93</b> -- empirical frequencies + two transition rates (purine AG, pyrimidine
 *       CT, allowed to differ) + shared transversion rate. 5 free params.</li>
 *   <li><b>GTR</b> -- empirical frequencies + all 6 rates free (5 estimated relative to
 *       GT fixed at 1, since only the ratios matter -- GtrModel normalizes overall
 *       scale into the branch lengths anyway). 8 free params.</li>
 * </ul>
 * Rate array order throughout this package: AC, AG, AT, CG, CT, GT.
 */
public final class NamedSubstitutionModels {

    private NamedSubstitutionModels() {
    }

    public static RateParameterization jc69() {
        return new RateParameterization("JC69", false, new double[0], new double[0], new double[0],
                theta -> new double[]{1, 1, 1, 1, 1, 1});
    }

    public static RateParameterization f81() {
        return new RateParameterization("F81", true, new double[0], new double[0], new double[0],
                theta -> new double[]{1, 1, 1, 1, 1, 1});
    }

    public static RateParameterization k80() {
        return new RateParameterization("K80", false, new double[]{2.0}, new double[]{0.01}, new double[]{100.0},
                theta -> new double[]{1, theta[0], 1, 1, theta[0], 1});
    }

    public static RateParameterization hky85() {
        return new RateParameterization("HKY85", true, new double[]{2.0}, new double[]{0.01}, new double[]{100.0},
                theta -> new double[]{1, theta[0], 1, 1, theta[0], 1});
    }

    public static RateParameterization tn93() {
        return new RateParameterization("TN93", true, new double[]{2.0, 2.0}, new double[]{0.01, 0.01}, new double[]{100.0, 100.0},
                theta -> new double[]{1, theta[0], 1, 1, theta[1], 1});
    }

    public static RateParameterization gtr() {
        return new RateParameterization("GTR", true,
                new double[]{1, 1, 1, 1, 1}, new double[]{0.01, 0.01, 0.01, 0.01, 0.01}, new double[]{100, 100, 100, 100, 100},
                theta -> new double[]{theta[0], theta[1], theta[2], theta[3], theta[4], 1.0});
    }

    /** All 6 models, simplest (fewest parameters) first -- a sensible order for model-selection scans. */
    public static List<RateParameterization> all() {
        return List.of(jc69(), f81(), k80(), hky85(), tn93(), gtr());
    }

    public static RateParameterization byName(String name) {
        for (RateParameterization m : all()) {
            if (m.modelName().equalsIgnoreCase(name)) return m;
        }
        throw new IllegalArgumentException("Unknown substitution model '" + name + "'; expected one of: jc69, f81, k80, hky85, tn93, gtr");
    }
}
