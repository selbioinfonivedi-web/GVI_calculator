package org.gvi.algorithms.mu;

/** Which regression {@link MuResult} was built from. */
public enum MuEstimationMethod {
    /** Distance from each dated sequence straight to the reference (Section 5.1's baseline description). */
    PAIRWISE_TO_REFERENCE,
    /** Root-to-tip distance along an actual Neighbor-Joining tree rooted at the reference -- accounts for shared ancestry between sequences instead of comparing each independently to one point. */
    TREE_ROOT_TO_TIP,
    /** Root-to-tip distance along an NJ topology whose branch lengths were re-optimized by maximum likelihood under GTR(+Gamma) -- see {@link EvolutionaryRateCalculator#computeGtrGammaAware}. Slower; opt-in. */
    GTR_GAMMA_ML_BRANCH_LENGTHS,
    /** Joint least-squares fit of the rate AND every internal node's date directly against every edge, under temporal-precedence constraints -- the same objective LSD2 (To et al. 2016) solves. See {@link EvolutionaryRateCalculator#computeLeastSquaresDating} / {@link LeastSquaresDatingEstimator}. */
    LEAST_SQUARES_DATING,
    /** Uncorrelated lognormal relaxed clock -- each branch gets its own rate instead of one shared rate, the same generative model BEAST2's UCLD clock uses, fit here by maximum likelihood rather than MCMC. See {@link EvolutionaryRateCalculator#computeRelaxedClock} / {@link RelaxedClockMlEstimator}. */
    RELAXED_CLOCK_ML
}
