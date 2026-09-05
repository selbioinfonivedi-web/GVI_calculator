package org.gvi.composite;

import java.util.List;

/** Composite GVI(t) (Section 5.9): weighted, normalized combination of whichever indices were available. */
public record GviResult(double gvi, List<GviComponent> components, List<IndexKey> excludedIndices,
                         double effectiveWeightSum, List<String> diagnostics) {

    /**
     * Fraction of the full weighting scheme that must have usable data before one dataset's GVI can be
     * compared against another's.
     * <p>
     * {@link #effectiveWeightSum} is already the share of the scheme that had data ({@code Sigma(w)} over
     * the contributing indices, out of 1.0), and the engine renormalizes the survivors to sum to 1 so the
     * result always lands in [0,1] and always <em>looks</em> like a GVI. Nothing consumed that number, so a
     * score computed from a sliver of the scheme was published in the same column as a nearly-complete one.
     * <p>
     * That is not a cosmetic problem, because the renormalization inverts the meaning of the score. A real
     * case from this project's corpus: enterotoxaemia failed quality gating on seven of nine indices and its
     * GVI was built from RI and GC alone -- {@code effectiveWeightSum = 0.057}, under 6% of the scheme. RI
     * happened to be 0.0, and RI carried 56% of the renormalized weight, so the reported GVI was 0.053. Next
     * to haemorrhagic septicaemia's 0.680 that reads as "much lower risk", when what actually happened is
     * that the data was too broken to score at all. A near-zero GVI from near-zero coverage is the most
     * dangerous output this tool can produce: it is indistinguishable, by value alone, from a confident
     * finding of low virulence.
     * <p>
     * 0.60 is a reporting threshold, not a scientific constant -- it is set just below the combined weight of
     * the four indices that dominate the scheme (Re, MB, dN/dS, mu, together ~0.73), so a score is called
     * comparable only when most of that backbone survived gating.
     */
    public static final double MIN_COMPARABLE_COVERAGE = 0.60;

    /**
     * True when enough of the weighting scheme had usable data for this score to be compared against
     * another dataset's. When false, the GVI is still a valid summary of the indices that survived, but it
     * is a different quantity from a fully-populated GVI and must not be ranked against one.
     */
    public boolean comparable() {
        return effectiveWeightSum >= MIN_COMPARABLE_COVERAGE;
    }

    /** Human-readable coverage statement, e.g. "6 of 9 indices, 78% of the weighting scheme". */
    public String coverageSummary() {
        int contributing = components.size();
        return String.format(java.util.Locale.ROOT, "%d of %d indices, %.0f%% of the weighting scheme",
                contributing, contributing + excludedIndices.size(), effectiveWeightSum * 100);
    }

    /**
     * beta_genomic(t) = beta0 * [1 + GVI_genomic(t) * scaleFactor], where GVI_genomic
     * excludes the Re component (Section 5.9 / "Use in Epidemic Forecasting").
     * <p>
     * <b>Why Re is excluded here specifically.</b> In any compartmental model, the
     * transmission rate and the reproduction number are related by definition --
     * {@code Re = beta * infectious_period}. Feeding a Re-inclusive index into a formula
     * that produces a new beta means the downstream model then derives its own Re from
     * that inflated beta, counting the same transmissibility signal twice: once as an
     * input to GVI, and again through the beta it raised. Re carries the single largest
     * configured weight, so the double-count was not marginal.
     * <p>
     * The composite GVI reported everywhere else is unchanged and still includes Re --
     * only the value feeding this beta equation drops it, and the remaining weights are
     * renormalized exactly as they are for any other excluded index.
     */
    public double betaGenomic(double beta0, double scaleFactor) {
        return beta0 * (1.0 + gviExcludingRe() * scaleFactor);
    }

    /**
     * The composite score recomputed over every component except Re, with the surviving
     * weights renormalized to sum to 1. Equals {@link #gvi()} when Re did not contribute.
     */
    public double gviExcludingRe() {
        double reWeight = components.stream()
                .filter(c -> c.key() == IndexKey.RE)
                .mapToDouble(GviComponent::effectiveWeight)
                .sum();
        if (reWeight <= 0.0) return gvi;

        double remaining = 1.0 - reWeight;
        if (remaining <= 0.0) return 0.0; // Re was the only contributing index; no genomic signal remains

        double sum = 0.0;
        for (GviComponent c : components) {
            if (c.key() == IndexKey.RE) continue;
            sum += (c.effectiveWeight() / remaining) * c.normalizedValue();
        }
        return sum;
    }

    /** True when Re contributed to {@link #gvi()}, so beta_genomic is computed from a different score. */
    public boolean betaExcludesRe() {
        return components.stream().anyMatch(c -> c.key() == IndexKey.RE && c.effectiveWeight() > 0.0);
    }
}
