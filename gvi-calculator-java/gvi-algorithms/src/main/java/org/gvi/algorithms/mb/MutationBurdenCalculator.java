package org.gvi.algorithms.mb;

import org.gvi.core.exception.GviComputationException;
import org.gvi.core.model.NucleotideSequence;
import org.gvi.core.model.VariantRecord;
import org.gvi.core.util.IupacUtil;

import java.util.ArrayList;
import java.util.List;

/**
 * Index 4 - Mutation Burden (Section 5.4): count of variant positions
 * relative to a reference (SNP = 1 event, indel of any length = 1 event),
 * excluding sites below the coverage threshold when depth information is
 * available.
 */
public final class MutationBurdenCalculator {

    private static final int DEFAULT_COVERAGE_THRESHOLD = 10;

    private final int coverageThreshold;

    public MutationBurdenCalculator() {
        this(DEFAULT_COVERAGE_THRESHOLD);
    }

    public MutationBurdenCalculator(int coverageThreshold) {
        this.coverageThreshold = coverageThreshold;
    }

    /** Preferred path: VCF variant calls carry per-site depth for a real coverage filter. */
    public MbResult computeFromVariants(String sequenceId, List<VariantRecord> variants) {
        int included = 0, excludedLowCoverage = 0, unknownCoverage = 0;
        for (VariantRecord v : variants) {
            if (v.depth() != null) {
                if (v.depth() < coverageThreshold) {
                    excludedLowCoverage++;
                    continue;
                }
            } else {
                unknownCoverage++;
            }
            included++;
        }
        List<String> diagnostics = new ArrayList<>();
        diagnostics.add(excludedLowCoverage + " variant(s) excluded for depth < " + coverageThreshold + "x");
        if (unknownCoverage > 0) {
            diagnostics.add(unknownCoverage + " variant(s) had no reported depth; included without coverage filtering");
        }
        var band = MbReferenceTable.classify(included);
        String category = band.pathogenExample() + " (" + band.timeline() + "); " + band.controlImplication();
        return new MbResult(sequenceId, included, MbResult.MbSource.VCF_WITH_COVERAGE_FILTER,
                excludedLowCoverage, unknownCoverage, category, diagnostics);
    }

    /**
     * Fallback path when no VCF/depth information is available: diff the
     * query against the reference directly in the alignment. Contiguous
     * runs of the same indel type collapse into a single event; ambiguous
     * or gap-padding positions are excluded from the count (not treated as
     * mutations), and there is no coverage filter (flagged in diagnostics).
     */
    public MbResult computeFromAlignment(NucleotideSequence reference, NucleotideSequence query) {
        if (reference.length() != query.length()) {
            throw new GviComputationException("Cannot compute mutation burden: '" + reference.getId() + "' and '"
                    + query.getId() + "' are not aligned to equal length");
        }
        String ref = reference.getSequence();
        String qry = query.getSequence();
        int mutations = 0;
        int excludedAmbiguous = 0;
        boolean inGapRun = false;
        char gapRunType = 0;

        for (int i = 0; i < ref.length(); i++) {
            char r = ref.charAt(i);
            char q = qry.charAt(i);
            boolean refGap = IupacUtil.isGap(r);
            boolean qryGap = IupacUtil.isGap(q);

            if (!refGap && !qryGap) {
                inGapRun = false;
                if (!IupacUtil.isComparable(r) || !IupacUtil.isComparable(q)) {
                    excludedAmbiguous++;
                    continue;
                }
                if (r != q) mutations++;
            } else if (refGap && qryGap) {
                inGapRun = false; // pure alignment padding, not a real indel
            } else {
                char type = refGap ? 'I' : 'D'; // insertion in query, or deletion from reference
                if (!(inGapRun && gapRunType == type)) {
                    mutations++;
                    inGapRun = true;
                    gapRunType = type;
                }
            }
        }

        List<String> diagnostics = new ArrayList<>();
        diagnostics.add("No VCF/depth data supplied -- coverage filter (>=" + coverageThreshold + "x) was NOT applied to this count");
        if (excludedAmbiguous > 0) {
            diagnostics.add(excludedAmbiguous + " ambiguous/N position(s) excluded from comparison");
        }
        var band = MbReferenceTable.classify(mutations);
        String category = band.pathogenExample() + " (" + band.timeline() + "); " + band.controlImplication();
        return new MbResult(query.getId(), mutations, MbResult.MbSource.ALIGNMENT_NO_COVERAGE_FILTER,
                0, excludedAmbiguous, category, diagnostics);
    }

    /**
     * Secondary/derived output (Section 4.3): approximate circulation time
     * implied by observed mutation burden and an independently estimated
     * substitution rate. Dimensionally: mu is substitutions/site/year, so
     * genome-wide substitutions/year = mu * genomeLength; weeks = (MB / that
     * rate) * 52.1775. Only meaningful when mu was estimated for the same
     * lineage -- callers must not present this as precise.
     */
    public static double estimateOutbreakAgeWeeks(int mutationBurden, double muSubPerSiteYear, int genomeLength) {
        if (muSubPerSiteYear <= 0 || genomeLength <= 0) {
            throw new GviComputationException("Cannot estimate outbreak age: mu and genome length must both be positive");
        }
        double substitutionsPerYear = muSubPerSiteYear * genomeLength;
        double years = mutationBurden / substitutionsPerYear;
        return years * 52.1775;
    }
}
