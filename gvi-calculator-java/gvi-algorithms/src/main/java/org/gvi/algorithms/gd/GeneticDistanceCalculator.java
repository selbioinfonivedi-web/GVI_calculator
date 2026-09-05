package org.gvi.algorithms.gd;

import org.gvi.core.exception.GviComputationException;
import org.gvi.core.model.NucleotideSequence;
import org.gvi.core.model.SequenceAlignment;
import org.gvi.core.util.IupacUtil;

import java.util.ArrayList;
import java.util.List;

/**
 * Index 6 - Genetic Distance (Section 5.6). Implements Hamming, Jukes-Cantor
 * (JC69), and Kimura 2-Parameter (K80). Every pairwise comparison uses
 * pairwise deletion (gap/N/ambiguous positions in either sequence are
 * excluded from both numerator and denominator, not treated as
 * differences).
 */
public final class GeneticDistanceCalculator {

    /** Site counts needed to evaluate any of the three distance formulas from one pass over the alignment. */
    private record SiteCounts(long comparableSites, long differences, long transitions, long transversions) {
    }

    public GdResult compute(NucleotideSequence a, NucleotideSequence b, GdMethod method) {
        if (a.length() != b.length()) {
            throw new GviComputationException("Cannot compute genetic distance: '" + a.getId() + "' (len="
                    + a.length() + ") and '" + b.getId() + "' (len=" + b.length() + ") are not aligned to equal length");
        }
        SiteCounts counts = countSites(a.getSequence(), b.getSequence());
        if (counts.comparableSites() == 0) {
            throw new GviComputationException("No comparable (unambiguous) sites between '" + a.getId() + "' and '" + b.getId() + "'");
        }

        double p = (double) counts.differences() / counts.comparableSites();
        List<String> diagnostics = new ArrayList<>();
        diagnostics.add(counts.comparableSites() + " of " + a.length() + " aligned positions were comparable (rest were gaps/ambiguous)");

        return switch (method) {
            case HAMMING -> {
                String cat = "N/A (reference bands defined for Jukes-Cantor distance only)";
                yield new GdResult(a.getId(), b.getId(), method, p, p, false, cat, diagnostics);
            }
            case JUKES_CANTOR -> buildJcResult(a.getId(), b.getId(), p, diagnostics);
            case KIMURA_2_PARAMETER -> buildK80Result(a.getId(), b.getId(), counts, diagnostics);
        };
    }

    private GdResult buildJcResult(String idA, String idB, double p, List<String> diagnostics) {
        double arg = 1.0 - (4.0 / 3.0) * p;
        if (arg <= 0.0) {
            diagnostics.add("Jukes-Cantor distance saturated (p-distance too high for the model); reporting raw p-distance instead");
            return new GdResult(idA, idB, GdMethod.JUKES_CANTOR, p, null, true, "Saturated / unreliable, use raw p-distance", diagnostics);
        }
        double d = org.gvi.core.util.MathUtil.stripNegativeZero(-0.75 * Math.log(arg));
        GdReferenceTable.Band band = GdReferenceTable.classify(d);
        String category = band.classification() + " (" + band.timeSinceDivergence() + ", confidence: " + band.confidence() + ")";
        return new GdResult(idA, idB, GdMethod.JUKES_CANTOR, p, d, false, category, diagnostics);
    }

    private GdResult buildK80Result(String idA, String idB, SiteCounts counts, List<String> diagnostics) {
        double bigP = (double) counts.transitions() / counts.comparableSites();
        double bigQ = (double) counts.transversions() / counts.comparableSites();
        double arg1 = 1.0 - 2.0 * bigP - bigQ;
        double arg2 = 1.0 - 2.0 * bigQ;
        if (arg1 <= 0.0 || arg2 <= 0.0) {
            double p = (double) counts.differences() / counts.comparableSites();
            diagnostics.add("Kimura 2-parameter distance saturated (transition/transversion proportions too high for the model); reporting raw p-distance instead");
            return new GdResult(idA, idB, GdMethod.KIMURA_2_PARAMETER, p, null, true, "Saturated / unreliable, use raw p-distance", diagnostics);
        }
        double d = org.gvi.core.util.MathUtil.stripNegativeZero(-0.5 * Math.log(arg1) - 0.25 * Math.log(arg2));
        double p = (double) counts.differences() / counts.comparableSites();
        return new GdResult(idA, idB, GdMethod.KIMURA_2_PARAMETER, p, d, false,
                "N/A (reference bands defined for Jukes-Cantor distance only)", diagnostics);
    }

    private SiteCounts countSites(String seqA, String seqB) {
        long comparable = 0, diffs = 0, ts = 0, tv = 0;
        int len = seqA.length();
        for (int i = 0; i < len; i++) {
            char x = seqA.charAt(i);
            char y = seqB.charAt(i);
            if (!IupacUtil.isComparable(x) || !IupacUtil.isComparable(y)) continue;
            comparable++;
            if (Character.toUpperCase(x) != Character.toUpperCase(y)) {
                diffs++;
                if (IupacUtil.isTransition(x, y)) ts++;
                else if (IupacUtil.isTransversion(x, y)) tv++;
            }
        }
        return new SiteCounts(comparable, diffs, ts, tv);
    }

    /** Distance from every query sequence to the alignment's designated reference (Section 6.3 use case: variant classification). */
    public List<GdResult> distancesToReference(SequenceAlignment alignment, GdMethod method) {
        NucleotideSequence ref = alignment.getReference();
        List<GdResult> results = new ArrayList<>();
        for (NucleotideSequence query : alignment.getQueries()) {
            results.add(compute(ref, query, method));
        }
        return results;
    }
}
