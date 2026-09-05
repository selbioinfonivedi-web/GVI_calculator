package org.gvi.cli;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

/**
 * Turns raw quality-gate exclusion messages into a grouped, actionable explanation.
 * <p>
 * The exclusions are the most useful thing this tool produces -- each names a specific, usually
 * fixable problem with the input -- but as emitted they are one long list of prose, in which five
 * indices failing for a single bad alignment reads as five separate problems. Grouping by cause
 * matches what the user actually has to do about it: one alignment to fix, not five indices to
 * chase.
 * <p>
 * Lives here rather than in the desktop UI so the CLI report, the desktop app and the web front end
 * all describe a failure the same way. A user who sees "the alignment is 38% gap characters" in one
 * place and a different phrasing elsewhere has to work out whether they are the same problem.
 */
public final class ExclusionExplainer {

    /** One cause, the indices it took out, and what to do about it. */
    public record Group(String cause, List<String> affectedIndices, String remedy) {
    }

    private ExclusionExplainer() {
    }

    /** Groups the {@code skipped} entries that represent composite exclusions. Order is stable. */
    public static List<Group> explain(List<String> skipped) {
        Map<String, List<String>> byCause = new LinkedHashMap<>();
        for (String s : skipped) {
            if (s == null || !s.contains("excluded from composite")) continue;
            int paren = s.indexOf(" (");
            String index = paren > 0 ? s.substring(0, paren) : s;
            byCause.computeIfAbsent(causeOf(s), c -> new ArrayList<>()).add(index);
        }
        List<Group> out = new ArrayList<>();
        byCause.forEach((cause, indices) -> out.add(new Group(cause, List.copyOf(indices), remedyFor(cause))));
        return out;
    }

    /** The human-readable cause behind one exclusion message. */
    public static String causeOf(String reason) {
        String r = reason == null ? "" : reason.toLowerCase(Locale.ROOT);
        if (r.contains("gap characters") || r.contains("mutation burden is")) {
            return "The alignment is not the closely-related, genuinely-aligned set these indices assume";
        }
        if (r.contains("date-randomization")) {
            return "Collection dates do not explain the divergence (date-randomization test failed)";
        }
        if (r.contains("temporal signal is too weak")) {
            return "Temporal signal too weak to read a molecular clock";
        }
        if (r.contains("implausible for")) {
            return "Fitted rate is implausible for this genome type";
        }
        if (r.contains("non-positive")) {
            return "Fitted substitution rate is not a physically meaningful quantity";
        }
        if (r.contains("lineages-through-time")) {
            return "Re came from the weak fallback estimator";
        }
        if (r.contains("neutrality") || r.contains("z-test")) {
            return "dN/dS above 1 is not statistically distinguishable from neutral evolution";
        }
        return "Other quality-gate failure";
    }

    /** What to actually do about a given cause. */
    public static String remedyFor(String cause) {
        String c = cause == null ? "" : cause.toLowerCase(Locale.ROOT);
        if (c.contains("not the closely-related")) {
            return "If the sequences have different coverage windows, enable auto mode or trim to the shared coverage "
                    + "window. This is usually partial coverage rather than misalignment, so re-running MAFFT will not "
                    + "help -- the aligner already placed them correctly. If they are genuinely too divergent, they are "
                    + "probably not one gene from one organism.";
        }
        if (c.contains("date-randomization")) {
            return "The alignment most likely pools independently-introduced lineages. The report's \"mu remedy\" notice "
                    + "lists the detected groups -- split by those and estimate a rate within each.";
        }
        if (c.contains("temporal signal too weak")) {
            return "Widen the sampling window. Sequences spanning under a couple of years rarely support a clock.";
        }
        if (c.contains("implausible for this genome type")) {
            return "Check the genome type is set correctly, and that the alignment holds one lineage. A DNA genome "
                    + "evolving at RNA-virus rates usually means pooled lineages rather than a genuinely fast genome.";
        }
        if (c.contains("physically meaningful")) {
            return "Check the collection dates parsed correctly, and that the alignment holds a single lineage.";
        }
        if (c.contains("weak fallback")) {
            return "Supply case-incidence data for a true epidemiological Re (Cori et al.), or use a dataset where the "
                    + "birth-death fit converges.";
        }
        if (c.contains("neutral")) {
            return "Confirm the reading frame and gene boundaries with a real GFF3. A ratio above 1 on few codons is at "
                    + "least as likely to be a frameshift artifact as real selection.";
        }
        return "See the full report for the complete reason.";
    }
}
