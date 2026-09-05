package org.gvi.algorithms.cai;

import org.gvi.core.exception.GviComputationException;
import org.gvi.core.exception.GviException;
import org.gvi.core.model.NucleotideSequence;
import org.gvi.core.util.IupacUtil;

import java.util.ArrayList;
import java.util.List;

/** GC Content Deviation half of Index 8 (Section 5.8): |GC%_observed - GC%_reference|. */
public final class GcContentCalculator {

    public GcResult compute(NucleotideSequence sequence, double referenceGcPercent) {
        String seq = sequence.getSequence();
        long gc = 0, unambiguous = 0;
        for (int i = 0; i < seq.length(); i++) {
            char c = Character.toUpperCase(seq.charAt(i));
            if (!IupacUtil.isUnambiguousBase(c)) continue;
            unambiguous++;
            if (c == 'G' || c == 'C') gc++;
        }
        if (unambiguous == 0) {
            throw new GviComputationException("Cannot compute GC content for '" + sequence.getId() + "': no unambiguous bases");
        }
        double observed = 100.0 * gc / unambiguous;
        double deviation = Math.abs(observed - referenceGcPercent);
        var band = GcReferenceTable.classify(deviation);
        String category = band.status() + " (" + band.interpretation() + "; " + band.evolutionaryAge() + ")";

        List<String> diagnostics = new ArrayList<>();
        diagnostics.add(unambiguous + " of " + seq.length() + " bases were unambiguous and used for GC%");

        return new GcResult(sequence.getId(), observed, referenceGcPercent, deviation, category, diagnostics, null);
    }

    /**
     * Same GC%/deviation computation, plus a mutational-signature scan
     * ({@link MutationalSignatureAnalyzer}) against the actual reference
     * sequence: is the compositional shift explained by known host
     * RNA-editing pressure (APOBEC3 C-&gt;T / ADAR A-&gt;G hotspot
     * enrichment) rather than the index's default "HGT suspected" reading?
     * Falls back to the plain computation (no signature scan, noted in
     * diagnostics) if the reference isn't aligned to equal length.
     */
    public GcResult compute(NucleotideSequence sequence, double referenceGcPercent, NucleotideSequence referenceSequence) {
        GcResult base = compute(sequence, referenceGcPercent);
        List<String> diagnostics = new ArrayList<>(base.diagnostics());
        String category = base.category();
        MutationalSignatureResult signature = null;

        try {
            signature = new MutationalSignatureAnalyzer().analyze(referenceSequence, sequence);
            diagnostics.addAll(signature.diagnostics());
            if (signature.anySignatureDetected()) {
                List<String> detected = new ArrayList<>();
                if (signature.apobecLikeC2T().significant()) detected.add("APOBEC3-like C->T");
                if (signature.adarLikeA2G().significant()) detected.add("ADAR-like A->G");
                category = category + " -- NOTE: compositional shift shows significant "
                        + String.join(" and ", detected)
                        + " hypermutation enrichment -- likely host-immune RNA-editing pressure, not necessarily HGT/reassortment";
            }
        } catch (GviException e) {
            diagnostics.add("Mutational-signature scan skipped: " + e.getMessage());
        }

        return new GcResult(base.sequenceId(), base.observedGcPercent(), base.referenceGcPercent(),
                base.deviationPercent(), category, diagnostics, signature);
    }
}
