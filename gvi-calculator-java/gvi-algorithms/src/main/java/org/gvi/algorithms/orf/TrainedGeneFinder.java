package org.gvi.algorithms.orf;

import org.gvi.algorithms.common.CodonUtil;
import org.gvi.algorithms.dnds.ml.CodonAlphabet;
import org.gvi.core.model.GeneAnnotation;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Self-training gene finder -- the same bootstrapping principle Prodigal (Hyatt et al. 2010,
 * "Prodigal: prokaryotic gene recognition and translation initiation site identification") uses to
 * call genes with no external training corpus: long open reading frames are exponentially unlikely
 * to occur by chance in a non-coding region regardless of species, so they can be used to train a
 * species-specific (here, genome-specific) coding-potential model directly from the input sequence
 * itself, then every candidate ORF is scored against that trained model instead of being kept or
 * dropped by length alone.
 * <p>
 * <b>What this does.</b>
 * <ol>
 *   <li>Enumerate every candidate ORF with {@link OrfFinder} (unchanged: 6-frame, stop-to-stop
 *       segmentation, first in-frame ATG defines each segment's start).</li>
 *   <li>Treat every candidate at least {@link #MIN_TRAINING_ORF_CODONS} codons long as a confident
 *       real gene and train a codon-usage frequency table from their internal codons (i.e. excluding
 *       each one's own start and stop codon, which are structural rather than compositional signal).</li>
 *   <li>Train a background ("what non-coding sequence looks like") model from the whole reference
 *       sequence's own single-nucleotide composition: {@code P(codon) = P(base1)*P(base2)*P(base3)},
 *       independence assumed -- a GC-content-only null model, i.e. exactly Prodigal's own generic
 *       starting model before its own self-training refines it further.</li>
 *   <li>Score every candidate (not just the training set) by its mean per-codon log-odds,
 *       {@code log(P_coding(codon) / P_background(codon))} averaged over its internal codons, and
 *       keep those that score positive -- i.e. more consistent with the trained coding model than
 *       with the genome's own bulk composition.</li>
 * </ol>
 * <p>
 * <b>What this deliberately does not do</b>, relative to genuine Prodigal: no higher-order
 * (context-dependent) Markov chain -- codon usage here is a simple 0th-order frequency table, not a
 * periodic Markov model conditioned on neighbouring bases; no ribosome-binding-site motif model, so
 * translation initiation site selection is unchanged from {@link OrfFinder}'s "first in-frame ATG"
 * convention rather than refined against Shine-Dalgarno signal; no dynamic programming across
 * overlapping candidate genes. This scores and filters candidate gene CALLS with a real, genome-self-trained
 * coding-potential model -- a real accuracy upgrade over a length-only heuristic -- without
 * attempting genuine ab initio start-site prediction.
 * <p>
 * Scoped to intron-free genomes (bacteria and viruses, per {@code --organism-class}): like {@link
 * OrfFinder}, there is no splice-site handling.
 */
public final class TrainedGeneFinder {

    /**
     * Prodigal-inspired floor for "confident enough to self-train from without external data" --
     * a run of this many uninterrupted sense codons is vanishingly unlikely by chance (roughly
     * (61/64)^270 for 90 codons under a uniform-random null), so treating it as a real gene to learn
     * this genome's codon bias from does not require already knowing which sequences are genes.
     */
    public static final int MIN_TRAINING_ORF_CODONS = 90;

    /** Below this many training-quality ORFs there is not enough signal to self-train reliably; fall back to the plain length heuristic. */
    private static final int MIN_TRAINING_SET_SIZE = 3;

    /** A candidate scores as coding when it fits the trained model better than the background model -- the natural zero-point of a log-odds ratio. */
    private static final double CODING_LOG_ODDS_THRESHOLD = 0.0;

    /** Laplace smoothing so a codon absent from a small training set is merely rare, not impossible (which would make any candidate using it score -Infinity). */
    private static final double PSEUDOCOUNT = 1.0;

    /** {@link OrfFinder}'s own relative-length heuristic, used only as the fallback when there isn't enough signal to self-train (see {@link #MIN_TRAINING_SET_SIZE}). */
    private static final double FALLBACK_RELATIVE_LENGTH_THRESHOLD = 0.3;

    public record ScoredGene(GeneAnnotation gene, double codingLogOdds) {
    }

    /**
     * @return genes this run's trained model calls as coding, longest first; empty only when {@link
     * OrfFinder} itself found no candidate at all (mirrors {@link OrfFinder#findOrfs}'s own contract).
     */
    public List<GeneAnnotation> findGenes(String referenceSequence) {
        List<ScoredGene> scored = scoreAll(referenceSequence);
        if (scored.isEmpty()) return List.of();
        if (Double.isNaN(scored.get(0).codingLogOdds())) {
            // Sentinel from scoreAll for "not enough training signal" (every entry is NaN in that
            // case, never a mix -- see scoreAll's javadoc): fall back to the plain length heuristic.
            return fallbackByRelativeLength(scored.stream().map(ScoredGene::gene).toList());
        }

        List<GeneAnnotation> kept = scored.stream()
                .filter(s -> s.codingLogOdds() > CODING_LOG_ODDS_THRESHOLD)
                .map(ScoredGene::gene)
                .sorted(Comparator.comparingLong(GeneAnnotation::length).reversed())
                .toList();
        if (!kept.isEmpty()) return kept;

        // Nothing cleared the coding-potential bar (e.g. a very short or compositionally flat
        // sequence) -- report the single best-scoring candidate rather than nothing, the same
        // "don't report an empty gene set when candidates exist" principle OrfFinder's own caller
        // applies when even the length heuristic finds nothing.
        return List.of(scored.stream().max(Comparator.comparingDouble(ScoredGene::codingLogOdds)).orElseThrow().gene());
    }

    /**
     * Every candidate ORF {@link OrfFinder} found, each scored against this genome's own
     * self-trained coding-potential model -- exposed separately from {@link #findGenes} so a caller
     * (tests, or a future diagnostics surface) can inspect why a given region was kept or rejected,
     * not just the final kept/rejected list.
     *
     * @return empty only when {@link OrfFinder} found no candidate at all; when there was not enough
     * long-ORF signal to self-train (see {@link #MIN_TRAINING_SET_SIZE}), returns every raw candidate
     * with {@link ScoredGene#codingLogOdds()} as {@link Double#NaN} rather than a trained score.
     */
    public List<ScoredGene> scoreAll(String referenceSequence) {
        List<GeneAnnotation> candidates = new OrfFinder(OrfFinder.DEFAULT_MIN_ORF_CODONS).findOrfs(referenceSequence);
        if (candidates.isEmpty()) return List.of();

        List<GeneAnnotation> trainingSet = candidates.stream()
                .filter(g -> g.length() / 3 >= MIN_TRAINING_ORF_CODONS)
                .toList();
        if (trainingSet.size() < MIN_TRAINING_SET_SIZE) {
            return candidates.stream().map(g -> new ScoredGene(g, Double.NaN)).toList();
        }

        Map<String, Double> codingFreq = trainCodingModel(referenceSequence, trainingSet);
        Map<String, Double> backgroundFreq = trainBackgroundModel(referenceSequence);

        List<ScoredGene> scored = new ArrayList<>(candidates.size());
        for (GeneAnnotation g : candidates) {
            scored.add(new ScoredGene(g, codingLogOdds(referenceSequence, g, codingFreq, backgroundFreq)));
        }
        return scored;
    }

    /** {@link OrfFinder}'s original filter, used only when there isn't enough long-ORF signal in this genome to self-train a coding model. */
    private List<GeneAnnotation> fallbackByRelativeLength(List<GeneAnnotation> candidates) {
        long longest = candidates.stream().mapToLong(GeneAnnotation::length).max().orElseThrow();
        return candidates.stream()
                .filter(g -> g.length() >= longest * FALLBACK_RELATIVE_LENGTH_THRESHOLD)
                .sorted(Comparator.comparingLong(GeneAnnotation::length).reversed())
                .toList();
    }

    /** Codon-usage frequency table trained on the training set's own internal (non-start, non-stop) codons, Laplace-smoothed over the 61 sense codons. */
    private Map<String, Double> trainCodingModel(String referenceSequence, List<GeneAnnotation> trainingSet) {
        Map<String, Double> counts = new HashMap<>();
        for (String codon : CodonAlphabet.SENSE_CODONS) counts.put(codon, PSEUDOCOUNT);

        for (GeneAnnotation g : trainingSet) {
            for (String codon : internalCodons(referenceSequence, g)) {
                counts.merge(codon, 1.0, Double::sum);
            }
        }
        return normalize(counts);
    }

    /** GC-content-only null model: each codon's background probability is the product of its 3 bases' genome-wide frequencies, assuming independence. */
    private Map<String, Double> trainBackgroundModel(String referenceSequence) {
        Map<Character, Long> baseCounts = new HashMap<>();
        for (char c : referenceSequence.toUpperCase().toCharArray()) {
            if (c == 'A' || c == 'C' || c == 'G' || c == 'T') baseCounts.merge(c, 1L, Long::sum);
        }
        long total = baseCounts.values().stream().mapToLong(Long::longValue).sum();
        Map<Character, Double> baseFreq = new HashMap<>();
        for (char base : new char[]{'A', 'C', 'G', 'T'}) {
            baseFreq.put(base, (baseCounts.getOrDefault(base, 0L) + PSEUDOCOUNT) / (total + 4 * PSEUDOCOUNT));
        }

        Map<String, Double> background = new HashMap<>();
        for (String codon : CodonAlphabet.SENSE_CODONS) {
            double p = 1.0;
            for (char c : codon.toCharArray()) p *= baseFreq.get(c);
            background.put(codon, p);
        }
        return background;
    }

    /** Mean per-codon log-odds of this candidate's internal codons under the trained coding model vs. the background model. */
    private double codingLogOdds(String referenceSequence, GeneAnnotation gene, Map<String, Double> codingFreq, Map<String, Double> backgroundFreq) {
        List<String> internal = internalCodons(referenceSequence, gene);
        if (internal.isEmpty()) return Double.NEGATIVE_INFINITY; // nothing but a start and stop codon: cannot be scored, so cannot pass
        double sum = 0.0;
        for (String codon : internal) {
            sum += Math.log(codingFreq.get(codon) / backgroundFreq.get(codon));
        }
        return sum / internal.size();
    }

    /** This gene's codons excluding its own first (start) and last (stop) codon, and any codon containing a non-ACGT base. */
    private List<String> internalCodons(String referenceSequence, GeneAnnotation gene) {
        String cds = CodonUtil.extractCds(referenceSequence, gene.start(), gene.end(), gene.strand());
        List<String> codons = CodonUtil.splitIntoCodons(cds).codons();
        List<String> internal = new ArrayList<>();
        for (int i = 1; i < codons.size() - 1; i++) {
            String codon = codons.get(i).toUpperCase();
            if (isStandardSenseCodon(codon)) internal.add(codon);
        }
        return internal;
    }

    private boolean isStandardSenseCodon(String codon) {
        return codon.length() == 3 && codon.chars().allMatch(c -> c == 'A' || c == 'C' || c == 'G' || c == 'T')
                && CodonAlphabet.SENSE_CODONS.contains(codon);
    }

    private Map<String, Double> normalize(Map<String, Double> counts) {
        double total = counts.values().stream().mapToDouble(Double::doubleValue).sum();
        Map<String, Double> freq = new HashMap<>();
        for (var entry : counts.entrySet()) freq.put(entry.getKey(), entry.getValue() / total);
        return freq;
    }
}
