package org.gvi.ui.theme;

import org.gvi.composite.IndexKey;

import java.util.EnumMap;
import java.util.Map;

/**
 * Short, static explainer text per index -- authored copy, not a computed value. Kept in its own class
 * specifically so it is never confused with anything {@link org.gvi.cli.PipelineResult} actually
 * measured: a metric card shows a real number from a real run, and separately, always visually distinct,
 * one of these fixed sentences describing what that number means in general.
 */
public final class MetricDescriptions {

    private static final Map<IndexKey, String> DESCRIPTIONS = new EnumMap<>(IndexKey.class);

    static {
        DESCRIPTIONS.put(IndexKey.MU, "Evolutionary rate: substitutions per site per year, estimated from "
                + "root-to-tip divergence (or maximum-likelihood dating) along a phylogenetic tree.");
        DESCRIPTIONS.put(IndexKey.RE, "Effective reproduction number: estimated from case-incidence data "
                + "(Cori et al. 2013), or from the tree alone (Euler-Lotka / BDSKY) when no incidence data is supplied.");
        DESCRIPTIONS.put(IndexKey.PI, "Nucleotide diversity (π): mean pairwise nucleotide differences per "
                + "site across the sequence set.");
        DESCRIPTIONS.put(IndexKey.MB, "Mutation burden: count of variant positions (each SNP or indel counts as "
                + "one event) relative to the reference sequence.");
        DESCRIPTIONS.put(IndexKey.DNDS, "Selection pressure (dN/dS): ratio of nonsynonymous to synonymous "
                + "substitution rates. Below 1 indicates purifying selection, above 1 indicates positive selection.");
        DESCRIPTIONS.put(IndexKey.GD, "Genetic distance: pairwise sequence divergence from the reference "
                + "(Jukes-Cantor, Kimura 2-parameter, or Hamming, depending on the selected method).");
        DESCRIPTIONS.put(IndexKey.CAI, "Codon adaptation index: how closely codon usage matches a reference "
                + "host codon-usage table (Sharp & Li 1987). Closer to 1 means better-adapted codon usage.");
        DESCRIPTIONS.put(IndexKey.GC, "GC content deviation from the reference, with an optional test for "
                + "APOBEC3/ADAR mutational-signature enrichment.");
        DESCRIPTIONS.put(IndexKey.RI, "Recombination index: a PHI test (Bruen, Bryant & Poss 2006) for "
                + "phylogenetic incompatibility signal consistent with recombination.");
    }

    private MetricDescriptions() {
    }

    public static String describe(IndexKey key) {
        return DESCRIPTIONS.getOrDefault(key, "");
    }
}
