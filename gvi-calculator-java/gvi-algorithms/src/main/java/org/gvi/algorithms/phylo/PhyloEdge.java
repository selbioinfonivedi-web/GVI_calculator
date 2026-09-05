package org.gvi.algorithms.phylo;

/**
 * One edge of an unrooted tree, connecting two node ids (as assigned by
 * {@link NeighborJoining}, or supplied directly when hand-building a small
 * tree, e.g. for tests). Length is in substitutions/site.
 */
public record PhyloEdge(int a, int b, double length) {
}
