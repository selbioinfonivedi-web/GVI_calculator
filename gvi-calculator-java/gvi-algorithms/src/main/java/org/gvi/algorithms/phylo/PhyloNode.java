package org.gvi.algorithms.phylo;

import java.util.ArrayList;
import java.util.List;

/**
 * One node of a rooted tree. {@code label} is non-null only for nodes
 * corresponding to an original input sequence ("taxa"); Neighbor-Joining's
 * internal branch points have {@code label == null}. {@code branchLength}
 * is the length of the edge connecting this node to its parent
 * (substitutions/site), 0 for the root.
 */
public final class PhyloNode {

    final int id;
    private final String label;
    private double branchLength;
    private final List<PhyloNode> children = new ArrayList<>();
    private PhyloNode parent;

    PhyloNode(int id, String label) {
        this.id = id;
        this.label = label;
    }

    /**
     * Test-support factory: builds a standalone, detached node for
     * constructing hand-specified trees in tests outside this package
     * (production trees are only ever built via {@link PhyloTreeFactory};
     * the normal constructor is package-private to enforce that).
     */
    public static PhyloNode createDetachedForTesting(int id, String label) {
        return new PhyloNode(id, label);
    }

    /** Test-support: attaches {@code child} to this node with the given branch length -- see {@link #createDetachedForTesting}. */
    public void attachChildForTesting(PhyloNode child, double branchLength) {
        child.setParent(this, branchLength);
        this.addChild(child);
    }

    void setParent(PhyloNode parent, double branchLength) {
        this.parent = parent;
        this.branchLength = branchLength;
    }

    void addChild(PhyloNode child) {
        children.add(child);
    }

    /**
     * Moves this node to be a child of {@code newParent}, detaching it from its
     * current parent's child list first -- used by topology-search moves (NNI)
     * that need to rearrange an already-built tree. This node's own
     * {@code branchLength} (the length of the edge to whichever node is its
     * parent) is left untouched; callers that want a fresh length for the new
     * edge context should call {@link #setBranchLength} afterward.
     */
    public void reparent(PhyloNode newParent) {
        if (parent != null) {
            parent.children.remove(this);
        }
        this.parent = newParent;
        newParent.children.add(this);
    }

    public boolean isTaxon() {
        return label != null;
    }

    public boolean isRoot() {
        return parent == null;
    }

    public String label() {
        return label;
    }

    public double branchLength() {
        return branchLength;
    }

    /** Mutates this node's branch length in place -- used by ML branch-length optimization after the tree is built. */
    public void setBranchLength(double branchLength) {
        this.branchLength = branchLength;
    }

    public List<PhyloNode> children() {
        return children;
    }

    public PhyloNode parent() {
        return parent;
    }

    /** Sum of branch lengths from the tree's root down to this node (substitutions/site). */
    public double rootToTipDistance() {
        double sum = 0;
        PhyloNode n = this;
        while (n.parent != null) {
            sum += n.branchLength;
            n = n.parent;
        }
        return sum;
    }
}
