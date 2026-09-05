package org.gvi.core.util;

/**
 * Guards O(n^2) pairwise computations (nucleotide diversity, genetic
 * distance matrices) against exhausting the JVM heap on large sequence
 * sets (Section 2 of the build spec: "never let the JVM hit OutOfMemoryError
 * unhandled"). Callers ask {@link #shouldSubsample(int, long)} before
 * starting a full pairwise loop and fall back to random pair subsampling if
 * it returns true.
 */
public final class MemoryGuard {

    /** Fraction of max heap we allow a single pairwise computation to plan on using. */
    private static final double SAFE_HEAP_FRACTION = 0.5;

    private MemoryGuard() {
    }

    /**
     * @param sequenceCount   number of sequences that would be compared pairwise
     * @param bytesPerPairwiseComparison estimated bytes of transient state per pair
     *                                    (e.g. diff-counting over the alignment length)
     */
    public static boolean shouldSubsample(int sequenceCount, long bytesPerPairwiseComparison) {
        long pairs = pairCount(sequenceCount);
        long estimatedBytes = pairs * bytesPerPairwiseComparison;
        long safeHeap = (long) (Runtime.getRuntime().maxMemory() * SAFE_HEAP_FRACTION);
        return estimatedBytes > safeHeap || sequenceCount > 2000;
    }

    public static long pairCount(int n) {
        return (long) n * (n - 1) / 2;
    }
}
