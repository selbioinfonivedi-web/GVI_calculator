package org.gvi.composite;

/** One index's effect on GVI when its weight is perturbed +/-N% (Section 5.9's sensitivity analysis). */
public record SensitivityResult(IndexKey key, double baseGvi, double gviAtLowWeight, double gviAtHighWeight) {

    public double spread() {
        return gviAtHighWeight - gviAtLowWeight;
    }
}
