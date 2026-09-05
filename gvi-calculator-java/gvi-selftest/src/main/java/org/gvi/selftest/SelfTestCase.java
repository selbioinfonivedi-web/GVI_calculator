package org.gvi.selftest;

/** One self-test outcome: which check, whether it passed, and (if not) why. */
public record SelfTestCase(String name, boolean passed, String detail) {
}
