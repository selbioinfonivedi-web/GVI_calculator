package org.gvi.selftest;

import java.util.List;

public record SelfTestReport(List<SelfTestCase> cases, boolean allPassed) {

    public static SelfTestReport of(List<SelfTestCase> cases) {
        boolean all = cases.stream().allMatch(SelfTestCase::passed);
        return new SelfTestReport(cases, all);
    }
}
