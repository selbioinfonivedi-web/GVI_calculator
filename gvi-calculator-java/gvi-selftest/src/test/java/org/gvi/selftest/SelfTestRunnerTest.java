package org.gvi.selftest;

import org.junit.jupiter.api.Test;

import static org.assertj.core.api.Assertions.assertThat;

class SelfTestRunnerTest {

    @Test
    void allBundledSelfTestsPassOnThisInstallation() {
        SelfTestReport report = new SelfTestRunner().runAll();
        for (SelfTestCase c : report.cases()) {
            assertThat(c.passed()).as("%s: %s", c.name(), c.detail()).isTrue();
        }
        assertThat(report.allPassed()).isTrue();
        assertThat(report.cases()).hasSize(13);
    }
}
