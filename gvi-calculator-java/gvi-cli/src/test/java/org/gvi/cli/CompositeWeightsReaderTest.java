package org.gvi.cli;

import org.gvi.composite.CompositeWeights;
import org.gvi.composite.IndexKey;
import org.gvi.core.exception.GviInputException;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;

import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.assertThatThrownBy;
import static org.assertj.core.api.Assertions.within;

class CompositeWeightsReaderTest {

    @TempDir
    Path tempDir;

    @Test
    void aCompleteSetOfWeightsIsUsedAsGivenAfterRenormalization() throws IOException {
        Path file = writeJson("{\"mu\": 0.1, \"Re\": 0.2, \"pi\": 0.1, \"MB\": 0.2, \"dN/dS\": 0.1, "
                + "\"GD\": 0.1, \"CAI\": 0.1, \"GC_Deviation\": 0.05, \"RI\": 0.05}");

        CompositeWeightsReader.Result result = CompositeWeightsReader.read(file);

        assertThat(result.warnings()).isEmpty();
        assertThat(result.weights().get(IndexKey.RE)).isCloseTo(0.2, within(1e-9));
        assertThat(result.weights().get(IndexKey.MU)).isCloseTo(0.1, within(1e-9));
        double sum = result.weights().asMap().values().stream().mapToDouble(Double::doubleValue).sum();
        assertThat(sum).isCloseTo(1.0, within(1e-9));
    }

    @Test
    void aPartialOverrideKeepsSpecDefaultsForEverythingElseThenRenormalizes() throws IOException {
        // Only bump Re; everything else should still reflect the spec's relative defaults after renormalization.
        Path file = writeJson("{\"Re\": 0.6}");

        CompositeWeightsReader.Result result = CompositeWeightsReader.read(file);
        CompositeWeights defaults = CompositeWeights.defaults();

        assertThat(result.warnings()).isEmpty();
        // Re's share of the total went up relative to its default share.
        double reShareBefore = defaults.get(IndexKey.RE);
        double reShareAfter = result.weights().get(IndexKey.RE);
        assertThat(reShareAfter).isGreaterThan(reShareBefore);
        // mu/pi ratio should be preserved (both untouched, so their relative proportion doesn't change).
        double ratioBefore = defaults.get(IndexKey.MU) / defaults.get(IndexKey.PI);
        double ratioAfter = result.weights().get(IndexKey.MU) / result.weights().get(IndexKey.PI);
        assertThat(ratioAfter).isCloseTo(ratioBefore, within(1e-9));
        double sum = result.weights().asMap().values().stream().mapToDouble(Double::doubleValue).sum();
        assertThat(sum).isCloseTo(1.0, within(1e-9));
    }

    @Test
    void anUnrecognizedLabelProducesAWarningNotAFailure() throws IOException {
        Path file = writeJson("{\"Re\": 0.3, \"not_a_real_index\": 0.5}");

        CompositeWeightsReader.Result result = CompositeWeightsReader.read(file);

        assertThat(result.warnings()).hasSize(1);
        assertThat(result.warnings().get(0)).contains("not_a_real_index");
        assertThat(result.weights().get(IndexKey.RE)).isGreaterThan(0);
    }

    @Test
    void aNegativeWeightIsIgnoredWithAWarning() throws IOException {
        Path file = writeJson("{\"Re\": -0.1}");

        CompositeWeightsReader.Result result = CompositeWeightsReader.read(file);

        assertThat(result.warnings()).hasSize(1);
        assertThat(result.warnings().get(0)).contains("negative");
        // Falls back to the spec default for Re rather than the rejected negative value.
        assertThat(result.weights().get(IndexKey.RE)).isCloseTo(
                CompositeWeights.defaults().get(IndexKey.RE), within(1e-9));
    }

    @Test
    void aMissingFileFailsClearly() {
        assertThatThrownBy(() -> CompositeWeightsReader.read(tempDir.resolve("does_not_exist.json")))
                .isInstanceOf(GviInputException.class);
    }

    @Test
    void malformedJsonFailsClearly() throws IOException {
        Path file = writeJson("{not valid json");
        assertThatThrownBy(() -> CompositeWeightsReader.read(file))
                .isInstanceOf(GviInputException.class);
    }

    private Path writeJson(String content) throws IOException {
        Path file = tempDir.resolve("weights.json");
        Files.writeString(file, content, StandardCharsets.UTF_8);
        return file;
    }
}
