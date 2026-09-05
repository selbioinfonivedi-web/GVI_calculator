package org.gvi.composite;

import org.gvi.core.spi.IndexResult;

import java.util.List;

record FakeIndexResult(String indexName, double primaryValue) implements IndexResult {

    @Override
    public String category() {
        return "test";
    }

    @Override
    public List<String> diagnostics() {
        return List.of();
    }
}
