package org.gvi.ui.nav;

/** One leaf item in the left navigation rail: a stable id, its display label, and which group it sits under. */
public record WorkflowStep(String id, String label, NavGroup group) {
}
