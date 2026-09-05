package org.gvi.algorithms.re;

import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.dataformat.yaml.YAMLFactory;
import org.gvi.core.exception.GviInputException;

import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

/**
 * Per-pathogen generation times, loaded from the bundled
 * {@code generation_times/generation_times.yaml} (or a user-supplied override).
 * <p>
 * Replaces a single hardcoded 5-day default that was applied to every organism.
 * Because the phylodynamic Re fallback converts a growth rate to Re through a
 * serial interval whose mean is this value, a wrong generation time does not
 * merely add noise -- it compresses real differences toward Re = 1. Across the
 * project's 16-pathogen corpus that default produced Re in 1.0007-1.0262 for
 * everything from a DNA arbovirus to a gut bacterium.
 * <p>
 * Three states are distinguished, and only the first yields a number:
 * <ul>
 *   <li>{@code re_applicable: true} with a positive {@code T_days} -- usable.</li>
 *   <li>{@code re_applicable: false} -- the organism is acquired environmentally,
 *       via toxin, or from food rather than through a chain of host-to-host
 *       transmission, so no serial interval exists and Re should not be
 *       reported at all.</li>
 *   <li>{@code re_applicable: unknown} -- a stub awaiting a real value.</li>
 * </ul>
 * The latter two raise {@link MissingGenerationTimeException} rather than
 * silently substituting a default.
 */
public final class GenerationTimeTable {

    private static final String BUNDLED_RESOURCE = "/generation_times/generation_times.yaml";

    /** One pathogen's entry. {@code tDays} is null unless {@code reApplicable} is TRUE. */
    public record Entry(String pathogenId, String displayName, Double tDays, String confidence,
                        Applicability reApplicable, boolean vectorBorne, String note) {
    }

    public enum Applicability { TRUE, FALSE, UNKNOWN }

    private final Map<String, Entry> entries;

    private GenerationTimeTable(Map<String, Entry> entries) {
        this.entries = entries;
    }

    /** Loads the table bundled inside the jar. */
    public static GenerationTimeTable bundled() {
        try (InputStream in = GenerationTimeTable.class.getResourceAsStream(BUNDLED_RESOURCE)) {
            if (in == null) {
                throw new GviInputException("Bundled generation-time table not found on the classpath at " + BUNDLED_RESOURCE);
            }
            return parse(new ObjectMapper(new YAMLFactory()).readTree(in));
        } catch (IOException e) {
            throw new GviInputException("Could not read the bundled generation-time table: " + e.getMessage());
        }
    }

    /** Loads a user-supplied override file with the same schema. */
    public static GenerationTimeTable fromFile(Path path) {
        try (InputStream in = Files.newInputStream(path)) {
            return parse(new ObjectMapper(new YAMLFactory()).readTree(in));
        } catch (IOException e) {
            throw new GviInputException("Could not read generation-time table '" + path + "': " + e.getMessage());
        }
    }

    private static GenerationTimeTable parse(JsonNode root) {
        JsonNode table = root == null ? null : root.get("generation_times");
        if (table == null || !table.isObject()) {
            throw new GviInputException("Generation-time table has no 'generation_times:' mapping at its top level");
        }
        Map<String, Entry> parsed = new LinkedHashMap<>();
        for (Iterator<Map.Entry<String, JsonNode>> it = table.fields(); it.hasNext(); ) {
            Map.Entry<String, JsonNode> field = it.next();
            String key = field.getKey().toLowerCase(Locale.ROOT);
            JsonNode v = field.getValue();

            Applicability applicability = parseApplicability(v.path("re_applicable"), key);
            JsonNode t = v.get("T_days");
            Double tDays = (t == null || t.isNull()) ? null : t.asDouble();

            // Invariants -- a malformed table should fail at load, not produce a wrong Re later.
            if (applicability == Applicability.TRUE && (tDays == null || tDays <= 0)) {
                throw new GviInputException("Generation-time table entry '" + key
                        + "' is marked re_applicable: true but has no positive T_days");
            }
            if (applicability != Applicability.TRUE && tDays != null) {
                throw new GviInputException("Generation-time table entry '" + key + "' has T_days set but is marked re_applicable: "
                        + applicability.name().toLowerCase(Locale.ROOT) + " -- a non-applicable entry must not carry a value");
            }

            parsed.put(key, new Entry(key, text(v, "display_name", key), tDays, text(v, "confidence", "none"),
                    applicability, v.path("vector_borne").asBoolean(false), text(v, "T_note", "")));
        }
        return new GenerationTimeTable(parsed);
    }

    private static Applicability parseApplicability(JsonNode n, String key) {
        if (n == null || n.isNull()) {
            throw new GviInputException("Generation-time table entry '" + key + "' is missing 're_applicable'");
        }
        if (n.isBoolean()) return n.asBoolean() ? Applicability.TRUE : Applicability.FALSE;
        String s = n.asText("").trim().toLowerCase(Locale.ROOT);
        return switch (s) {
            case "true", "yes" -> Applicability.TRUE;
            case "false", "no" -> Applicability.FALSE;
            case "unknown" -> Applicability.UNKNOWN;
            default -> throw new GviInputException("Generation-time table entry '" + key
                    + "' has an unrecognized re_applicable value '" + n.asText() + "' (expected true, false or unknown)");
        };
    }

    private static String text(JsonNode n, String field, String fallback) {
        JsonNode v = n.get(field);
        return (v == null || v.isNull()) ? fallback : v.asText(fallback);
    }

    public java.util.Optional<Entry> find(String pathogenId) {
        if (pathogenId == null) return java.util.Optional.empty();
        return java.util.Optional.ofNullable(entries.get(pathogenId.toLowerCase(Locale.ROOT)));
    }

    /** Ids that currently yield a usable generation time. */
    public List<String> usableIds() {
        List<String> out = new ArrayList<>();
        entries.forEach((k, v) -> {
            if (v.reApplicable() == Applicability.TRUE) out.add(k);
        });
        return out;
    }

    /**
     * The generation time a run should use, given what the caller supplied.
     * <p>
     * Precedence: an explicit value always wins; otherwise a pathogen id is looked up here and
     * fails loudly rather than substituting a default; with neither, the caller's default is
     * carried forward.
     * <p>
     * This lives on the table rather than in a front end because it used to live in {@code GviCli},
     * which meant only the command line ever consulted the table. The web layer built its
     * {@code PipelineConfig} directly and silently ran every request on the 5-day default, while
     * reporting that the table had been used. Any front end that resolves a generation time must
     * reach the same answer, so the resolution belongs beside the data it reads.
     *
     * @param explicit          whether the caller supplied a generation time directly
     * @param explicitValueDays that value, used when {@code explicit}
     * @param pathogenId        table key, consulted only when not explicit
     * @param fallbackDays      used when neither an explicit value nor a pathogen id is given
     */
    public static double resolve(boolean explicit, double explicitValueDays,
                                 String pathogenId, double fallbackDays,
                                 GenerationTimeTable table) {
        if (explicit) return explicitValueDays;
        if (pathogenId == null || pathogenId.isBlank()) return fallbackDays;
        return table.requireGenerationTimeDays(pathogenId);
    }

    /** Every entry, ordered by id -- for a UI that must show what is and is not usable. */
    public List<Entry> allEntries() {
        List<Entry> out = new ArrayList<>(entries.values());
        out.sort(java.util.Comparator.comparing(Entry::pathogenId));
        return out;
    }

    public int size() {
        return entries.size();
    }

    /**
     * The generation time for {@code pathogenId}, or a {@link MissingGenerationTimeException}
     * explaining specifically why there isn't one -- never a silent default.
     */
    public double requireGenerationTimeDays(String pathogenId) {
        Entry e = find(pathogenId).orElseThrow(() -> new MissingGenerationTimeException(
                "No generation-time entry for --pathogen-id '" + pathogenId + "'. Add one to generation_times.yaml, "
                        + "or pass --generation-time-days explicitly. Ids with a usable value: " + String.join(", ", usableIds())));

        return switch (e.reApplicable()) {
            case TRUE -> e.tDays();
            case FALSE -> throw new MissingGenerationTimeException(
                    "Re is not applicable to '" + e.displayName() + "' (--pathogen-id " + e.pathogenId() + "): " + e.note()
                            + " Because there is no host-to-host transmission chain, there is no serial interval and no "
                            + "meaningful reproduction number -- run without the Re index rather than supplying a placeholder.");
            case UNKNOWN -> throw new MissingGenerationTimeException(
                    "Generation time for '" + e.displayName() + "' (--pathogen-id " + e.pathogenId()
                            + ") is still a stub in generation_times.yaml. Populate T_days and set re_applicable, or pass "
                            + "--generation-time-days explicitly for this run.");
        };
    }
}
