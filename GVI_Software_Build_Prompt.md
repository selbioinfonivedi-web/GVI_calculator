

## 0. Role & Mandate

You are acting as the technical lead of a 20-person multidisciplinary engineering
team building a production-grade, standalone desktop application in **Java**.
You must simulate and satisfy the concerns of every discipline on this team as
you design, implement, and test the software:

| # | Role | Count | Responsibility |
|---|------|-------|-----------------|
| 1 | Technical Lead / Architect | 1 | Overall architecture, module contracts, code review gate |
| 2 | Senior Java Engineers (core platform) | 4 | I/O engine, memory management, packaging, crash-resilience |
| 3 | Bioinformatics Algorithm Engineers | 4 | Implement the 8 index algorithms exactly per spec (Section 5) |
| 4 | Biostatisticians | 3 | Confidence intervals, hypothesis tests (PHI test, bootstrap), weight calibration, sensitivity analysis |
| 5 | Virology / Molecular Biology Domain Reviewers | 2 | Validate biological plausibility of outputs, reference value tables, gene-region annotations |
| 6 | UI/UX Engineer | 1 | JavaFX desktop UI + report design |
| 7 | QA / Test Engineers | 2 | Unit/integration/self-test framework, fuzz testing, benchmark datasets |
| 8 | DevOps / Release Engineer | 1 | Cross-platform packaging (jpackage), CI, reproducible builds |
| 9 | Technical Writer | 1 | User manual, algorithm reference doc, validation report |
| 10 | Project Coordinator | 1 | Milestones, acceptance sign-off |


## 1. Project Objective

Build a **standalone, offline, cross-platform desktop application** (Windows,
macOS, L**Java** that:inux) written in 

1. Ingests raw **genomic sequence data** (FASTA/multi-FASTA alignments, VCF
   variant calls, optional Newick phylogenetic trees) and **temporal
   metadata** (collection dates, and optionally case-incidence time series).
2. Computes **each of the 8 component indices independently** (user can run
   any subset) **and** the **composite Genomic Virulence Index (GVI)**.
3. Never crashes regardless of malformed input, missing fields, huge files,
   or low memory — it degrades gracefully with actionable error messages.
4. Ships with a **built-in self-test / self-diagnostic suite** that runs on
   first launch and on demand, validating the installation against bundled
   reference datasets with known expected outputs.
5. Runs identically via GUI (interactive) and CLI (headless/batch, scriptable
   for HPC or pipeline integration).

---

## 2. Non-Functional Requirements ("runs anywhere without crashing")

- **Language/Runtime**: Java 17 LTS (or newer LTS). No native/JNI
  dependencies — pure JVM bytecode only, so behavior is identical across OSes.
- **Packaging**: Use `jpackage` to produce self-contained native installers
  (`.msi`/`.exe` for Windows, `.dmg`/`.pkg` for macOS, `.deb`/`.rpm`/AppImage
  for Linux) that **bundle their own JRE** — the end user must not need to
  install Java separately, and must not be able to hit a "wrong JVM version"
  crash.
- **Global crash containment**:
  - Install a top-level `Thread.setDefaultUncaughtExceptionHandler` and a
    JavaFX `Application`-level exception handler that logs the full stack
    trace to a rotating log file and shows a user-friendly dialog — the
    application must **never** terminate silently or dump a raw stack trace
    to the end user.
  - No bare `System.exit()` calls on error paths. All I/O and parsing errors
    are caught, wrapped in typed exceptions (`GviInputException`,
    `GviComputationException`), and surfaced to the UI/CLI as actionable
    messages (e.g., "Line 4021 of input.fasta contains non-IUPAC character
    'X' — skipped" rather than a `StringIndexOutOfBoundsException`).
- **Memory safety for large genomic datasets**:
  - Stream-parse FASTA/VCF files (do not load multi-GB files fully into
    heap). Use buffered readers with bounded look-ahead.
  - Before any O(n²) operation (e.g., pairwise π, GD matrices), estimate
    memory requirement from sequence count × length, compare against
    `Runtime.getRuntime().maxMemory()`, and if it would exceed a safe
    threshold, automatically switch to a **chunked/streaming pairwise
    algorithm** or warn the user and offer subsampling — never let the JVM
    hit `OutOfMemoryError` unhandled.
  - Default heap (`-Xmx`) set conservatively in the launcher, but
    user-configurable via a settings panel.
- **Input robustness**: a schema-validation layer runs before any
  computation and rejects/flags malformed records individually (never aborts
  the whole batch for one bad record), with a summary report of skipped
  records.
- **Determinism**: identical input + identical config (weights, RNG seed for
  bootstrap) must always produce identical output, for reproducibility and
  for the self-test suite to work.
- **No telemetry / fully offline**: no network calls required for core
  computation (this is a scientific/epidemiological tool that may run in
  air-gapped public-health environments).
- **Test coverage gate**: ≥85% line coverage on all algorithm modules
  (enforced in CI via JaCoCo), required before merge.

---

## 3. Recommended Technology Stack

| Concern | Library | Why |
|---|---|---|
| Sequence I/O (FASTA/FASTQ) | BioJava (`biojava-core`, `biojava-alignment`) | Mature, pure-Java, IUPAC-aware |
| VCF parsing | HTSJDK | Standard for variant call files |
| Statistics (regression, bootstrap, distributions) | Apache Commons Math3 | Pure Java, no native deps |
| CLI parsing | picocli | Clean annotated CLI, works headless |
| Desktop UI | JavaFX | Cross-platform, bundles cleanly with jpackage |
| Charts | JavaFX Charts or JFreeChart | Diversity plots, GVI trend, dN/dS gene maps |
| JSON/CSV export | Jackson, Apache Commons CSV | Structured report output |
| Logging | SLF4J + Logback (rotating file appender) | Crash forensics without user-visible noise |
| Testing | JUnit 5, AssertJ, JaCoCo | Unit + coverage |
| Packaging | `jpackage` (JDK built-in) | Self-contained native installers |
| Build | Maven or Gradle (multi-module) | One module per index + core + UI + CLI |

Multi-module layout suggestion:
```
gvi-core/          # data models, I/O, validation, exceptions
gvi-algorithms/    # one package per index: mu, re, pi, mb, dnds, gd, ri, cai
gvi-composite/      # GVI weighting engine, β integration, calibration
gvi-selftest/       # bundled benchmark datasets + expected outputs + runner
gvi-cli/            # picocli entrypoint, headless batch mode
gvi-ui/             # JavaFX desktop app
gvi-report/         # PDF/CSV/JSON report generation
```

---

## 4. Input Data Specifications

### 4.1 Genomic input
- **Multi-FASTA alignment** (`.fasta`/`.fa`/`.aln`): aligned nucleotide
  sequences, one reference + N query sequences. Must handle gaps (`-`),
  ambiguity codes (IUPAC: N, R, Y, W, S, K, M, B, D, H, V), mixed case.
- **VCF** (`.vcf`/`.vcf.gz`): variant calls relative to a reference, with
  per-site depth (`DP`) for coverage-threshold filtering (default ≥10×,
  configurable).
- **Reference genome** (`.fasta`) + optional **gene/CDS annotation** (GFF3)
  — required for codon-level metrics (dN/dS, CAI) so ORFs can be extracted
  in the correct reading frame.
- **Optional Newick tree** (`.nwk`) — enables tree-informed distance and
  (optionally) coalescent-based Re estimation.
- **Host reference codon usage table** (bundled defaults: human, and a
  small set of common host organisms; user-loadable custom table) — required
  for CAI.

### 4.2 Temporal input
- **Metadata CSV**: one row per sequence, columns: `sequence_id`,
  `collection_date` (ISO-8601), optional `location`, optional
  `host`. Required to align genomic records with time for μ, MB-rate, and
  temporal signal regression.
- **Case-incidence time series CSV** (optional, for Re): `date`,
  `new_cases`. If absent, Re falls back to the phylodynamic/root-to-tip
  estimator (Section 5.2).

### 4.3 Validation rules (enforced before computation)
- Sequence IDs in FASTA must have a matching row in metadata CSV (warn +
  exclude if missing, never crash).
- Dates must parse as valid calendar dates; reject future dates beyond
  today with a warning.
- Alignment: all sequences compared pairwise must be equal length (report
  and exclude any that are not, rather than throwing an unhandled exception).

---

## 5. Core Algorithms — Implement Exactly As Specified

Every index must be computable **standalone** (given only the inputs it
needs) and **as part of the composite GVI**. Each module returns a typed
result object containing the raw value, the 0–1 normalized value used in the
composite, the reference-category classification (from the tables below),
and (where applicable) a confidence interval.

### 5.1 Index 1 — Evolutionary Rate (μ)

**Formula**:
```
μ = (total substitutions observed) / (number of sequenced sites × evolutionary time interval in years)
```

**Implementation**: primary estimator is **root-to-tip regression** (linear
regression of genetic distance from each sequence to a reference/root
against that sequence's sampling date — the TempEst/BEAST "clock signal"
approach), since raw substitution counts alone conflate divergence time.
1. For each sequence, compute genetic distance to the earliest/reference
   sequence using Jukes-Cantor (Section 5.6).
2. Regress distance (y) against decimal collection date (x) using
   Apache Commons Math `SimpleRegression`.
3. μ = slope of the regression line (substitutions/site/year).
4. Report R² of the regression as a "temporal signal strength" diagnostic —
   if R² < 0.3, flag the μ estimate as unreliable in the output (do not
   silently present a bad estimate as precise).

**Reference classification table** (bundle as config, not hardcoded magic
numbers):
| μ range | Category |
|---|---|
| <1e-4 | High-fidelity (measles/polio-like) |
| 1e-4–1e-3 | Moderate RNA (dengue/WNV/Zika-like) |
| 1e-3–1e-2 | Low-fidelity RNA (influenza/SARS-CoV-2/HIV-like) |
| >1e-2 | Ultra-high error (VSV-like) |

**β integration**: `β_adjusted = β₀ × (1 + (μ/μ_ref) × weight_μ)`

### 5.2 Index 2 — Effective Reproduction Number (Re)

Implement **two estimators**, selected automatically by data availability
(document which was used in the output):

**(a) Incidence-based (Cori et al. 2013, "EpiEstim" method)** — used when a
case-incidence time series is provided:
```
Re(t) = Σ_{s=t-τ..t} I(s)  /  Σ_{s=t-τ..t} Σ_{u=1..s} I(s-u) · w(u)
```
where `I` is daily incidence, `w(u)` is the serial-interval distribution
(default: bundled pathogen-agnostic Gamma(mean=5, sd=2) days, user-
overridable), and τ is a sliding window (default 7 days). Report the
posterior mean and 95% credible interval using the Cori Gamma-Poisson
conjugate posterior (closed form — implement directly, no MCMC needed).

**(b) Phylodynamic/root-to-tip fallback** — used when only genomic +
temporal data are available (no incidence series): estimate a **coalescent
exponential-growth rate** from the branch-length/tip-date regression
residual structure (Pybus & Rambaut-style approach: lineage-through-time
slope), then convert growth rate `r` to Re via
`Re = 1 + r × generation_time` (generation_time configurable, default from
Section 5.1's μ-linked pathogen class or user input).

**Reference table**: bundle the doc's Re category table (`<0.8` Controlled …
`>2.0` Pandemic) as config, with associated doubling-time formula
`doubling_time = ln(2) / r`.

**β integration**: `β ∝ Re` directly — `β_factor = Re / Re_baseline`.

### 5.3 Index 3 — Nucleotide Diversity (π)

**Formula**:
```
π = (Σ pairwise differences) / (number of pairwise comparisons × alignment length)
```
Implementation:
- For n aligned sequences of length L, compute Hamming distance for every
  pair (i,j), i<j, excluding gap/gap and N/any comparisons from both the
  numerator and the effective L for that pair (pairwise deletion).
- π = (Σ over all pairs of [differences / valid-sites-in-pair]) / C(n,2).
- **Complexity guard** (per Section 2): if n > ~2000, use random pair
  subsampling (default 100,000 pairs, user-configurable) with a bootstrap
  CI instead of full O(n²), and clearly label the result as an estimate.
- Report 95% CI via bootstrap resampling of sequences (biostatistician-
  owned module), default 1,000 replicates, seeded RNG for reproducibility.

**Reference table**: bundle doc's π thresholds (<0.0005 acute outbreak …
>0.02 hyperendemic).

### 5.4 Index 4 — Mutation Burden (MB)

**Formula**: count of positions where the query differs from reference
(SNP = 1, indel of any length = 1 event), **excluding sites below the
coverage threshold** (default ≥10×, from VCF `DP` field).

Implementation:
- Prefer VCF input (has DP for proper coverage filtering); fall back to
  alignment-based diffing against the reference row if VCF unavailable
  (coverage filter then unavailable — flag this in the report).
- Group contiguous deleted/inserted bases into a single indel event (do not
  count each base of a multi-base indel separately).
- Derived metric: `outbreak_age_weeks ≈ MB / μ` (only computed when a μ
  estimate is available for the same lineage), surfaced as a secondary
  output, not the primary MB value.

**Reference table**: bundle doc's MB thresholds (0–5 within-household … 40–100 pandemic variant).

### 5.5 Index 5 — dN/dS Ratio (ω)

**Method**: implement the **Nei–Gojobori (1986) counting method** (the
standard, tractable, non-ML approach appropriate for a deterministic,
crash-proof standalone tool):
1. Translate reference and query CDS in the correct frame (from GFF3
   annotation) into codons.
2. For each codon in the reference, compute the number of synonymous (S)
   and nonsynonymous (N) sites using the standard degeneracy-based
   fractional-site method (average over the 3 possible single-nucleotide
   paths at each position using the standard genetic code table — bundle
   the codon table as a resource, not inline literals).
3. For each reference/query codon pair, count observed synonymous (Sd) and
   nonsynonymous (Nd) differences, correcting for multiple substitutions
   per codon by enumerating all possible mutational pathways between the two
   codons and averaging Sd/Nd across pathways (standard Nei-Gojobori
   correction).
4. Sum S, N, Sd, Nd across all codons in the gene/region.
5. `pS = Sd/S`, `pN = Nd/N`.
6. Apply **Jukes-Cantor correction** to pS and pN individually:
   `dS = -3/4 · ln(1 - 4/3 · pS)`, `dN = -3/4 · ln(1 - 4/3 · pN)`.
7. `dN/dS = dN / dS`.
8. Compute per-gene/per-region if annotation defines multiple ORFs (e.g.,
   report separately for Spike, RdRp, N, etc. when annotation available),
   since interpretation is gene-specific (Section 5.5 reference table).

**Reference table**: bundle doc's gene-specific dN/dS ranges (RNA
polymerase <0.3 … Envelope 1.5–2.5) as annotatable per-region config, matched
by gene name in the GFF3.

**β integration**: dN/dS > 1 in transmissibility-linked genes directly
raises predicted Re/β contribution — implement as a configurable per-gene
multiplier in the composite module (Section 5.9).

### 5.6 Index 6 — Genetic Distance (GD)

Implemented all three, selectable by user/config:
- **Hamming**: `p` = proportion of differing (non-gap, non-N) sites.
- **Jukes-Cantor (JC69)**: `d_JC = -3/4 · ln(1 - 4/3 · p)`
- **Kimura 2-Parameter (K80)**: distinguish transitions (A↔G, C↔T) from
  transversions (A↔C, A↔T, G↔C, G↔T):
  `d_K80 = -1/2·ln(1 - 2P - Q) - 1/4·ln(1 - 2Q)`
  where P = transition proportion, Q = transversion proportion.
- Guard against the log-argument going ≤0 (saturation at high divergence) —
  return "distance saturated / unreliable, use raw p-distance" rather than
  `NaN`/crash on `ln` of a non-positive number.

**Reference table**: bundle doc's JC-distance time-since-divergence
thresholds (0–0.001 <1wk … >0.05 different strain/serotype).

### 5.7 Index 7 — Recombination Index (RI)

**Method**: implement the **PHI test (Pairwise Homoplasy Index, Bruen,
Bryant & Poss 2006)** as the statistical detection method (well-defined,
implementable without MCMC, gives a p-value a biostatistician can defend):
1. Slide a window across the alignment; compute the incompatibility (per-
   site homoplasy) matrix among informative sites.
2. Compute the observed mean pairwise homoplasy statistic (PHI).
3. Generate a null distribution via permutation of site order (default
   1,000 permutations, seeded RNG) and compute empirical p-value.
4. `RI = proportion of statistically significant recombinant windows` (or,
   operationally, `RI = detected breakpoints / genome length` when
   breakpoint scanning is enabled) — implement both the doc's two
   operational definitions (`recombinant isolates / total isolates` when
   isolate-level classification is available, and `breakpoints per genome`
   from the PHI scan) and let the user pick which is reported.

**Reference table**: bundle doc's RI thresholds (<0.02 single-lineage …
>0.30 reassortment-prone).

### 5.8 Index 8 — Codon Adaptation Index (CAI) & GC Content Deviation

**CAI** (Sharp & Li 1987):
1. Build/load a reference host codon-usage table (relative synonymous codon
   usage, RSCU) — bundle defaults, allow custom upload.
2. For each codon, `w(codon) = RSCU(codon) / RSCU(max synonymous codon)`.
3. For a gene of L codons (excluding stop codons and single-codon amino
   acids per convention), `CAI = exp( (1/L) · Σ ln w(codon_i) )` — i.e., the
   geometric mean of the w-values, computed in log-space for numerical
   stability.

**GC Content Deviation**:
```
GC%_observed = (count G + count C) / total non-ambiguous bases × 100
GC_Deviation = |GC%_observed - GC%_reference|
```

**Reference tables**: bundle doc's CAI bands (>0.80 optimized … <0.50
unadapted) and GC deviation bands (<1% native, 1–3% minor drift, >3%
anomaly/HGT).

### 5.9 Composite GVI

```
GVI(t) = w1·μ_norm + w2·Re_norm + w3·π_norm + w4·MB_norm
       + w5·(dN/dS)_norm + w6·GD_norm + w7·CAI_norm + w8·RI_norm
```
- Each raw index value is **min-max normalized to [0,1]** against a
  configurable reference range (defaults = the reference tables in Section
  5.1–5.8; user-overridable), before weighting.
- Default weights (bundled as editable config, Σw=1, midpoints of doc
  ranges):
  `μ=0.135, Re=0.30, π=0.10, MB=0.215, dN/dS=0.15, GD=0.10, CAI+GC=0.055, RI=0.035`
  (Note: CAI+GC counted as one combined weighted sub-block per the source
  spec — implement as `w7a·CAI_norm + w7b·GC_norm` summing to the CAI+GC
  weight, sub-split configurable.)
- `β_genomic(t) = β₀ × [1 + GVI(t) × scale_factor]` (scale_factor
  user-configurable, default 2.0), exposed as an optional output for users
  feeding this into an external SEIR/epidemic model — **not** itself a full
  epidemic simulator (out of scope; GVI outputs a value/time-series a
  downstream model consumes).
- **Missing-index handling**: if a user has data for only some indices
  (e.g., no incidence series so Re unavailable), the composite must
  **re-normalize the remaining weights to sum to 1** and clearly flag in
  the report which indices were excluded and why — never silently treat a
  missing index as zero.
- **Weight calibration module** (biostatistician-owned, optional/advanced):
  fit weights against an observed epidemic curve via least-squares or
  Approximate Bayesian Computation (ABC) using Apache Commons Math
  optimizers; hold-out validation (80/20 split), report RMSE, target <15%
  of mean per the source spec.
- **Sensitivity analysis module**: vary each weight ±20% one at a time,
  recompute GVI, report resulting spread — surfaced as a tornado chart in
  the UI report.

---

## 6. Independent + Combined Execution Model

- Each of the 8 algorithm modules implements a common
  `GviIndexCalculator<TInput, TResult>` interface with `compute()`,
  `validateInputs()`, and `getRequiredInputTypes()`.
- The CLI/UI lets the user select **any subset** of indices to run
  standalone (each producing its own report section) and/or the full
  composite GVI (auto-running only the indices for which input was
  supplied, per the missing-weight handling in 5.9).
- No module may depend on another module's implementation directly — only
  on the shared `gvi-core` data model — so any index can be tested, run, and
  reported in complete isolation.

---

## 7. Output & Reporting

- **Per-index report**: raw value, normalized value, reference-category
  classification, confidence interval (where applicable), diagnostic flags
  (e.g., low temporal signal, saturated distance, insufficient data).
- **Composite GVI report**: weighted breakdown (stacked bar of each index's
  contribution), GVI(t) time series if multiple timepoints supplied,
  β_adjusted / β_genomic output, sensitivity tornado chart.
- **Export formats**: CSV (machine-readable, one row per sequence/sample),
  JSON (full structured result incl. metadata/provenance), PDF (human-
  readable report with charts, for the UI users).
- **Provenance**: every report records input file hashes, software version,
  weight/config values used, and timestamp — for scientific reproducibility.

---

## 8. Testing & Validation Strategy

### 8.1 Unit tests (QA + algorithm engineers)
- Each formula tested against **hand-calculated or literature-published
  reference values** (e.g., known π from a published Dengue dataset, known
  dN/dS from a PAML tutorial dataset, known JC/K80 distances from textbook
  examples). Cite the source dataset/paper in the test class Javadoc.
- Edge cases per module: empty input, single sequence, all-identical
  sequences (π=0, GD=0), fully saturated divergence (JC log domain guard),
  sequences with only Ns, non-multiple-of-3 CDS length, missing GFF3
  annotation, VCF with no DP field, incidence series with zero cases,
  single-timepoint data (μ/Re temporal regression undefined — must return a
  clear "insufficient temporal range" result, not crash).

### 8.2 Integration tests
- Full pipeline: raw FASTA+metadata in → composite GVI report out, against
  2–3 bundled synthetic "known-answer" datasets constructed by the
  bioinformatics team with hand-verified expected outputs for every index.

### 8.3 Built-in Self-Test / Self-Diagnostic Mode
- Bundle the benchmark datasets from 8.2 as resources inside the
  application (not external downloads — must work fully offline/air-gapped).
- `--self-test` CLI flag and a "Run Diagnostics" UI menu item execute the
  full pipeline against bundled datasets and compare outputs to
  stored expected values within a defined numerical tolerance (e.g.,
  ±1e-6 for deterministic formulas, ±5% for bootstrap/permutation-based
  values given the seeded RNG should still make these exact — assert exact
  equality where the RNG is seeded).
- Reports PASS/FAIL per index with diffs, plus environment info (OS, JVM
  version, available memory) — this is the mechanism that lets a user
  confirm the software works correctly on *this specific machine* before
  trusting it on real data.
- Also run automatically, silently, once on first launch after install;
  surface a non-blocking warning banner if any check fails.

### 8.4 Fuzz / robustness testing
- QA generates a corpus of deliberately malformed inputs (truncated files,
  binary garbage renamed to `.fasta`, mismatched lengths, huge files,
  Unicode in headers, missing required CSV columns) and asserts the
  application **never crashes** — it either processes with warnings or
  rejects with a clear typed error, always exits/returns cleanly.

### 8.5 Statistical validation (biostatisticians)
- Verify bootstrap CIs achieve nominal coverage on simulated data with
  known true π/dN/dS.
- Verify PHI test false-positive rate on simulated non-recombinant
  alignments matches the nominal significance level.
- Verify Cori Re estimator reproduces published EpiEstim example outputs
  (cross-check against the R `EpiEstim` package's documented example
  dataset, since Cori et al.'s method is the same reference implementation
  most public-health agencies already trust).

### 8.6 Domain review (virologists/life scientists)
- Sign-off checklist: do reference-category boundaries and biological
  interpretations in the UI/report match the source definitions (Section
  5) verbatim; are gene-specific dN/dS annotations correctly mapped from
  GFF3 gene names to the reference table's gene categories.

---

## 9. Documentation Deliverables

1. **User Manual** — install, load data, run individual/composite indices,
   read reports.
2. **Algorithm Reference** — every formula in Section 5 with citations,
   reproduced from this spec, kept in sync with code.
3. **Validation Report** — results of Section 8's benchmark/self-test suite,
   signed off by biostatistics + domain review roles.
4. **Developer Guide** — module architecture, how to add a 9th index later
   without touching existing modules (interface-driven extensibility).

---

## 10. Development Phases

1. **Phase 1 — Foundation**: `gvi-core` data models, streaming I/O,
   validation layer, exception hierarchy, logging, crash-handler scaffolding.
2. **Phase 2 — Algorithms**: implement all 8 index modules in
   `gvi-algorithms`, each with unit tests against literature benchmarks
   (Section 8.1), independently runnable via a temporary test CLI.
3. **Phase 3 — Composite Engine**: `gvi-composite` weighting/normalization,
   missing-data re-normalization, calibration, sensitivity analysis.
4. **Phase 4 — Self-Test Framework**: `gvi-selftest` bundled datasets +
   runner, wired into both CLI and first-launch UI check.
5. **Phase 5 — CLI**: `gvi-cli` full headless batch mode with picocli,
   scriptable exit codes for pipeline integration.
6. **Phase 6 — UI**: `gvi-ui` JavaFX desktop app, charts, report viewer.
7. **Phase 7 — Reporting**: `gvi-report` PDF/CSV/JSON generation with
   provenance metadata.
8. **Phase 8 — Packaging & Release**: `jpackage` installers for all 3 OSes,
   fuzz-testing pass, final domain sign-off, documentation freeze.

---

## 11. Acceptance Criteria

- [ ] All 8 indices computable standalone from minimal required inputs only.
- [ ] Composite GVI computable from any combination of available indices
      with correct weight re-normalization for missing ones.
- [ ] Zero unhandled exceptions across the full fuzz-testing corpus.
- [ ] Self-test suite passes on a clean install with no network access, on
      Windows, macOS, and Linux.
- [ ] ≥85% test coverage on `gvi-algorithms` and `gvi-composite`.
- [ ] All formulas in the shipped Algorithm Reference doc match the
      implementation exactly (spot-checked by biostatistics + domain
      reviewers against this spec).
- [ ] Application launches and runs the self-test successfully with default
      JVM heap settings on a machine with as little as 2GB available RAM
      (using the streaming/chunking guards from Section 2), and reports a
      clear "reduce dataset size" message rather than crashing if a
      genuinely too-large dataset is loaded.
