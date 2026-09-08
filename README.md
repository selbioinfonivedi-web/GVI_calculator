# GVI Calculator — Genomic Virulence Index

A standalone, fully offline implementation of the eight-index Genomic Virulence Index
framework. Takes an aligned multi-FASTA and produces a single scalar, GVI(t) ∈ [0,1],
by computing eight named genomic indices, scaling each to a reference range, and
combining them under a weighted sum.

It runs entirely in one JVM process. No MAFFT, no IQ-TREE, no LSD2, no network.

```
GVI(t) = w₁·μ + w₂·Re + w₃·π + w₄·MB + w₅·dN/dS + w₆·GD + w₇·CAI + w₈·RI     (Σw = 1)
```

---

## Quick start

Requires **JDK 17**. Build with Maven:

```bash
cd gvi-calculator-java
mvn install
```

### Web interface (the supported surface)

```bash
java -jar gvi-calculator-java/gvi-web/target/gvi-calculator-web.jar
# → http://127.0.0.1:8080
```

Binds loopback only. It accepts uploads and runs analyses with no authentication,
so it is a local analyst's tool rather than a service; `--host` allows a deliberate
deployment behind a reverse proxy and warns when used.

### Command line

```bash
java -jar gvi-calculator-java/gvi-cli/target/gvi-calculator.jar \
    --fasta aligned.fasta \
    --metadata metadata.csv \
    --gff annotation.gff3 \
    --organism-class virus --genome-type rna \
    --pathogen-id fmd \
    --json result.json
```

`--self-test` runs 13 diagnostic checks inside the delivered binary.

### Desktop application

Parked. The web tool is the supported interface; the JavaFX module is excluded from
the default build but still compiles:

```bash
mvn -Pdesktop install
```

---

## Repository layout

| Path | Contents |
|---|---|
| `gvi-calculator-java/` | The Maven reactor — six modules built by default, plus the parked desktop module |
| `pathogen_data/` | The 16-dataset corpus. Three files per dataset: `aligned.fasta`, `metadata.csv`, a GFF3 |
| `gvi_results_final/` | Current corpus results — **the baseline the regression diffs against** |
| `run_corpus.sh` | Runs every dataset end to end and rebuilds the summary |
| `summarize_corpus.py` | Builds `ALL_PATHOGENS_SUMMARY.csv` from a results directory |
| `Genomic_Indices_Detailed_Definitions.docx` | The source specification |
| `.github/workflows/ci.yml` | Build, self-test, corpus regression, release on tag |
| `gui_concepts/` | Interface design explorations |

Every push runs the four CI jobs; the corpus job diffs all 16 datasets against
`gvi_results_final/` and fails on any drift.

### Modules

```
gvi-core          model, readers, genetic code, trimming, exceptions
   └── gvi-algorithms      every estimator, reference table, phylogenetic primitive
          └── gvi-composite     weighting, normalisation, renormalisation, sensitivity
                 ├── gvi-selftest    the 13 checks shipped in the binary
                 └── gvi-cli         picocli entry point + the pipeline orchestrator
                        ├── gvi-web        HTTP server, JSON API, browser client
                        └── gvi-ui         JavaFX desktop (parked behind -Pdesktop)
```

Every front end calls the same `GviPipeline.run(PipelineConfig)` through the canonical
constructor. This is enforced rather than assumed: the desktop UI once diverged from the
CLI by calling a back-compatible constructor that defaulted organism class and genome
type, so identical input scored differently depending on which front end ran it.

---

## Testing

```bash
cd gvi-calculator-java
mvn test                      # 322 tests
java -jar gvi-cli/target/gvi-calculator.jar --self-test
bash ../run_corpus.sh /tmp/corpus_check    # all 16 datasets end to end
```

Three layers, because they catch different things:

- **Unit tests** prove code does what it was written to do.
- **Ground-truth recovery tests** simulate data from known parameters and check what
  comes back. μ recovers 0.002 to within 0.6%; dN/dS recovers 0.3 and 2.5 to within 6%;
  the birth–death propagator agrees with direct Monte Carlo to four decimal places.
- **The corpus run** catches what neither can. Two real defects reached stored results
  while every test stayed green: CAI was silently skipped on all seven datasets that
  could compute it, and stored files drifted behind the source on RI and μ. Both were
  obvious in a corpus diff.

`gvi_results_final/` is committed for exactly that reason. Regenerating it casually
defeats the purpose — a baseline change should carry the diff and a statement of which
scientific decision moved the number.

---

## Reading a result

Two fields decide whether a score means anything:

- **`effectiveWeightSum`** — the share of the weighting scheme that had usable data.
- **`comparable`** — true when that share is at least 0.60.

The composite renormalises whichever indices survive quality gating, so a score built
from two indices lands in [0,1] and looks exactly like one built from nine. In this
corpus, enterotoxaemia scored 0.053 from two indices at 5.7% coverage — beside
haemorrhagic septicaemia's 0.680 that reads as "much lower risk", when in fact the data
was too broken to score. **A near-zero GVI from near-zero coverage is indistinguishable,
by value alone, from a confident finding of low virulence.**

Nine of the sixteen corpus datasets currently clear the floor.

---

## Known limitations

- **The weights are asserted, not fitted.** They are the midpoints of ranges in the
  specification, rescaled from a 1.09 sum to 1.0. Nothing has been calibrated against an
  observed outcome. Report the ±20% sensitivity sweep alongside any score.
- **Re depends heavily on the generation time.** δ = 365.25 / generation time scales the
  whole birth–death process, so a value wrong by a factor of *k* moves Re by roughly the
  same factor — and Re carries the largest weight (0.2752). Supply `--pathogen-id` or
  `--generation-time-days`; the bundled table currently has a usable value for only 9 of
  its 27 entries.
- **Case-incidence data is the accurate path for Re.** The Cori estimator recovers a
  known Re to within 0.02; tree-shape inference is far weaker. No corpus dataset
  currently supplies incidence data.
- **π and GD saturate their normalisation ceilings** on inter-serotype alignments, so
  they contribute the maximum the scheme allows and carry no discriminating information.
- **The indices are not orthogonal.** μ, π, MB, GD and GC intercorrelate at r ≥ 0.89
  while carrying 52.98% of the scheme.

---

## Licence

This project is released under the **MIT Licence** — see [LICENSE](LICENSE).
Copyright © 2026 Pooja B.N.

Two things bundled here carry their own terms and are **not** covered by the MIT
licence above:

| Component | Licence | Where |
|---|---|---|
| IBM Plex typefaces | SIL Open Font License 1.1 | `gvi-calculator-java/gvi-ui/src/main/resources/fonts/LICENSE-IBM-Plex.txt` |
| Codon usage tables | Kazusa Codon Usage Database (Nakamura, Gojobori & Ikemura 2000) — cite on use | `gvi-calculator-java/gvi-algorithms/src/main/resources/codon_usage/` |

Runtime dependencies (Apache Commons, Jackson, picocli, SLF4J, Logback, JavaFX,
JUnit, AssertJ) are Apache-2.0, MIT, EPL or GPL+CE as declared in their own
distributions; none are redistributed in source form here.

### Citing

If this tool contributes to published work, cite the specification it implements
alongside the methods it uses — Cori et al. 2013 (Re), Stadler 2010 (birth–death
sampling), Nei & Gojobori (dN/dS), Bruen, Bryant & Poss 2006 (PHI), Sharp & Li 1987
(CAI), Saitou & Nei 1987 (neighbour-joining), Yang 1994 (discrete gamma),
Felsenstein 1985 (bootstrap), To et al. 2016 (least-squares dating).
