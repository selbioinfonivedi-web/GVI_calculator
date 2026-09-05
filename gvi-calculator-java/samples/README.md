# Example input files

Ready-to-use examples showing the exact format each optional input expects.
`example_metadata.csv` uses the same sequence IDs as `samples.fasta`
(`reference`, `query1`, `query2`, `query3`), so all four files work together
in one real run:

```bash
java -jar gvi-cli/target/gvi-calculator.jar \
  --fasta samples/samples.fasta \
  --metadata samples/example_metadata.csv \
  --incidence samples/example_incidence.csv \
  --weights samples/example_weights.json \
  --indices mu,re
```

- **`samples.fasta`** — 4 short toy sequences (1 reference + 3 queries), for exercising the pipeline without needing real data.
- **`example_metadata.csv`** — `sequence_id,collection_date,location,host`. Collection dates unlock the temporal indices (mu, Re) when a FASTA's own headers don't already encode a date.
- **`example_incidence.csv`** — `date,new_cases`. A synthetic 30-day case-count series (a real outbreak-shaped curve: growth, then plateau) for the Cori et al. Re estimator. Real incidence data comes from public health surveillance (WHO situation reports, national health ministry dashboards, ProMED, HealthMap), not this tool.
- **`example_weights.json`** — a partial custom composite-weights override (`{"Re": 0.45, "MB": 0.25, "dN/dS": 0.10}`). Any of the 9 index labels not mentioned (`mu, Re, pi, MB, dN/dS, GD, CAI, GC_Deviation, RI`) keeps its spec-default weight, and the full set is renormalized to sum to 1.0 -- see `--help` on `--weights` for the complete rule.

All three were verified against the real pipeline (`mvn test`; also run manually against `samples.fasta` while writing them) -- these are working examples, not just illustrative snippets.
