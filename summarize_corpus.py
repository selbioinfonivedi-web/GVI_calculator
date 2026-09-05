#!/usr/bin/env python3
"""Builds the cross-pathogen comparison table from a directory of GVI result JSONs.

The previous summary put every pathogen's GVI in one column with nothing to say how much evidence
backed each number -- so enterotoxaemia's 0.053 (built from 2 of 9 indices, under 6% of the weighting
scheme) sat beside haemorrhagic septicaemia's 0.680 (6 of 9) as though the two were the same kind of
measurement. They are not comparable, and a reader sorting that column would rank them anyway.

This version carries coverage and a COMPARABLE flag next to every score, and refuses to rank the
non-comparable ones.
"""
import csv
import json
import pathlib
import sys

KEYS = ["MU", "RE", "PI", "MB", "DNDS", "GD", "RI", "CAI", "GC"]


def main(outdir: str) -> int:
    out = pathlib.Path(outdir)
    rows = []
    for path in sorted(out.glob("*.json")):
        try:
            data = json.loads(path.read_text())
        except json.JSONDecodeError:
            continue
        gvi = data.get("dataset_gvi")
        if gvi is None:
            continue

        scored = {c["key"]: c for c in gvi.get("components", [])}
        indices = {}
        indices.update(data.get("population_indices") or {})
        indices.update(data.get("dataset_indices") or {})

        coverage = gvi.get("effectiveWeightSum", 0.0)
        row = {
            "pathogen": path.stem,
            "GVI": round(gvi.get("gvi", 0.0), 6),
            "indices_scored": len(scored),
            "coverage_pct": round(coverage * 100, 1),
            # 0.60 mirrors GviResult.MIN_COMPARABLE_COVERAGE. A score below it summarises a much smaller
            # body of evidence than a full one and must not be ranked against it.
            "comparable": "yes" if coverage >= 0.60 else "NO",
        }
        for key in KEYS:
            entry = indices.get(key)
            row[key] = "" if entry is None else round(entry.get("primaryValue", 0.0), 8)
            row[key + "_scored"] = "" if entry is None else ("yes" if key in scored else "no")
        rows.append(row)

    if not rows:
        print("No result JSONs found in", out, file=sys.stderr)
        return 1

    # Comparable scores first and ranked; non-comparable ones listed after, alphabetically, so the file
    # never implies an ordering between the two groups.
    rows.sort(key=lambda r: (r["comparable"] != "yes", -r["GVI"] if r["comparable"] == "yes" else 0,
                             r["pathogen"]))

    fields = ["pathogen", "GVI", "comparable", "coverage_pct", "indices_scored"]
    for key in KEYS:
        fields += [key, key + "_scored"]

    target = out / "ALL_PATHOGENS_SUMMARY.csv"
    with target.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    comparable = [r for r in rows if r["comparable"] == "yes"]
    print(f"\n{len(rows)} datasets -- {len(comparable)} comparable, {len(rows) - len(comparable)} not")
    print(f"{'pathogen':28} {'GVI':>8}  {'cov':>5}  {'idx':>4}  comparable")
    for r in rows:
        print(f"{r['pathogen']:28} {r['GVI']:8.4f}  {r['coverage_pct']:4.0f}%  {r['indices_scored']:4}  {r['comparable']}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else "."))
