#!/usr/bin/env bash
# Re-runs every dataset in the corpus through the current build and rebuilds the comparison summary.
#
# There was no script for this: the stored result.json files were produced ad-hoc and have repeatedly
# drifted behind the source (RI and mu both changed under them). Anything quoted as "current behaviour"
# should come from a fresh run of this script, not from a file on disk.
#
# --organism-class and --genome-type are set per dataset and are NOT optional metadata. Organism class
# gates whether host-relative CAI is meaningful; genome type sets which ceiling mu is normalised and
# judged against. Getting either wrong silently changes the score.
set -uo pipefail

export JAVA_HOME="${JAVA_HOME:-$HOME/opt/jdk-17.0.20.1+1}"
export PATH="$JAVA_HOME/bin:$PATH"

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
JAR="$ROOT/gvi-calculator-java/gvi-cli/target/gvi-calculator.jar"
OUT="${1:-$ROOT/gvi_results_$(date +%F)}"

[ -f "$JAR" ] || { echo "Build first: (cd gvi-calculator-java && mvn install -DskipTests)"; exit 1; }
mkdir -p "$OUT"

# name|dir|organism-class|genome-type|gff(optional)|host-codon-table(viruses only)
DATASETS=(
  "african_swine_fever|pathogen_data/viruses/african_swine_fever|virus|dna|annotation.gff3|pig"
  "classical_swine_fever|pathogen_data/viruses/classical_swine_fever|virus|rna|annotation.gff3|pig"
  "foot_and_mouth_disease|pathogen_data/viruses/foot_and_mouth_disease|virus|rna|annotation.gff3|cattle"
  "blue_tongue|pathogen_data/viruses/blue_tongue|virus|rna|annotation.gff3|sheep"
  "lumpy_skin_disease|pathogen_data/viruses/lsd|virus|dna|annotation.gff3|cattle"
  "sheep_goat_pox|pathogen_data/viruses/sheep_goat_pox|virus|dna|annotation.gff3|sheep"
  "ppr|pathogen_data/viruses/ppr|virus|rna|annotation.gff3|goat"
  "anthrax|pathogen_data/bacteria/anthrax|bacterium|dna|genes.gff3|"
  "haemorrhagic_septicaemia|pathogen_data/bacteria/haemorrhagic_septicaemia|bacterium|dna|genes.gff3|"
  "enterotoxaemia|pathogen_data/bacteria/enterotoxaemia|bacterium|dna|genes.gff3|"
  "blackleg|pathogen_data/bacteria/blackleg|bacterium|dna|genes.gff3|"
  "black_quarter|pathogen_data/bacteria/black_quarter|bacterium|dna|genes.gff3|"
  "babesiosis|pathogen_data/parasites/babesiosis|parasite|dna||"
  "theileriosis|pathogen_data/parasites/theileriosis|parasite|dna||"
  "trypanosomiasis|pathogen_data/parasites/trypanosomiasis|parasite|dna||"
  "fascioliasis|pathogen_data/parasites/fascioliasis|parasite|dna|cox1.gff3|"
)

echo "Writing results to $OUT"
for entry in "${DATASETS[@]}"; do
  IFS='|' read -r name dir cls genome gff host <<< "$entry"
  fasta="$ROOT/$dir/aligned.fasta"
  [ -f "$fasta" ] || { echo "  SKIP $name (no aligned.fasta)"; continue; }

  args=(--fasta "$fasta" --organism-class "$cls" --genome-type "$genome" --trim-to-covered
        --json "$OUT/$name.json" --csv "$OUT/$name.csv")
  [ -f "$ROOT/$dir/metadata.csv" ] && args+=(--metadata "$ROOT/$dir/metadata.csv")
  [ -n "$gff" ] && [ -f "$ROOT/$dir/$gff" ] && args+=(--gff "$ROOT/$dir/$gff")
  # CAI needs a host codon table. Without one the index is skipped and the remaining
  # weights are renormalized, so the run still prints a GVI -- just a different one.
  # The host is taken from each dataset's own metadata host column.
  [ -n "$host" ] && args+=(--codon-usage-species "$host")

  printf '  %-26s ' "$name"
  if java -jar "$JAR" "${args[@]}" > "$OUT/$name.report.txt" 2>&1; then
    gvi=$(grep -m1 '^GVI = ' "$OUT/$name.report.txt" || echo "GVI = ?")
    echo "$gvi"
  else
    echo "FAILED (see $name.report.txt)"
  fi
done

python3 "$ROOT/summarize_corpus.py" "$OUT" && echo "Summary: $OUT/ALL_PATHOGENS_SUMMARY.csv"
