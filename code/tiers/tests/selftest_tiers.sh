#!/usr/bin/env bash
# Plumbing test for the three-level calls at species, genus and family rank, on synthetic
# data (numbers mean nothing). Runs prep_split for all four split modes, a many-hit vsearch
# search per split, best_per_taxon at each rank, conformal_tiers, and the apply mode.
#
#   bash code/tiers/tests/selftest_tiers.sh [workdir]
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
CODE=$(cd "$HERE/../.." && pwd)
T=$(cd "$HERE/.." && pwd)
VSEARCH=${VSEARCH:-vsearch}
W=${1:-$(mktemp -d)}; mkdir -p "$W"; cd "$W"
echo "[selftest] working in $W"
python3 "$CODE/protax/tests/make_synthetic_unite.py" --out synth > /dev/null

for mode in loo lso lgo lfo; do
  python3 "$CODE/prep_split.py" -i synth/unite_synth.fasta -o "out_$mode" --mode "$mode" --seed 0 > /dev/null
done
best() {   # split rank out
  "$VSEARCH" --usearch_global "out_$1/queries.fasta" --db "out_$1/refs.fasta" --id 0.5 \
    --maxaccepts 0 --maxrejects 0 --quiet --userfields query+target+id+qcov+tcov+ql+tl --userout /dev/stdout \
  | python3 "$T/best_per_taxon.py" --taxonomy "out_$1/taxonomy.tsv" --rank "$2" \
      --queries "out_$1/queries.fasta" --out "$3"
}
best loo s loo_species.tsv;  best lso s lso_species.tsv     # species rank
best lso g lso_genus.tsv;    best lgo g lgo_genus.tsv       # genus rank
best lgo f lgo_family.tsv;   best lfo f lfo_family.tsv      # family rank

python3 "$T/conformal_tiers.py" --rank s --known-best loo_species.tsv --known-taxonomy out_loo/taxonomy.tsv \
  --heldout-best lso_species.tsv --splits 20 --out-prefix results/tiers_species > /dev/null
python3 "$T/conformal_tiers.py" --rank g --known-best lso_genus.tsv --known-taxonomy out_lso/taxonomy.tsv \
  --heldout-best lgo_genus.tsv --splits 20 --out-prefix results/tiers_genus
python3 "$T/conformal_tiers.py" --rank f --known-best lgo_family.tsv --known-taxonomy out_lgo/taxonomy.tsv \
  --heldout-best lfo_family.tsv --splits 20 --out-prefix results/tiers_family > /dev/null

# apply mode: label some "soil" ASVs at genus rank against the LSO reference
"$VSEARCH" --usearch_global synth/soil_asvs.fasta --db out_lso/refs.fasta --id 0.5 --maxaccepts 0 \
  --maxrejects 0 --quiet --userfields query+target+id+qcov+tcov+ql+tl --userout /dev/stdout \
| python3 "$T/best_per_taxon.py" --taxonomy out_lso/taxonomy.tsv --rank g --queries synth/soil_asvs.fasta \
    --out asv_genus.tsv
python3 "$T/conformal_tiers.py" --rank g --known-best lso_genus.tsv --known-taxonomy out_lso/taxonomy.tsv \
  --apply-best asv_genus.tsv --splits 5 --out-prefix results/tiers_asvs | tail -2
echo "[selftest] PASS: every stage ran; outputs in $W"
