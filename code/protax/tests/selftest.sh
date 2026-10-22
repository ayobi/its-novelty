#!/usr/bin/env bash
# End-to-end plumbing test on synthetic data: prep_split -> vsearch identity hits ->
# PROTAX model build -> PROTAX -> scores -> comparison -> soil cross-tab, plus a build
# without the PROTAX taxonomy. Run it once before the real benchmark to check that
# vsearch, Perl and the PROTAX-fungi scripts are wired up. Takes about a minute.
#
#   PROTAX_SCRIPTS=protaxfungi/protaxscripts PROTAX_PARAMS=protaxfungi/model \
#     bash code/protax/tests/selftest.sh [workdir]
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
CODE=$(cd "$HERE/../.." && pwd)                 # repo code/ (prep_split.py lives here)
PX=$(cd "$HERE/.." && pwd)                      # code/protax/
: "${PROTAX_SCRIPTS:?set PROTAX_SCRIPTS to protaxfungi/protaxscripts}"
: "${PROTAX_PARAMS:?set PROTAX_PARAMS to protaxfungi/model (params.its2.level2..7)}"
export PROTAX_SCRIPTS VSEARCH=${VSEARCH:-vsearch} THREADS=${THREADS:-2}
W=${1:-$(mktemp -d)}; mkdir -p "$W"; cd "$W"
echo "[selftest] working in $W"

python3 "$HERE/make_synthetic_unite.py" --out synth
for mode in lgo lso; do
  python3 "$CODE/prep_split.py" -i synth/unite_synth.fasta -o "out_${mode}" --mode "$mode" --seed 0 > "prep_${mode}.log"
  "$VSEARCH" --usearch_global "out_${mode}/queries.fasta" --db "out_${mode}/refs.fasta" --id 0.5 \
    --maxaccepts 1 --top_hits_only --userfields query+target+id --userout "${mode}.hits" --quiet
  python3 "$PX/build_protax_model.py" --split-dir "out_${mode}" --protax-taxonomy synth/official_taxonomy \
    --params "$PROTAX_PARAMS" --out "model_${mode}"
  bash "$PX/run_protax.sh" "model_${mode}" "out_${mode}/queries.fasta" "protax_${mode}"
  python3 "$PX/protax_scores.py" --protax-dir "protax_${mode}" --level 6 \
    --taxonomy "model_${mode}/taxonomy" --out "protax_${mode}/scores_genus.tsv"
done

python3 "$PX/compare_protax_conformal.py" \
  --known-queries out_lso/queries.fasta --novel-queries out_lgo/queries.fasta \
  --known-hits lso.hits --novel-hits lgo.hits \
  --known-protax protax_lso/scores_genus.tsv --novel-protax protax_lgo/scores_genus.tsv \
  --mixes 0.05,0.2,0.5 --mix-size 60 --draws 50 --splits 20 --out-prefix results/selftest

# soil branch: full synthetic pool as reference, conformal p-values from the LSO knowns
python3 "$PX/build_protax_model.py" --unite-fasta synth/unite_synth.fasta \
  --protax-taxonomy synth/official_taxonomy --params "$PROTAX_PARAMS" --out model_pool
bash "$PX/run_protax.sh" model_pool synth/soil_asvs.fasta protax_asvs
python3 "$PX/protax_scores.py" --protax-dir protax_asvs --out protax_asvs/scores_genus.tsv
"$VSEARCH" --usearch_global synth/soil_asvs.fasta --db out_lso/refs.fasta --id 0.5 --maxaccepts 1 \
  --top_hits_only --userfields query+target+id --userout asvs.hits --quiet
python3 - <<'EOF'
import numpy as np
cal = {}
for l in open("lso.hits"):
    q, t, s = l.split("\t")[:3]; cal[q] = max(cal.get(q, 0), float(s))
c = np.sort([-v for v in cal.values()])
hit = {}
for l in open("asvs.hits"):
    q, t, s = l.split("\t")[:3]; hit[q] = max(hit.get(q, 0), float(s))
ids = [l[1:].split()[0] for l in open("synth/soil_asvs.fasta") if l.startswith(">")]
with open("asv_pvalues.tsv", "w") as o:
    o.write("asv_id\tbest_sim\tp_value\n")
    for a in ids:
        s = -hit.get(a, 0.0); p = (1 + len(c) - np.searchsorted(c, s, side="left")) / (len(c) + 1)
        o.write(f"{a}\t{-s}\t{p:.4f}\n")
EOF
python3 "$PX/crosstab_soil.py" --flags asv_pvalues.tsv --protax protax_asvs/scores_genus.tsv \
  --out results/soil_crosstab.tsv

# the no-taxonomy fallback must at least run
python3 "$PX/build_protax_model.py" --split-dir out_lgo --params "$PROTAX_PARAMS" --out model_lgo_fallback > /dev/null
bash "$PX/run_protax.sh" model_lgo_fallback out_lgo/queries.fasta protax_lgo_fallback > /dev/null
python3 "$PX/protax_scores.py" --protax-dir protax_lgo_fallback --out protax_lgo_fallback/scores.tsv > /dev/null
echo "[selftest] PASS: every stage ran; outputs in $W"
