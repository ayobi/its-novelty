#!/usr/bin/env bash
# Regenerate every benchmark and soil number in the revision through one pipeline:
#   exhaustive vsearch (--maxaccepts 0 --maxrejects 0, --id 0.5), hits kept only when the
#   alignment covers >= 80% of the shorter sequence, queries without hits kept (identity 0).
#
# Produces, under $OUT/results:
#   protax_vs_conformal_*      PROTAX-fungi comparison (identity side refreshed) + figure
#   tiers_{species,genus,genus_species_present,family}_*   three-level calls + figure
#   genus_pairs.tsv            near-identical sequences filed under other genus names
#   soil_itsx_*                soil calls, verdict table, and table3.txt (GlobalFungi recurrence)
#
# Each step is skipped, with a message, when its inputs are missing. Searches are cached in
# $OUT/hits, so a rerun only repeats the analyses. Takes about half an hour on 32 threads.
#
#   bash code/refresh_numbers.sh
#   env: CODE SPLITS NEWSPLITS PROTAX SOIL GFDIR THREADS OUT (defaults below)
set -euo pipefail
CODE=${CODE:-$(cd "$(dirname "$0")" && pwd)}
SPLITS=${SPLITS:-$HOME/Downloads}          # out_lso_its2, out_lgo_its2 (the paper's splits)
NEWSPLITS=${NEWSPLITS:-$HOME/tiers_run}    # out_loo_its2, out_lfo_its2
PROTAX=${PROTAX:-$HOME/protax_run}         # protax_lso/, protax_lgo/ with scores_genus.tsv
SOIL=${SOIL:-$CODE/../data/asvs_its2_itsx.fasta}   # ITSxRust ITS2 of the soil ASVs
GFDIR=${GFDIR:-$HOME/Downloads}            # gf_*_vs_all.hits and dada2_out/ for gf_recurrence.sh
THREADS=${THREADS:-32}
OUT=${OUT:-refresh}
T=$CODE/tiers; P=$CODE/protax
ABUND=$CODE/../data/asv_abundance.tsv

command -v vsearch > /dev/null || { echo "vsearch is not on the PATH (conda activate base?)"; exit 1; }
mkdir -p "$OUT/hits" "$OUT/best" "$OUT/results"
say() { echo "[$(date +%H:%M:%S)] $*"; }
have() { for f in "$@"; do [ -s "$f" ] || { say "skipping: missing $f"; return 1; }; done; }

search() {   # name queries refs
  if [ -s "$OUT/hits/$1.tsv.gz" ]; then say "search $1: cached"; return; fi
  say "search $1: $(grep -c '>' "$2") queries against $(grep -c '>' "$3") references"
  vsearch --usearch_global "$2" --db "$3" --id 0.5 --maxaccepts 0 --maxrejects 0 --threads "$THREADS" \
    --quiet --userfields query+target+id+qcov+tcov+ql+tl --userout /dev/stdout | gzip > "$OUT/hits/$1.tsv.gz.part"
  mv "$OUT/hits/$1.tsv.gz.part" "$OUT/hits/$1.tsv.gz"
}
best() {     # name rank taxonomy queries
  zcat "$OUT/hits/$1.tsv.gz" | python3 "$T/best_per_taxon.py" --taxonomy "$3" --rank "$2" \
    --queries "$4" --out "$OUT/best/$1_$2.tsv"
}
besthits() { # best table -> (query, taxon, identity) rows for compare_protax_conformal.py
  awk -F'\t' 'NR > 1 && $2 != "-"' "$1" > "$2"
}

LSO=$SPLITS/out_lso_its2; LGO=$SPLITS/out_lgo_its2; LOO=$NEWSPLITS/out_loo_its2; LFO=$NEWSPLITS/out_lfo_its2

# ---------------------------------------------------------------- searches and per-taxon bests
have "$LSO/queries.fasta" "$LGO/queries.fasta" || exit 1
search lso "$LSO/queries.fasta" "$LSO/refs.fasta"; best lso g "$LSO/taxonomy.tsv" "$LSO/queries.fasta"
best lso s "$LSO/taxonomy.tsv" "$LSO/queries.fasta"
search lgo "$LGO/queries.fasta" "$LGO/refs.fasta"; best lgo g "$LGO/taxonomy.tsv" "$LGO/queries.fasta"
best lgo f "$LGO/taxonomy.tsv" "$LGO/queries.fasta"
if have "$LOO/queries.fasta"; then
  search loo "$LOO/queries.fasta" "$LOO/refs.fasta"; best loo s "$LOO/taxonomy.tsv" "$LOO/queries.fasta"
  best loo g "$LOO/taxonomy.tsv" "$LOO/queries.fasta"
fi
if have "$LFO/queries.fasta"; then
  search lfo "$LFO/queries.fasta" "$LFO/refs.fasta"; best lfo f "$LFO/taxonomy.tsv" "$LFO/queries.fasta"
fi

# ---------------------------------------------------------------- PROTAX-fungi comparison
if have "$PROTAX/protax_lso/scores_genus.tsv" "$PROTAX/protax_lgo/scores_genus.tsv"; then
  say "PROTAX comparison"
  besthits "$OUT/best/lso_g.tsv" "$OUT/best/lso_best.hits"; besthits "$OUT/best/lgo_g.tsv" "$OUT/best/lgo_best.hits"
  python3 "$P/compare_protax_conformal.py" \
    --known-queries "$LSO/queries.fasta" --novel-queries "$LGO/queries.fasta" \
    --known-hits "$OUT/best/lso_best.hits" --novel-hits "$OUT/best/lgo_best.hits" \
    --known-protax "$PROTAX/protax_lso/scores_genus.tsv" --novel-protax "$PROTAX/protax_lgo/scores_genus.tsv" \
    --out-prefix "$OUT/results/protax_vs_conformal" > /dev/null
  if [ -s "$P/make_protax_figure.py" ]; then
    python3 "$P/make_protax_figure.py" --results "$OUT/results/protax_vs_conformal" \
      --out "$OUT/results/fig_protax_comparison" > /dev/null
  fi
fi

# ---------------------------------------------------------------- three-level calls
tiers() { python3 "$T/conformal_tiers.py" "$@" > /dev/null; }
say "three-level calls"
tiers --rank g --known-best "$OUT/best/lso_g.tsv" --known-taxonomy "$LSO/taxonomy.tsv" \
  --heldout-best "$OUT/best/lgo_g.tsv" \
  $( [ -s "$OUT/best/loo_g.tsv" ] && echo --transfer-best "$OUT/best/loo_g.tsv" --transfer-taxonomy "$LOO/taxonomy.tsv" ) \
  --transfer-label "species present" --out-prefix "$OUT/results/tiers_genus"
if [ -s "$OUT/best/loo_s.tsv" ]; then
  tiers --rank s --known-best "$OUT/best/loo_s.tsv" --known-taxonomy "$LOO/taxonomy.tsv" \
    --heldout-best "$OUT/best/lso_s.tsv" --out-prefix "$OUT/results/tiers_species"
  tiers --rank g --known-best "$OUT/best/loo_g.tsv" --known-taxonomy "$LOO/taxonomy.tsv" \
    --heldout-best "$OUT/best/lgo_g.tsv" --out-prefix "$OUT/results/tiers_genus_species_present"
fi
if [ -s "$OUT/best/lfo_f.tsv" ]; then
  tiers --rank f --known-best "$OUT/best/lgo_f.tsv" --known-taxonomy "$LGO/taxonomy.tsv" \
    --heldout-best "$OUT/best/lfo_f.tsv" --out-prefix "$OUT/results/tiers_family"
fi
if [ -s "$OUT/results/tiers_species_tiers.tsv" ] && [ -s "$OUT/results/tiers_family_tiers.tsv" ]; then
  python3 "$T/make_tiers_figure.py" --species "$OUT/results/tiers_species" --genus "$OUT/results/tiers_genus" \
    --family "$OUT/results/tiers_family" --alpha 0.10 --out "$OUT/results/fig_tiers" > /dev/null
fi
python3 "$T/genus_pairs.py" --best "$OUT/best/lso_g.tsv" --taxonomy "$LSO/taxonomy.tsv" \
  --out "$OUT/results/genus_pairs.tsv" | tee "$OUT/results/genus_pairs_summary.txt"

# ---------------------------------------------------------------- soil
if have "$SOIL"; then
  say "soil ASVs"
  cat "$LGO/refs.fasta" "$LGO/queries.fasta" > "$OUT/best/backbone_its2.fasta"
  search soil "$SOIL" "$OUT/best/backbone_its2.fasta"; best soil g "$LGO/taxonomy.tsv" "$SOIL"
  for a in 0.05 0.10 0.20; do
    python3 "$T/conformal_tiers.py" --rank g --known-best "$OUT/best/lso_g.tsv" --known-taxonomy "$LSO/taxonomy.tsv" \
      --apply-best "$OUT/best/soil_g.tsv" --apply-alpha "$a" --splits 1 --out-prefix "$OUT/results/soil_itsx_a$a" \
      | { grep -E "^Applied|single absence" || true; } | tee -a "$OUT/results/soil_summary.txt"
  done
  python3 "$T/apply_to_verdicts.py" --apply "$OUT/results/soil_itsx_a0.10_apply_a0.1.tsv" \
    --out "$OUT/results/frey_asv_pvalues_itsx.tsv" | tee -a "$OUT/results/soil_summary.txt"
  if have "$GFDIR/gf_unknown_vs_all.hits" "$GFDIR/gf_known_vs_all.hits" "$ABUND"; then
    V=$(cd "$OUT/results" && pwd)/frey_asv_pvalues_itsx.tsv; A=$(cd "$(dirname "$ABUND")" && pwd)/$(basename "$ABUND")
    (cd "$GFDIR" && ARM=summary VERDICTS="$V" ABUND="$A" bash "$CODE/gf_recurrence.sh") > "$OUT/results/table3.txt"
    say "Table 3 -> $OUT/results/table3.txt"
  fi
fi

say "done. Key summaries:"
ls "$OUT"/results/*summary*.txt "$OUT/results/table3.txt" 2> /dev/null | sed 's/^/  /' || true
