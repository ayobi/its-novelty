#!/usr/bin/env bash
# Run PROTAX-fungi (v1 classifier, per-rank parameters) on a query FASTA against a model
# directory from build_protax_model.py.
#
# Follows UNITE's own PlutoF deployment (github.com/TU-NHM/protax_fungi_plutof_pub,
# run_protax.sh): vsearch --usearch_global --id 0.75 --maxaccepts 1000 --strand both,
# vsearch --sintax, sintax2sa.pl, then classify.pl level by level from the Fungi node.
# Three deliberate differences:
#   PROB_TH  default 1e-4 (PlutoF: 0.01). It only truncates what classify.pl prints; the
#            probabilities are unchanged. A lower cut stops many small dropped "known"
#            branches from inflating the absence probability.
#   CHUNKS   queries are classified in parallel chunks and concatenated in input order.
#   SINTAX_SEED  default 1. PlutoF leaves vsearch's SINTAX seed random, so its results
#            vary slightly between runs; a fixed seed makes this harness reproducible.
#
# usage:  PROTAX_SCRIPTS=path/to/protaxfungi/protaxscripts \
#           bash run_protax.sh MODEL_DIR QUERY.fasta OUT_DIR
# env:    VSEARCH=vsearch  REGION=its2  THREADS=4  CHUNKS=$THREADS  PROB_TH=1e-4
#         MAXLEVEL=6 (genus; 7 adds species)  SINTAX_SEED=1
set -euo pipefail

if [ $# -ne 3 ]; then
  sed -n '2,19p' "$0"; exit 1
fi
MODEL=$(cd "$1" && pwd); QUERY=$2; OUT=$3
PROTAX=${PROTAX_SCRIPTS:?set PROTAX_SCRIPTS to protaxfungi/protaxscripts}
PROTAX=$(cd "$PROTAX" && pwd)
VSEARCH=${VSEARCH:-vsearch}
REGION=${REGION:-its2}
THREADS=${THREADS:-4}
CHUNKS=${CHUNKS:-$THREADS}
PROB_TH=${PROB_TH:-1e-4}
MAXLEVEL=${MAXLEVEL:-6}
START=$(cat "$MODEL/startnode")

for f in classify.pl sintax2sa.pl testsample2init.pl nameprob.pl; do
  [ -s "$PROTAX/$f" ] || { echo "[protax] missing $PROTAX/$f" >&2; exit 1; }
done
mkdir -p "$OUT"
OUT=$(cd "$OUT" && pwd)
echo "[protax] model=$MODEL query=$QUERY out=$OUT region=$REGION prob_th=$PROB_TH levels=2..$MAXLEVEL"

# databases, built once per model
[ -s "$MODEL/$REGION.udb" ] || \
  "$VSEARCH" --makeudb_usearch "$MODEL/$REGION.fa" --output "$MODEL/$REGION.udb" --quiet
[ -s "$MODEL/sintax$REGION.udb" ] || \
  "$VSEARCH" --makeudb_usearch "$MODEL/sintax${REGION}train.fa" --output "$MODEL/sintax$REGION.udb" --quiet

# queries: first header token is the id everywhere downstream
awk '/^>/{print $1; next}{print}' "$QUERY" > "$OUT/query.fa"
grep '^>' "$OUT/query.fa" | cut -c2- > "$OUT/query.ids"
NQ=$(wc -l < "$OUT/query.ids")

"$VSEARCH" --usearch_global "$OUT/query.fa" --db "$MODEL/$REGION.udb" --id 0.75 \
  --maxaccepts 1000 --strand both --userfields query+target+id \
  --userout "$OUT/query.m8" --threads "$THREADS" --quiet
# SINTAX bootstraps with a random generator; vsearch seeds it randomly unless told
# otherwise, which makes PROTAX's probabilities vary run to run. Fix the seed, one thread.
"$VSEARCH" --sintax "$OUT/query.fa" --db "$MODEL/sintax$REGION.udb" \
  --tabbedout "$OUT/query.sintax" --strand both --randseed "${SINTAX_SEED:-1}" --threads 1 --quiet
# a query with no SINTAX prediction must reach sintax2sa.pl as a bare id
awk -F'\t' '$2=="" {print $1; next} {print $1"\t"$2}' "$OUT/query.sintax" > "$OUT/query.sintax.clean"
perl "$PROTAX/sintax2sa.pl" "$MODEL/taxonomy.ascii7" "$OUT/query.sintax.clean" > "$OUT/query.sasintax"

# classify in parallel chunks
rm -rf "$OUT/chunks"; mkdir -p "$OUT/chunks"
split -n "l/$CHUNKS" -d -a 4 "$OUT/query.ids" "$OUT/chunks/ids."
classify_chunk() {
  local ids=$1 dir
  dir="$(dirname "$1")/c.$(basename "$1" | sed 's/^ids\.//')"
  [ -s "$ids" ] || return 0
  mkdir -p "$dir"
  awk 'NR==FNR{k[$1]=1; next} ($1 in k)' "$ids" "$OUT/query.m8" > "$dir/query.m8"
  awk 'NR==FNR{k[$1]=1; next} ($1 in k)' "$ids" "$OUT/query.sasintax" > "$dir/query.sasintax"
  perl "$PROTAX/testsample2init.pl" "$START" "$ids" > "$dir/query1.logprob"
  for LEVEL in $(seq 2 "$MAXLEVEL"); do
    perl "$PROTAX/classify.pl" "$dir/query$((LEVEL-1)).logprob" "$MODEL/tax$LEVEL" \
      "$MODEL/ref.tax$LEVEL" "$MODEL/rseqs$LEVEL" "$MODEL/params.$REGION.level$LEVEL" \
      "$dir/query.m8" "$dir/query.sasintax" 0 "$PROB_TH" "$dir/query$LEVEL.logprob" 0 \
      2>> "$dir/classify.log"
  done
}
export -f classify_chunk
export OUT MODEL PROTAX START MAXLEVEL REGION PROB_TH
ls "$OUT"/chunks/ids.* | xargs -P "$CHUNKS" -I{} bash -c 'classify_chunk "$@"' _ {}

for LEVEL in $(seq 1 "$MAXLEVEL"); do
  : > "$OUT/query$LEVEL.logprob"
  for d in "$OUT"/chunks/c.*; do
    [ -d "$d" ] && cat "$d/query$LEVEL.logprob" >> "$OUT/query$LEVEL.logprob"
  done
  n=$(wc -l < "$OUT/query$LEVEL.logprob")
  [ "$n" -eq "$NQ" ] || { echo "[protax] level $LEVEL: $n lines for $NQ queries" >&2; exit 1; }
  if [ "$LEVEL" -ge 2 ]; then
    perl "$PROTAX/nameprob.pl" "$MODEL/taxonomy" "$OUT/query$LEVEL.logprob" > "$OUT/query$LEVEL.nameprob"
  fi
done
cat "$OUT"/chunks/c.*/classify.log > "$OUT/classify.log" 2>/dev/null || true
rm -rf "$OUT/chunks"
echo "[protax] classified $NQ queries; logprob/nameprob for levels 2..$MAXLEVEL in $OUT/"
