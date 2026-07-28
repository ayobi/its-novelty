#!/bin/bash
# gf_recurrence.sh -- cross-study recurrence of soil ASVs in the GlobalFungi ITS2
# variant library, both arms, against the full ASV set.
#
#   ARM=named bash gf_recurrence.sh              # seconds, checks the invocation
#   nohup bash gf_recurrence.sh > gf_recurrence.out 2>&1 &      # both arms
#   ARM=summary bash gf_recurrence.sh            # re-summarise existing hits
#
# Replaces the command-line analysis that produced gf_recur_cov.hits, which had no
# script on disk. Deliberate differences:
#   * database is all 5,521 harmonised ASVs, not the 1,100-record flagged subset, so
#     flagged and non-flagged share one search and one denominator
#   * --id 0.94, since the paper reports 97% and 95% only
#   * --maxaccepts 20 --maxrejects 0 set explicitly. With vsearch's default of one
#     accepted hit, a global variant matching both a flagged and a non-flagged ASV
#     reports only one of them, biasing the very comparison this run exists to make
#   * coverage from vsearch's qcov (percentage of the query covered) rather than
#     alnlen/ql, which counts gaps and can exceed 100%
#   * the named arm gets the same userfields, so both arms share one criterion
#
# Recurrence is ANY global variant matching an ASV at the stated identity over at least
# 90% of the global variant. Not best-hit: the question is whether the sequence occurs
# elsewhere, not whether its closest global relative clears the threshold.

set -euo pipefail

DB="${DB:-dada2_out/asvs_its2.fasta}"
UNK="${UNK:-GlobalFungi_ITS_variants/gf_unknown_trim.fa}"
KNOWN="${KNOWN:-GlobalFungi_ITS_variants/gf_its2_known_trim.fa}"
VERDICTS="${VERDICTS:-frey_asv_pvalues.tsv}"
ABUND="${ABUND:-asv_abundance.tsv}"
THREADS="${THREADS:-32}"
ID="${ID:-0.94}"
ARM="${ARM:-both}"
COV="${COV:-90}"
FIELDS=query+target+id+alnlen+ql+tl+qcov

for f in "$DB" "$VERDICTS" "$ABUND"; do
  [ -s "$f" ] || { echo "missing or empty: $f"; exit 1; }
done
echo "database: $(grep -c '^>' "$DB") ASVs"

run () {   # run <queries> <out.hits> <out.log>
  [ -s "$1" ] || { echo "missing or empty: $1"; exit 1; }
  echo "[$(date +%H:%M:%S)] searching $1  ($(grep -c '^>' "$1") queries)"
  vsearch --usearch_global "$1" --db "$DB" --id "$ID" \
          --strand both --maxaccepts 20 --maxrejects 0 \
          --threads "$THREADS" --userfields "$FIELDS" \
          --userout "$2" --log "$3"
  echo "[$(date +%H:%M:%S)] $2: $(wc -l < "$2") hit rows"
}

summarise () {   # summarise <hits> <arm label>
  [ -s "$1" ] || { echo "no hits file $1, skipping summary"; return; }
  awk -F'\t' -v arm="$2" -v cov="$COV" '
    FNR==1 { fileno++ }                        # positional, not by filename
    fileno==1 && FNR>1 { prev[$1]=$2+0; na++; next }
    fileno==2 && FNR>1 { flag[$1]=$4+0; nv++; next }
    fileno==3 {
      if (NF < 7) { badfmt++; next }
      if ($7+0 < cov) next
      if ($3+0 >= 97) r97[$2]=1
      if ($3+0 >= 95) r95[$2]=1
    }
    END {
      if (!na || !nv) {
        printf "ERROR: read %d abundance and %d verdict rows; check argument order\n",
               na+0, nv+0 > "/dev/stderr"
        exit 1
      }
      if (badfmt) {
        printf "ERROR: %d hit rows have fewer than 7 fields, so carry no qcov. This is\n",
               badfmt > "/dev/stderr"
        printf "       the old six-column gf_recur_cov.hits, not output of this script.\n" \
               > "/dev/stderr"
        exit 1
      }
      for (a in flag) {
        f = flag[a]; p = (prev[a] >= 2) ? "multi" : "single"
        n[f]++;     h97[f]     += (a in r97); h95[f]     += (a in r95)
        n[f"_"p]++; h97[f"_"p] += (a in r97); h95[f"_"p] += (a in r95)
      }
      ord[1]="1";        lbl["1"]="flagged"
      ord[2]="1_multi";  lbl["1_multi"]="  in >=2 samples"
      ord[3]="1_single"; lbl["1_single"]="  in 1 sample"
      ord[4]="0";        lbl["0"]="non-flagged"
      ord[5]="0_multi";  lbl["0_multi"]="  in >=2 samples"
      ord[6]="0_single"; lbl["0_single"]="  in 1 sample"
      printf "\n=== %s arm: any hit, >=%s%% coverage of the global variant ===\n", arm, cov
      printf "  %-20s %6s %8s %8s %8s %8s\n", "group", "n", ">=97", "%", ">=95", "%"
      for (i = 1; i <= 6; i++) {
        g = ord[i]; if (!(g in n)) continue
        printf "  %-20s %6d %8d %7.1f%% %8d %7.1f%%\n", lbl[g], n[g], h97[g]+0,
               100*(h97[g]+0)/n[g], h95[g]+0, 100*(h95[g]+0)/n[g]
      }
    }' "$ABUND" "$VERDICTS" "$1"
}

case "$ARM" in
  named)   run "$KNOWN" gf_known_vs_all.hits gf_known_vs_all.log
           summarise gf_known_vs_all.hits "named representative" ;;
  unnamed) run "$UNK" gf_unknown_vs_all.hits gf_unknown_vs_all.log
           summarise gf_unknown_vs_all.hits "unnamed variant" ;;
  summary) summarise gf_unknown_vs_all.hits "unnamed variant"
           summarise gf_known_vs_all.hits "named representative" ;;
  both)    run "$KNOWN" gf_known_vs_all.hits gf_known_vs_all.log
           summarise gf_known_vs_all.hits "named representative"
           run "$UNK" gf_unknown_vs_all.hits gf_unknown_vs_all.log
           summarise gf_unknown_vs_all.hits "unnamed variant" ;;
  *)       echo "ARM must be named, unnamed, summary or both"; exit 1 ;;
esac

echo
echo "note: the named arm searches species-hypothesis representatives, the unnamed arm"
echo "individual variants. The two rates are not symmetric and must not be reported"
echo "side by side without saying so."
