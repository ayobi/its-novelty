#!/usr/bin/env python3
"""Turn conformal_tiers.py --apply output into the verdict table gf_recurrence.sh reads
(asv_id, best_sim, p_value, novel_a05, novel_a10, novel_a20), using the single absence call.

  python apply_to_verdicts.py --apply results/soil_itsx_a0.10_apply_a0.1.tsv --out frey_asv_pvalues_itsx.tsv
"""
import argparse, csv

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--apply", required=True)
ap.add_argument("--out", required=True)
a = ap.parse_args()
n = flagged = 0
with open(a.apply) as fh, open(a.out, "w") as o:
    o.write("asv_id\tbest_sim\tp_value\tnovel_a05\tnovel_a10\tnovel_a20\n")
    for r in csv.DictReader(fh, delimiter="\t"):
        p, b = float(r["absence_p_single"]), float(r["best_identity"])
        o.write(f"{r['query_id']}\t{'NA' if b == 0 else f'{b:.1f}'}\t{p:.4f}\t"
                f"{int(p <= 0.05)}\t{int(p <= 0.10)}\t{int(p <= 0.20)}\n")
        n += 1; flagged += p <= 0.05
print(f"[verdicts] {n} ASVs, {flagged} ({flagged / max(n, 1):.1%}) absent at alpha 0.05 -> {a.out}")
