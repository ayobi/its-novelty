#!/usr/bin/env python3
"""Assign conformal p-values to Frey ASVs against the ITS2 known-calibration null.
Auto-novel (p=0) for ASVs with no reference hit >=50%."""
import numpy as np, sys

def best_scores(path):
    best = {}
    for line in open(path):
        p = line.rstrip("\n").split("\t")
        if len(p) < 3: continue
        try: s = float(p[2])
        except ValueError: continue
        if p[0] not in best or s > best[p[0]]: best[p[0]] = s
    return best

def pvalues(calib_novelty, test_novelty):
    c = np.sort(np.asarray(calib_novelty, float)); n = len(c)
    ge = n - np.searchsorted(c, np.asarray(test_novelty, float), side="left")
    return (1 + ge) / (n + 1)

# calibration null = ALL known (LSO) best-sims, as novelty scores (negated)
known = best_scores("lso_its2.hits")
calib_s = np.array([-v for v in known.values()])

# all 5521 ASV IDs from the fasta; hits give sims, missing = auto-novel
asv_ids = [l[1:].split()[0] for l in open("dada2_out/asvs_its2.fasta") if l.startswith(">")]
asv_hit = best_scores("frey_asvs.hits")

with open("frey_asv_pvalues.tsv", "w") as o:
    o.write("asv_id\tbest_sim\tp_value\tnovel_a05\tnovel_a10\tnovel_a20\n")
    for a in asv_ids:
        if a in asv_hit:
            sim = asv_hit[a]; p = float(pvalues(calib_s, [-sim])[0])
        else:
            sim = float("nan"); p = 0.0   # no hit >=50% -> auto-novel
        o.write(f"{a}\t{sim if sim==sim else 'NA'}\t{p:.4f}"
                f"\t{int(p<=0.05)}\t{int(p<=0.10)}\t{int(p<=0.20)}\n")

n = sum(1 for _ in open("frey_asv_pvalues.tsv")) - 1
import collections
rows = [l.split("\t") for l in open("frey_asv_pvalues.tsv").read().splitlines()[1:]]
for col,name in [(3,"a=0.05"),(4,"a=0.10"),(5,"a=0.20")]:
    nov = sum(int(r[col]) for r in rows)
    print(f"  novel @ {name}: {nov}/{n} ({100*nov/n:.1f}%)")
