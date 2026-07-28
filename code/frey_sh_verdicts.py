#!/usr/bin/env python3
"""Frey-side of cross-study C2: map each ASV to its best-hit SH, attach the
novelty verdict, collapse to one verdict per SH."""
import csv, collections

# accession -> SH
acc2sh = {}
for line in open("acc2sh.tsv"):
    a, sh = line.rstrip("\n").split("\t")
    acc2sh[a] = sh

# ASV -> best-hit accession
asv2acc = {}
for line in open("frey_asvs.hits"):
    p = line.rstrip("\n").split("\t")
    asv2acc[p[0]] = p[1]

# ASV -> novelty flag at a=0.05  (col index 3 in frey_asv_pvalues.tsv)
asv_novel = {}
rows = list(csv.reader(open("frey_asv_pvalues.tsv"), delimiter="\t"))
hdr = rows[0]
for r in rows[1:]:
    asv_novel[r[0]] = int(r[3])   # novel_a05

# collapse ASVs -> SH.  An SH's verdict: novel only if ALL its ASVs are novel?
# No — an SH that any ASV calls KNOWN is known (it has a real reference match).
# So SH is 'known' if any hitting ASV is known; novel only if all novel.
sh_flags = collections.defaultdict(list)
no_hit_novel = 0
for asv, nov in asv_novel.items():
    acc = asv2acc.get(asv)
    if acc is None:            # no reference hit -> auto-novel, no SH to map
        no_hit_novel += 1
        continue
    sh = acc2sh.get(acc)
    if sh is None:
        continue               # accession not in lookup (shouldn't happen)
    sh_flags[sh].append(nov)

with open("frey_sh_verdicts.tsv", "w") as o:
    o.write("SH\tn_asvs\tn_novel\tsh_verdict\n")
    for sh, flags in sh_flags.items():
        # SH known (0) if any ASV known; novel (1) only if all ASVs novel
        verdict = 0 if any(f == 0 for f in flags) else 1
        o.write(f"{sh}\t{len(flags)}\t{sum(flags)}\t{verdict}\n")

n_sh = len(sh_flags)
n_novel_sh = sum(1 for f in sh_flags.values() if all(x==1 for x in f))
print(f"Frey ASVs with a mappable SH: {sum(len(v) for v in sh_flags.values())}")
print(f"ASVs with no hit (auto-novel, no SH): {no_hit_novel}")
print(f"distinct SH hit by Frey: {n_sh}")
print(f"  of which novel-verdict: {n_novel_sh} ({100*n_novel_sh/n_sh:.1f}%)")
