#!/usr/bin/env python3
"""Cross-study prevalence: for each Frey SH, count GlobalFungi samples where it occurs.
Then compare prevalence between Frey's novel-verdict vs known-verdict SHs."""
import gzip, csv
import numpy as np

# Frey SH verdicts
verdict = {}
for r in csv.DictReader(open("frey_sh_verdicts.tsv"), delimiter="\t"):
    verdict[r["SH"]] = int(r["sh_verdict"])
frey_shs = set(verdict)
print(f"Frey SHs to look up: {len(frey_shs)}")

gf = gzip.open("GlobalFungi_5_SH_abundance_ITS2.txt.gz", "rt")
header = gf.readline().rstrip("\n").split("\t")
# map column index -> SH, but only for columns that are Frey SHs
col_sh = {i: sh for i, sh in enumerate(header) if sh in frey_shs}
print(f"of those, present as GlobalFungi columns: {len(col_sh)}")

prevalence = {sh: 0 for sh in col_sh.values()}
nsamp = 0
for line in gf:
    nsamp += 1
    parts = line.rstrip("\n").split("\t")
    for i, sh in col_sh.items():
        if i < len(parts) and parts[i] not in ("0", ""):
            prevalence[sh] += 1
gf.close()
print(f"GlobalFungi samples scanned: {nsamp}")

# split prevalence by Frey verdict
nov = [prevalence[sh] for sh in col_sh.values() if verdict[sh] == 1]
kno = [prevalence[sh] for sh in col_sh.values() if verdict[sh] == 0]
nov = np.array(nov); kno = np.array(kno)

print(f"\n=== cross-study prevalence: novel-verdict vs known-verdict SHs ===")
print(f"  novel-verdict SHs: n={len(nov):4d}  median GF prevalence={np.median(nov):8.1f}  mean={nov.mean():8.1f}")
print(f"  known-verdict SHs: n={len(kno):4d}  median GF prevalence={np.median(kno):8.1f}  mean={kno.mean():8.1f}")

# Mann-Whitney U (novel < known expected)
try:
    from scipy.stats import mannwhitneyu
    u, p = mannwhitneyu(nov, kno, alternative="less")
    print(f"  Mann-Whitney U (novel < known): U={u:.0f}  p={p:.3e}")
except ImportError:
    print("  (scipy not available; medians above still tell the story)")

# how many novel-verdict SHs are GF-absent entirely (prevalence 0)?
print(f"\n  novel-verdict SHs absent from GlobalFungi: {int((nov==0).sum())}/{len(nov)}")
print(f"  known-verdict SHs absent from GlobalFungi: {int((kno==0).sum())}/{len(kno)}")

# write per-SH for later
with open("frey_sh_gf_prevalence.tsv","w") as o:
    o.write("SH\tverdict\tgf_prevalence\n")
    for sh in col_sh.values():
        o.write(f"{sh}\t{verdict[sh]}\t{prevalence[sh]}\n")
