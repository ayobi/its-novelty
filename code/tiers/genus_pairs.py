#!/usr/bin/env python3
"""Near-identical sequences filed under different genus names.

For held-out-species queries (their genus stays in the reference), count those whose best hit
reaches --min-identity but lies in another genus, and list the genus pairs involved. With the
query's own species removed, such a hit is a close relative carrying a different genus name:
names for sexual and asexual stages of one fungus, segregate genera, or misidentified references.

  python genus_pairs.py --best lso_genus_best.tsv --taxonomy out_lso_its2/taxonomy.tsv \\
      --min-identity 97 --out results/genus_pairs.tsv
"""
import argparse
from collections import Counter, defaultdict

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--best", required=True, help="best_per_taxon output at genus rank")
    ap.add_argument("--taxonomy", required=True, help="taxonomy.tsv giving each query's genus")
    ap.add_argument("--min-identity", type=float, default=97.0)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    genus = {}
    with open(a.taxonomy) as fh:
        head = fh.readline().rstrip("\n").split("\t"); i, j = head.index("seq_id"), head.index("g")
        for line in fh:
            f = line.rstrip("\n").split("\t"); genus[f[i]] = f[j]
    best = defaultdict(dict)
    with open(a.best) as fh:
        fh.readline()
        for line in fh:
            q, t, s = line.rstrip("\n").split("\t")
            if t != "-":
                best[q][t] = float(s)
    close = other = 0; pairs = Counter(); own_ids = []
    for q, d in best.items():
        if not d or q not in genus:
            continue
        top, s = max(d.items(), key=lambda x: x[1])
        if s < a.min_identity:
            continue
        close += 1
        if top != genus[q]:
            other += 1; pairs[tuple(sorted((genus[q], top)))] += 1; own_ids.append(d.get(genus[q], 0.0))
    own_ids.sort()
    with open(a.out, "w") as o:
        o.write("genus_a\tgenus_b\tqueries\n")
        for (x, y), n in pairs.most_common():
            o.write(f"{x}\t{y}\t{n}\n")
    print(f"[pairs] {close} queries with a best hit >= {a.min_identity:g}%; {other} ({other / max(close, 1):.1%}) "
          f"in another genus, across {len(pairs)} genus pairs; median identity to their own genus "
          f"{own_ids[len(own_ids) // 2] if own_ids else float('nan'):.1f}% -> {a.out}")

if __name__ == "__main__":
    main()
