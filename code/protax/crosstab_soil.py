#!/usr/bin/env python3
"""Cross-tabulate PROTAX genus-rank absence calls against the conformal verdicts for the
soil ASVs (data/frey_asv_pvalues.tsv), at each PROTAX threshold.

Run PROTAX on the ASVs with a model built from the same reference the ASV identity hits
were computed against (build_protax_model.py --unite-fasta for a full pool, or
--split-dir for a split's refs), so the two calls see the same sequences.

  python crosstab_soil.py --flags data/frey_asv_pvalues.tsv \
      --protax protax_asvs/scores_genus.tsv --alpha 0.05 --out results/soil_crosstab.tsv
"""
import argparse

def read_tsv(path, key, cols):
    out = {}
    with open(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        ix = {h: i for i, h in enumerate(head)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            out[f[ix[key]]] = {c: f[ix[c]] for c in cols}
    return out

def kappa(a, b, c, d):
    n = a + b + c + d
    po = (a + d) / n
    pe = ((a + b) * (a + c) + (c + d) * (b + d)) / n ** 2
    return float("nan") if pe == 1 else (po - pe) / (1 - pe)

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--flags", required=True, help="asv_id, p_value columns (conformal)")
    ap.add_argument("--protax", required=True, help="protax_scores.py output for the ASVs")
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--thresholds", default="0.5,0.9,0.95")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    fl = read_tsv(a.flags, "asv_id", ["p_value"])
    pr = read_tsv(a.protax, "query_id", ["p_absent", "top_taxon", "top_prob"])
    ids = [x for x in fl if x in pr]
    if len(ids) < len(fl):
        print(f"[crosstab] {len(fl) - len(ids)} ASVs missing from the PROTAX table")
    conf = {x: float(fl[x]["p_value"]) <= a.alpha for x in ids}
    print(f"[crosstab] {len(ids)} ASVs; conformal absent at alpha={a.alpha:g}: {sum(conf.values())}")
    print(f"   {'PROTAX t':>8} {'both':>6} {'conf only':>10} {'PROTAX only':>12} {'neither':>8} {'agree':>7} {'kappa':>7}")
    for t in [float(x) for x in a.thresholds.split(",")]:
        pa = {x: float(pr[x]["p_absent"]) >= t for x in ids}
        A = sum(conf[x] and pa[x] for x in ids); B = sum(conf[x] and not pa[x] for x in ids)
        C = sum(pa[x] and not conf[x] for x in ids); D = len(ids) - A - B - C
        print(f"   {t:8.2f} {A:6d} {B:10d} {C:12d} {D:8d} {(A + D) / len(ids):7.3f} {kappa(A, B, C, D):7.3f}")
    with open(a.out, "w") as o:
        o.write("asv_id\tconformal_p\tconformal_absent\tprotax_p_absent\tprotax_top_taxon\tprotax_top_prob\n")
        for x in ids:
            o.write(f"{x}\t{fl[x]['p_value']}\t{int(conf[x])}\t{pr[x]['p_absent']}\t"
                    f"{pr[x]['top_taxon']}\t{pr[x]['top_prob']}\n")
    print(f"[crosstab] per-ASV table -> {a.out}")

if __name__ == "__main__":
    main()
