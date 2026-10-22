#!/usr/bin/env python3
"""PROTAX probability that a query's taxon at a given rank is absent from the reference.

classify.pl runs level by level. At level L it prints, per query, the level-L taxa that
have reference sequences, and for each level-(L-1) node it started from a "<node>,unk"
entry: the probability that the query sits below that node but in a taxon without
reference sequences, whether named in the taxonomy or not. Mass sent to unk at an
earlier level is not carried forward. So at rank L:

  P(absent at L) = 1 - (sum of level-L entries with references)   upper bound
                 =     (sum of unk entries at levels 2..L)         lower bound

The bounds differ only by mass classify.pl did not print (run_protax.sh PROB_TH); both are
written, with p_absent the first. A query with no level-L entries at all is absent with
probability 1.

  python protax_scores.py --protax-dir protax_lgo --level 6 \
      --taxonomy model_lgo_its2/taxonomy --out protax_lgo/scores_genus.tsv
"""
import argparse, math, os

def read_logprob(path):
    out = {}
    with open(path) as fh:
        for line in fh:
            f = line.split()
            if f:
                out[f[0]] = [(f[i], float(f[i + 1])) for i in range(1, len(f) - 1, 2)]
    return out

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--protax-dir", required=True, help="run_protax.sh output directory")
    ap.add_argument("--level", type=int, default=6, help="6 = genus (default), 5 = family, 7 = species")
    ap.add_argument("--taxonomy", help="model taxonomy file, to name the top taxon")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    ids = [l.split()[0] for l in open(os.path.join(a.protax_dir, "query.ids")) if l.strip()]
    levels = {L: read_logprob(os.path.join(a.protax_dir, f"query{L}.logprob")) for L in range(2, a.level + 1)}
    names = {}
    if a.taxonomy:
        with open(a.taxonomy) as fh:
            for line in fh:
                f = line.split()
                if len(f) >= 4:
                    names[f[0]] = f[3]
    n_trunc = 0
    with open(a.out, "w") as o:
        o.write("query_id\tp_absent\tp_absent_lower\ttruncated\ttop_taxon\ttop_prob\n")
        for q in ids:
            known = [(n, math.exp(lp)) for n, lp in levels[a.level].get(q, []) if not n.endswith(",unk")]
            unk = sum(math.exp(lp) for L in levels for n, lp in levels[L].get(q, []) if n.endswith(",unk"))
            p_known = min(1.0, sum(p for _, p in known))
            p_abs, p_low = 1.0 - p_known, min(1.0, unk)
            trunc = max(0.0, p_abs - p_low)
            n_trunc += trunc > 0.01
            top, tp = max(known, key=lambda x: x[1]) if known else ("", 0.0)
            o.write(f"{q}\t{p_abs:.6g}\t{p_low:.6g}\t{trunc:.3g}\t{names.get(top, top)}\t{tp:.6g}\n")
    print(f"[scores] {len(ids)} queries, level {a.level}; {n_trunc} with >0.01 unprinted mass "
          f"(lower PROB_TH if that number is large) -> {a.out}")

if __name__ == "__main__":
    main()
