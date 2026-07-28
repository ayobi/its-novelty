#!/usr/bin/env python3
"""Score rank recovery for held-out queries against a reference search.

Input: a hits file with 3 tab-separated columns (query_id, target_id, score),
where higher score = closer, and taxonomy.tsv (seq_id p c o f g s). For each
query it reports whether the top hit's family (default rank) matches the query's
true family, and writes a per-query TSV carrying the top-hit score so the same
file feeds the novelty-AUROC step later.

Any search tool that can emit (query, target, score) works:
  vsearch   --userfields query+target+id        (score = % identity)
  blastn    -outfmt "6 qseqid sseqid pident"     (score = % identity)
  taxotagger export (see the query script)       (score = cosine similarity)
"""
import argparse

RANK_COL = {"phylum": "p", "class": "c", "order": "o",
            "family": "f", "genus": "g", "species": "s"}

def load_tax(path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {n: i for i, n in enumerate(header)}
        rows = {}
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 2:
                continue
            rows[p[idx["seq_id"]]] = p
    return rows, idx

def query_ids(fasta):
    ids = []
    with open(fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                ids.append(line[1:].strip().split()[0])
    return ids

def best_hits(path):
    best = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            f = line.split("\t")
            if len(f) < 3:
                continue
            q, t = f[0], f[1]
            try:
                s = float(f[2])
            except ValueError:
                continue
            if q not in best or s > best[q][1]:
                best[q] = (t, s)
    return best

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--taxonomy", required=True)
    ap.add_argument("--hits", required=True)
    ap.add_argument("--queries", required=True,
                    help="queries.fasta, defines the full query set incl. no-hit")
    ap.add_argument("--rank", default="family", choices=list(RANK_COL))
    ap.add_argument("--out", required=True, help="per-query TSV output")
    a = ap.parse_args()

    tax, idx = load_tax(a.taxonomy)
    rcol = idx[RANK_COL[a.rank]]
    qids = query_ids(a.queries)
    best = best_hits(a.hits)

    n = len(qids)
    n_hit = n_correct = n_nohit = n_missing = 0
    with open(a.out, "w") as out:
        out.write(f"query_id\ttrue_{a.rank}\tpred_{a.rank}\ttop_target\tscore\tcorrect\n")
        for q in qids:
            trow = tax.get(q)
            true_r = trow[rcol] if trow else ""
            if not true_r:
                n_missing += 1
            if q in best:
                t, s = best[q]
                ttrow = tax.get(t)
                pred_r = ttrow[rcol] if ttrow else ""
                ok = int(bool(true_r) and pred_r == true_r)
                n_hit += 1
                n_correct += ok
                out.write(f"{q}\t{true_r}\t{pred_r}\t{t}\t{s:.4f}\t{ok}\n")
            else:
                n_nohit += 1
                out.write(f"{q}\t{true_r}\t\t\t\t0\n")

    print(f"[score] rank             : {a.rank}")
    print(f"[score] queries          : {n}")
    print(f"[score] with a hit       : {n_hit}")
    print(f"[score] no hit           : {n_nohit}")
    if n_missing:
        print(f"[score] query not in tax : {n_missing}")
    print(f"[score] correct          : {n_correct}")
    print(f"[score] recovery / all   : {n_correct/n:.4f}" if n else "[score] recovery / all   : n/a")
    print(f"[score] recovery / hit   : {n_correct/n_hit:.4f}" if n_hit else "[score] recovery / hit   : n/a")

if __name__ == "__main__":
    main()
