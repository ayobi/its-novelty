#!/usr/bin/env python3
"""Best identity of each query to each reference taxon at one rank, from a many-hit vsearch run.

The three-level calls need every query's best identity to every candidate taxon (genus,
family or species), not only its single top hit. Run vsearch exhaustively (every hit at or
above the floor) and pipe its output straight in, so the large hit list never hits disk:

  vsearch --usearch_global out_lso_its2/queries.fasta --db out_lso_its2/refs.fasta --id 0.5 \\
      --maxaccepts 0 --maxrejects 0 --threads 32 --quiet \\
      --userfields query+target+id+qcov+tcov+ql+tl --userout /dev/stdout \\
  | python best_per_taxon.py --taxonomy out_lso_its2/taxonomy.tsv --rank g \\
      --queries out_lso_its2/queries.fasta --out lso_genus_best.tsv

On the real reference a capped search (--maxaccepts 2000) was no faster and missed about 9%
of the genera a query matched at 60% identity or more, which would make sets too small.

Hits must cover at least --min-cover percent (default 80) of the shorter of the two sequences.
vsearch's identity ignores terminal gaps, so an exhaustive search also finds short end-to-end
overlaps scored at high identity; for amplicons that keep conserved 5.8S or LSU flanks these
dominate the "best" hit. The coverage rule needs the userfields above.

Output: query_id, taxon, best_identity, one row per query and taxon it hit; a query with no
hit at all gets a single row with taxon "-" and identity 0, so every query is accounted for.
If a capped search is used anyway, pass its --maxaccepts and the run reports how many queries
reached it; for those, taxa beyond the cap are missing, which shrinks sets but leaves the
coverage guarantee intact, since calibration and testing use the same truncated lists.
"""
import argparse, sys
from collections import defaultdict

RANK_COL = {"s": "s", "g": "g", "f": "f", "o": "o"}

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--hits", default="-", help="vsearch userout (query, target, id); '-' = stdin")
    ap.add_argument("--taxonomy", required=True, help="the split's taxonomy.tsv (seq_id, p, c, o, f, g, s)")
    ap.add_argument("--rank", required=True, choices=list(RANK_COL), help="s, g, f or o")
    ap.add_argument("--queries", required=True, help="query FASTA, so queries without hits are kept")
    ap.add_argument("--maxaccepts", type=int, default=None, help="the vsearch --maxaccepts used")
    ap.add_argument("--min-cover", type=float, default=80.0,
                    help="minimum percent of the shorter sequence the alignment must cover (default 80)")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    taxon = {}
    with open(a.taxonomy) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        i_id, i_r = head.index("seq_id"), head.index(RANK_COL[a.rank])
        for line in fh:
            f = line.rstrip("\n").split("\t")
            taxon[f[i_id]] = f[i_r]
    queries = [l[1:].split()[0] for l in open(a.queries) if l.startswith(">")]

    best = defaultdict(dict)
    nhits = defaultdict(int)
    unknown_targets = nlines = short_cover = 0
    fh = sys.stdin if a.hits == "-" else open(a.hits)
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) < 3:
            continue
        nlines += 1
        if len(f) < 7:
            sys.exit("[best] hits need the fields query+target+id+qcov+tcov+ql+tl, to check alignment coverage")
        q, t = f[0], f[1]
        try:
            s = float(f[2]); qcov, tcov, ql, tl = float(f[3]), float(f[4]), int(f[5]), int(f[6])
        except ValueError:
            continue
        if (tcov if tl <= ql else qcov) < a.min_cover:
            short_cover += 1
            continue
        tx = taxon.get(t)
        if tx is None:
            unknown_targets += 1
            continue
        nhits[q] += 1
        if s > best[q].get(tx, -1.0):
            best[q][tx] = s
    if nlines == 0:
        sys.exit("[best] no hits were read at all: did the vsearch search run? (is vsearch on the PATH?)")
    with open(a.out, "w") as o:
        o.write("query_id\ttaxon\tbest_identity\n")
        for q in queries:
            if q in best:
                for tx, s in sorted(best[q].items(), key=lambda x: -x[1]):
                    o.write(f"{q}\t{tx}\t{s:.2f}\n")
            else:
                o.write(f"{q}\t-\t0\n")
    capped = sum(1 for q in queries if a.maxaccepts and nhits.get(q, 0) >= a.maxaccepts)
    sizes = sorted(len(best[q]) for q in queries if q in best) or [0]
    print(f"[best] {len(queries)} queries; {sum(1 for q in queries if q not in best)} without hits; "
          f"taxa per query median {sizes[len(sizes) // 2]}, max {sizes[-1]}; "
          f"{capped} reached maxaccepts; {short_cover} of {nlines} alignments dropped for covering under "
          f"{a.min_cover:g}% of the shorter sequence" + (f"; {unknown_targets} hits to targets missing from taxonomy"
                                            if unknown_targets else "") + f" -> {a.out}")

if __name__ == "__main__":
    main()
