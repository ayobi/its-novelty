#!/usr/bin/env python3
"""
make_flagged_fasta.py -- write the flagged-novel ASV set, defined by the conformal
p-value rather than by the identity cutoff the procedure derives.

  python code/make_flagged_fasta.py --alpha 0.05

Those two definitions are not interchangeable. An earlier version of this analysis built
the flagged set by thresholding identity below the derived 64.8% cutoff, which omitted 43
ASVs whose p-values were at or below 0.05 while their identities fell between 64.8 and
65.2. The p-value is the actual procedure; this script exists so the deposited FASTA and
the deposited verdict table cannot disagree again.
"""
import argparse

ap = argparse.ArgumentParser()
ap.add_argument("--verdicts", default="data/frey_asv_pvalues.tsv")
ap.add_argument("--asvs", default="data/asvs_its2.fasta")
ap.add_argument("--alpha", type=float, default=0.05,
                help="0.05, 0.10 or 0.20; selects the matching novel_a* column")
ap.add_argument("--out", default="data/flagged_novel_asvs.fa")
a = ap.parse_args()

col = "novel_a%02d" % round(a.alpha * 100)
with open(a.verdicts) as fh:
    header = fh.readline().rstrip("\n").split("\t")
    if col not in header:
        raise SystemExit("no column %r in %s; have %s" % (col, a.verdicts, header))
    i_id, i_flag = header.index("asv_id"), header.index(col)
    flagged = {r.split("\t")[i_id] for r in fh
               if r.strip() and r.split("\t")[i_flag].strip() == "1"}
print("flagged at alpha=%.2f: %d ASVs" % (a.alpha, len(flagged)))

written, keep = 0, False
with open(a.asvs) as fh, open(a.out, "w") as out:
    for line in fh:
        if line.startswith(">"):
            keep = line[1:].split()[0] in flagged
            written += keep
        if keep:
            out.write(line)
print("wrote %d of %d to %s" % (written, len(flagged), a.out))
if written != len(flagged):
    raise SystemExit("MISMATCH: %d flagged ids absent from the FASTA"
                     % (len(flagged) - written))
