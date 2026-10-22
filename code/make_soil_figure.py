#!/usr/bin/env python3
"""Soil figure.

A  what the 5,521 ASVs were, by ASV and by reads (fungal ITS2 recovered by ITSxRust; the
   rest sorted by their 5.8S end), and how the fungal ITS2 ASVs were called at alpha = 0.10
B  ASVs called absent (alpha = 0.05) against the other fungal ITS2 ASVs: share recurring
   among GlobalFungi's unnamed variants, and share matching a named species-hypothesis
   representative (>= 97% identity over >= 90% of the GlobalFungi sequence)

  python make_soil_figure.py --kept asvs_its2_itsx.fasta --skipped asvs_itsx_skipped.fasta \\
      --abundance data/asv_abundance.tsv --calls refresh/results/soil_itsx_a0.10_apply_a0.1.tsv \\
      --table3 refresh/results/table3.txt --out fig_soil
"""
import argparse, csv, re
from collections import Counter
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ORIGIN = [("fungal ITS2", "#8C6D31"), ("plant-like 5.8S end", "#4DAF4A"),
          ("other 5.8S end", "#CC79A7"), ("no 5.8S end found", "#BBBBBB")]
CALLS = [("identified", "#009E73"), ("unresolved", "#D9D9D9"), ("absent", "#0072B2")]

def fasta(path):
    d, h = {}, None
    for l in open(path):
        if l.startswith(">"):
            h = l[1:].strip().split("|")[0].split()[0]; d[h] = []
        elif h:
            d[h].append(l.strip())
    return {k: "".join(v) for k, v in d.items()}

def origin(seq):
    if "AAGTTGCGCCC" in seq:
        return "plant-like 5.8S end"
    if re.search(r"TTGCGCC", seq):
        return "other 5.8S end"
    return "no 5.8S end found"

def table3(path):
    """(unnamed %, named %) at >= 97% for the flagged and non-flagged groups."""
    arm, out = None, {}
    for line in open(path):
        if "unnamed variant arm" in line: arm = "unnamed"
        elif "named representative arm" in line: arm = "named"
        f = line.split()
        if arm and f and f[0] in ("flagged", "non-flagged"):
            out[(arm, f[0])] = (int(f[1]), float(f[3].rstrip("%")))
    return out

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--kept", required=True); ap.add_argument("--skipped", required=True)
    ap.add_argument("--abundance", required=True); ap.add_argument("--calls", required=True)
    ap.add_argument("--table3", required=True); ap.add_argument("--out", required=True)
    a = ap.parse_args()
    kept, skipped = fasta(a.kept), fasta(a.skipped)
    reads = {r["asv_id"]: int(r["total_reads"]) for r in csv.DictReader(open(a.abundance), delimiter="\t")}
    cls = {x: "fungal ITS2" for x in kept}
    cls.update({x: origin(s) for x, s in skipped.items()})
    n_asv, n_reads = len(cls), sum(reads.get(x, 0) for x in cls)
    by_asv = Counter(cls.values())
    by_reads = Counter()
    for x, c in cls.items():
        by_reads[c] += reads.get(x, 0)
    calls = Counter(r["tier"] for r in csv.DictReader(open(a.calls), delimiter="\t"))
    n_calls = sum(calls.values())
    t3 = table3(a.table3)

    plt.rcParams.update({"font.size": 8, "legend.fontsize": 7, "xtick.labelsize": 7, "ytick.labelsize": 7})
    fig, ax = plt.subplots(1, 2, figsize=(6.85, 2.9), gridspec_kw={"width_ratios": [1.45, 1]})
    rows = [(f"ASVs\n({n_asv:,})", [(c, by_asv[c] / n_asv) for c, _ in ORIGIN], dict(ORIGIN)),
            (f"reads\n({n_reads / 1e6:.1f} M)", [(c, by_reads[c] / n_reads) for c, _ in ORIGIN], dict(ORIGIN)),
            (f"fungal ITS2\nASVs, \u03b1 = 0.10\n({n_calls:,})", [(c, calls[c] / n_calls) for c, _ in CALLS], dict(CALLS))]
    for i, (lab, parts, cols) in enumerate(rows):
        y = len(rows) - 1 - i - (0.35 if i == 2 else 0)
        left = 0.0
        for c, v in parts:
            ax[0].barh(y, v, left=left, color=cols[c], height=0.62, edgecolor="white", lw=0.6)
            if v >= 0.07:
                ax[0].text(left + v / 2, y, f"{v:.0%}", ha="center", va="center", fontsize=6.8,
                           color="white" if cols[c] in ("#8C6D31", "#0072B2", "#009E73", "#4DAF4A") else "black")
            left += v
        ax[0].text(-0.02, y, lab, ha="right", va="center", fontsize=7)
    ax[0].set_xlim(0, 1); ax[0].set_ylim(-0.95, 3.45)
    ax[0].set_yticks([]); ax[0].set_xticks([0, 0.25, 0.5, 0.75, 1]); ax[0].set_xticklabels(["0", "25%", "50%", "75%", "100%"])
    for s in ("left", "right", "top"):
        ax[0].spines[s].set_visible(False)
    h1 = [plt.Rectangle((0, 0), 1, 1, color=c) for _, c in ORIGIN]
    h2 = [plt.Rectangle((0, 0), 1, 1, color=c) for _, c in CALLS]
    leg = ax[0].legend(h1, [c for c, _ in ORIGIN], loc="upper left", bbox_to_anchor=(0.0, 1.06), ncol=2, frameon=False)
    ax[0].add_artist(leg)
    ax[0].legend(h2, [c for c, _ in CALLS], loc="lower left", bbox_to_anchor=(0.0, -0.36), ncol=3, frameon=False)
    ax[0].set_title("A", loc="left", fontweight="bold", pad=22)

    groups = ["recur among\nunnamed variants", "match a named\nspecies hypothesis"]
    absent = [t3[("unnamed", "flagged")][1], t3[("named", "flagged")][1]]
    other = [t3[("unnamed", "non-flagged")][1], t3[("named", "non-flagged")][1]]
    x = np.arange(2); w = 0.36
    b1 = ax[1].bar(x - w / 2, absent, w, color="#0072B2", label=f"called absent ({t3[('unnamed', 'flagged')][0]:,})")
    b2 = ax[1].bar(x + w / 2, other, w, color="#BBBBBB", label=f"other fungal ASVs ({t3[('unnamed', 'non-flagged')][0]:,})")
    for bars in (b1, b2):
        for b in bars:
            ax[1].text(b.get_x() + b.get_width() / 2, b.get_height() + 1.5, f"{b.get_height():.0f}%", ha="center",
                       va="bottom", fontsize=6.8)
    ax[1].set_xticks(x); ax[1].set_xticklabels(groups); ax[1].set_ylim(0, 112)
    ax[1].set_ylabel("share of ASVs (%)")
    ax[1].legend(loc="upper center", bbox_to_anchor=(0.5, -0.24), ncol=1, frameon=False)
    for s in ("right", "top"):
        ax[1].spines[s].set_visible(False)
    ax[1].set_title("B", loc="left", fontweight="bold", pad=22)
    fig.tight_layout(w_pad=3)
    fig.savefig(a.out + ".pdf", bbox_inches="tight"); fig.savefig(a.out + ".png", dpi=300, bbox_inches="tight")
    print(f"[figure] {n_asv} ASVs: " + ", ".join(f"{c} {by_asv[c]} ({by_reads[c] / n_reads:.1%} of reads)" for c, _ in ORIGIN)
          + f" -> {a.out}.pdf/.png")

if __name__ == "__main__":
    main()
