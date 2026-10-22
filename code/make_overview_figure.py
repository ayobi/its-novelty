#!/usr/bin/env python3
"""Figure 1: how the calls work.

A  best-hit identity of known-genus queries (held-out species: their genus is in the
   reference) and held-out-genus queries, with the cutoff the absence call derives at
   alpha = 0.05 and the familiar 97% for comparison (real data)
B  three-level calls on three illustrative queries: candidate genera within d points of
   the best hit form the set; one genus = identified, several = unresolved, a best hit
   below the absence cutoff = absent (schematic; d and the cutoff are the paper's values)

  python make_overview_figure.py --perquery refresh/results/protax_vs_conformal_perquery.tsv \\
      --conformal refresh/results/protax_vs_conformal_conformal.tsv --out fig_overview
"""
import argparse, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

KNOWN, HELD = "#0072B2", "#E69F00"

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--perquery", required=True, help="compare_protax_conformal.py _perquery.tsv")
    ap.add_argument("--conformal", required=True, help="compare_protax_conformal.py _conformal.tsv")
    ap.add_argument("--width", type=float, default=15.2, help="set width d at alpha = 0.10 (soil calibration)")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    rows = list(csv.DictReader(open(a.perquery), delimiter="\t"))
    k = np.array([float(r["identity"]) for r in rows if r["population"] == "known"])
    h = np.array([float(r["identity"]) for r in rows if r["population"] == "novel"])
    conf = {(r["score"], float(r["alpha"])): r for r in csv.DictReader(open(a.conformal), delimiter="\t")}
    c05 = float(conf[("identity", 0.05)]["median_cutoff"])
    fa05 = float(conf[("identity", 0.05)]["false_absence_mean"]); de05 = float(conf[("identity", 0.05)]["detection_mean"])

    plt.rcParams.update({"font.size": 8, "legend.fontsize": 7, "xtick.labelsize": 7, "ytick.labelsize": 7})
    fig, ax = plt.subplots(1, 2, figsize=(6.85, 3.1), gridspec_kw={"width_ratios": [1.15, 1]})

    # ---- A: real identity distributions
    bins = np.arange(50, 101, 2)
    for x, col, lab in ((k, KNOWN, f"genus present, species removed (n = {len(k):,})"),
                        (h, HELD, f"genus removed (n = {len(h):,})")):
        hit = x[x > 0]
        w = np.full(len(hit), 1 / len(x))
        ax[0].hist(hit, bins=bins, weights=w, histtype="stepfilled", alpha=0.35, color=col)
        ax[0].hist(hit, bins=bins, weights=w, histtype="step", lw=1.1, color=col, label=lab)
        ax[0].bar(46, np.mean(x == 0), width=2.6, color=col, alpha=0.35, edgecolor=col, lw=1.1)
    top = ax[0].get_ylim()[1]
    ax[0].set_ylim(0, top * 1.5)
    ax[0].axvline(c05, color="black", lw=1.1, ymax=0.72)
    ax[0].axvline(97, color="0.45", lw=1, ls=":", ymax=0.72)
    tr = ax[0].get_xaxis_transform()
    ax[0].text(c05, 0.735, f"cutoff at \u03b1 = 0.05\n({c05:.1f}%)", transform=tr, ha="center", va="bottom", fontsize=6.4)
    ax[0].text(97, 0.735, "97%", transform=tr, ha="center", va="bottom", fontsize=6.4, color="0.35")
    ax[0].text(48.6, top * 0.98, f"called absent:\n{fa05:.0%} of present\n{de05:.0%} of removed",
               ha="left", va="top", fontsize=6.1)
    ax[0].set_xticks([46, 50, 60, 70, 80, 90, 100]); ax[0].set_xticklabels(["no\nhit", "50", "60", "70", "80", "90", "100"])
    ax[0].set_xlim(43.5, 101)
    ax[0].set_xlabel("best-hit identity to the reference (%)"); ax[0].set_ylabel("share of queries")
    ax[0].legend(loc="upper left", bbox_to_anchor=(0.0, 1.0), frameon=False, handlelength=1.4, fontsize=6.6)
    ax[0].set_title("A", loc="left", fontweight="bold")

    # ---- B: schematic of the three-level calls
    d, cut = a.width, c05
    queries = [("query 1", [("A", 98)], [("B", 81), ("C", 76)], "identified: genus A"),
               ("query 2", [("A", 91), ("B", 88), ("C", 85)], [("D", 72)], "unresolved: A, B or C"),
               ("query 3", [], [("A", 64), ("B", 58)], "absent from the reference")]
    for i, (name, inside, outside, verdict) in enumerate(queries):
        y = len(queries) - 1 - i
        best = max(s for _, s in inside + outside)
        ax[1].plot([50, 100], [y, y], color="0.85", lw=0.8, zorder=0)
        if inside:
            ax[1].add_patch(Rectangle((best - d, y - 0.22), d, 0.44, color="#BBBBBB", alpha=0.5, lw=0, zorder=1))
        for g, s in inside:
            ax[1].plot(s, y, "o", ms=7, color="#009E73", zorder=3); ax[1].text(s, y + 0.26, g, ha="center", fontsize=7)
        for g, s in outside:
            ax[1].plot(s, y, "o", ms=7, mfc="white", mec="0.4", zorder=3); ax[1].text(s, y + 0.26, g, ha="center", fontsize=7, color="0.4")
        ax[1].text(49, y, name, ha="right", va="center", fontsize=7)
        ax[1].text(101, y, verdict, ha="left", va="center", fontsize=7,
                   color="#0072B2" if not inside else ("#009E73" if len(inside) == 1 else "0.3"))
    ax[1].axvline(cut, color="black", lw=1.1)
    ax[1].text(cut - 0.8, 2.62, "absence\ncutoff", ha="right", va="top", fontsize=6.5)
    ax[1].annotate("", xy=(98 - d, 2.5), xytext=(98, 2.5), arrowprops=dict(arrowstyle="<->", lw=0.7, color="0.3"))
    ax[1].text(98 - d / 2, 2.62, f"width d = {d:.0f} points", ha="center", va="bottom", fontsize=6.5, color="0.3")
    ax[1].set_xlim(40, 100); ax[1].set_ylim(-0.6, 3.05)
    ax[1].set_yticks([]); ax[1].set_xlabel("identity to each reference genus (%)")
    for s in ("left", "right", "top"):
        ax[1].spines[s].set_visible(False)
    ax[1].set_title("B", loc="left", fontweight="bold", pad=14)
    fig.tight_layout(w_pad=2.5)
    fig.subplots_adjust(right=0.80)
    fig.savefig(a.out + ".pdf"); fig.savefig(a.out + ".png", dpi=300)
    print(f"[figure] wrote {a.out}.pdf/.png (cutoff {c05:.1f}%, width {d:.1f})")

if __name__ == "__main__":
    main()
