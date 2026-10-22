#!/usr/bin/env python3
"""Manuscript figure for the PROTAX-fungi comparison, from compare_protax_conformal.py output.

A  Operating points. ROC curves (share of known-genus queries called absent against
   share of held-out-genus queries called absent) for best-hit identity and PROTAX's
   absence probability, with where each decision rule lands: PROTAX's own probability
   cutoffs, and the conformal call at alpha = 0.05 and 0.10 on either score.
B  Lists. Share of absence calls that are wrong, against the share of held-out-genus
   queries in the sample: PROTAX's own expected-error selection (stated: 5% wrong),
   PROTAX >= 0.95, and the conformal call at alpha = 0.05 on either score.

  python make_protax_figure.py --results results/protax_vs_conformal --out figures/fig_protax
"""
import argparse, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BLUE, RED = "#2166ac", "#b2182b"

def tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))

def roc(known, held):
    """False-absence rate and detection for every cutoff, most novel first."""
    t = np.unique(np.r_[known, held])[::-1]
    fa = np.array([np.mean(known >= x) for x in t])
    de = np.array([np.mean(held >= x) for x in t])
    return np.r_[0, fa, 1], np.r_[0, de, 1]

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", required=True, help="--out-prefix used for compare_protax_conformal.py")
    ap.add_argument("--out", required=True, help="output prefix (.pdf and .png are written)")
    ap.add_argument("--xmax", type=float, default=0.30, help="x range of panel A")
    a = ap.parse_args()

    pq = tsv(a.results + "_perquery.tsv")
    held = np.array([r["population"] == "novel" for r in pq])
    ident = np.array([float(r["identity"]) for r in pq])
    pabs = np.array([float(r["p_absent_protax"]) for r in pq])
    conf = {(r["score"], float(r["alpha"])): r for r in tsv(a.results + "_conformal.tsv")}
    thr = {float(r["threshold"]): r for r in tsv(a.results + "_protax_thresholds.tsv")}
    mix = tsv(a.results + "_mixtures.tsv")

    plt.rcParams.update({"font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
                         "legend.fontsize": 7, "xtick.labelsize": 7, "ytick.labelsize": 7})
    fig, ax = plt.subplots(1, 2, figsize=(6.85, 3.0))

    # ---- A: operating points
    x, y = roc(-ident[~held], -ident[held]); ax[0].plot(x, y, color=BLUE, lw=1.2, label="best-hit identity")
    x, y = roc(pabs[~held], pabs[held]); ax[0].plot(x, y, color=RED, lw=1.2, label="PROTAX P(absent)")
    for al in (0.05, 0.10):
        ax[0].axvline(al, color="0.75", lw=0.6, ls=":", zorder=0)
        for score, col, mk in (("identity", BLUE, "o"), ("PROTAX", RED, "s")):
            r = conf[(score, al)]
            ax[0].plot(float(r["false_absence_mean"]), float(r["detection_mean"]), mk, color=col, ms=5,
                       mec="black", mew=0.6, zorder=3)
    for t in (0.95, 0.99):
        r = thr[t]
        fa, de = float(r["false_absence_known"]), float(r["detection_heldout"])
        ax[0].plot(fa, de, "^", mfc="white", mec=RED, mew=1.2, ms=6, zorder=3)
        ax[0].annotate(f"P \u2265 {t:g}", (fa, de), xytext=(5, -10), textcoords="offset points", fontsize=7, color=RED)
    ax[0].plot([], [], "o", color=BLUE, mec="black", mew=0.6, ms=5, label="conformal call on identity")
    ax[0].plot([], [], "s", color=RED, mec="black", mew=0.6, ms=5, label="conformal call on PROTAX")
    ax[0].plot([], [], "^", mfc="white", mec=RED, mew=1.2, ms=6, label="PROTAX probability cutoff")
    for al in (0.05, 0.10):
        ax[0].text(al, 0.785, f"\u03b1 = {al:g}", fontsize=6.5, color="0.45", ha="center", va="top",
                   bbox=dict(facecolor="white", edgecolor="none", pad=0.5))
    ax[0].set_xlim(0, a.xmax); ax[0].set_ylim(0, 0.8)
    ax[0].set_xlabel("share of known-genus queries called absent")
    ax[0].set_ylabel("share of held-out-genus queries called absent")
    ax[0].legend(loc="lower right", frameon=False)
    ax[0].set_title("A", loc="left", fontweight="bold")

    # ---- B: lists
    lines = [("PROTAX expected FDR<=0.05", RED, ":", "^", "PROTAX's own list, stated 5% wrong"),
             ("PROTAX P>=0.95", RED, "-", "^", "PROTAX P(absent) \u2265 0.95"),
             ("conformal PROTAX a=0.05", RED, "--", "s", "conformal call on PROTAX, \u03b1 = 0.05"),
             ("conformal identity a=0.05", BLUE, "-", "o", "conformal call on identity, \u03b1 = 0.05")]
    for m, col, ls, mk, lab in lines:
        r = sorted((x for x in mix if x["method"] == m), key=lambda x: float(x["pi"]))
        pi = [float(x["pi"]) for x in r]
        ax[1].errorbar(pi, [float(x["fdr_mean"]) for x in r], yerr=[1.96 * float(x["fdr_se"]) for x in r],
                       color=col, ls=ls, marker=mk, ms=3.5, lw=1.1, capsize=0, label=lab,
                       mfc="white" if mk == "^" else col)
    ax[1].axhline(0.05, color="0.5", lw=0.7, ls=":")
    ax[1].text(0.27, 0.065, "5%", fontsize=7, color="0.4")
    ticks = sorted({float(x["pi"]) for x in mix})
    ax[1].set_xscale("log"); ax[1].set_xticks(ticks); ax[1].set_xticklabels([f"{t:g}" for t in ticks]); ax[1].minorticks_off()
    ax[1].set_ylim(0, 1)
    ax[1].set_xlabel("share of held-out-genus queries in the sample")
    ax[1].set_ylabel("share of absence calls that are wrong")
    ax[1].legend(loc="lower left", bbox_to_anchor=(0.0, 0.1), frameon=False)
    ax[1].set_title("B", loc="left", fontweight="bold")

    fig.tight_layout()
    fig.savefig(a.out + ".pdf"); fig.savefig(a.out + ".png", dpi=300)
    print(f"[figure] wrote {a.out}.pdf and {a.out}.png")

if __name__ == "__main__":
    main()
