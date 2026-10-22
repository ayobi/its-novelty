#!/usr/bin/env python3
"""Manuscript figure: three-level calls across ranks at one error budget alpha.

A  queries whose taxon is in the reference: species (species present), genus with the
   species missing (the calibration population), genus with the species present (same
   calls, applied to easier queries), family with the genus missing
B  queries whose taxon is absent from the reference, at each rank

Reads the _tiers.tsv tables written by conformal_tiers.py.

  python make_tiers_figure.py --species results/tiers_species --genus results/tiers_genus \\
      --family results/tiers_family --alpha 0.10 --out figures/fig_tiers
"""
import argparse, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

COLS = {"identified, right": "#009E73", "identified, wrong": "#E69F00",
        "unresolved, includes the true taxon": "#BBBBBB", "unresolved, misses the true taxon": "#666666",
        "unresolved, could be any taxon": "#EDEDED", "called absent": "#0072B2"}

def table(prefix):
    t = {}
    with open(prefix + "_tiers.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            t[(float(r["alpha"]), r["measure"])] = float(r["mean"])
    return t

def known(t, al, pre):
    g = lambda k: t[(al, f"{pre}_{k}")]
    return [g("id_ok"), g("id_wrong"), g("unres_ok") - g("vacuous"), g("unres_miss"), g("vacuous"), g("absent")], \
        g("unres_size_median")

def held(t, al):
    g = lambda k: t[(al, f"held_{k}")]
    return [0.0, g("identified"), 0.0, g("unresolved") - g("vacuous"), g("vacuous"), g("absent")], g("unres_size_median")

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--species", required=True); ap.add_argument("--genus", required=True)
    ap.add_argument("--family", required=True)
    ap.add_argument("--alpha", type=float, default=0.10)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    S, G, F = table(a.species), table(a.genus), table(a.family)
    al = a.alpha
    left = [("species\n(species in\nreference)", *known(S, al, "known")),
            ("genus\n(species\nmissing)", *known(G, al, "known")),
            ("genus\n(species in\nreference)", *known(G, al, "transfer")),
            ("family\n(genus\nmissing)", *known(F, al, "known"))]
    right = [("species", *held(S, al)), ("genus", *held(G, al)), ("family", *held(F, al))]

    plt.rcParams.update({"font.size": 8, "legend.fontsize": 7, "xtick.labelsize": 7, "ytick.labelsize": 7})
    fig, ax = plt.subplots(1, 2, figsize=(6.85, 3.3), gridspec_kw={"width_ratios": [4, 3]})
    labs = list(COLS)
    for axx, bars, title in ((ax[0], left, "A  taxon in the reference"), (ax[1], right, "B  taxon absent from the reference")):
        x = np.arange(len(bars))
        bottom = np.zeros(len(bars))
        for i, lab in enumerate(labs):
            v = np.array([b[1][i] for b in bars])
            axx.bar(x, v, bottom=bottom, color=COLS[lab], width=0.68, edgecolor="0.55", linewidth=0.3,
                    label=lab if axx is ax[0] else None)
            bottom += v
        for i, b in enumerate(bars):
            axx.text(i, 1.015, f"{b[2]:.0f}", ha="center", va="bottom", fontsize=6.5, color="0.35")
        axx.set_xticks(x); axx.set_xticklabels([b[0] for b in bars])
        axx.set_ylim(0, 1); axx.set_ylabel("share of queries")
        axx.set_title(title, loc="left", fontsize=9, fontweight="bold", pad=14)
    fig.legend(loc="lower center", ncol=3, frameon=False, handlelength=1.2)
    fig.tight_layout(rect=(0, 0.13, 1, 1))
    fig.savefig(a.out + ".pdf"); fig.savefig(a.out + ".png", dpi=300)
    print(f"[figure] wrote {a.out}.pdf and {a.out}.png (alpha = {al:g})")

if __name__ == "__main__":
    main()
