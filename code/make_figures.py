#!/usr/bin/env python3
"""Regenerate manuscript figures from disk. Outputs fig1-4 .pdf.
Run from ~/Downloads. Requires matplotlib, numpy."""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 11, "svg.fonttype": "none",
                     "axes.spines.top": False, "axes.spines.right": False})

# ---- Fig 1: leakage / four-method recovery -------------------------------
def fig1():
    labels = ["Embed", "Align", "CNN", "LCA", r"Embed$_{\rm LGO}$"]
    vals   = [69.0, 41.0, 37.0, 38.9, 23.0]   # TODO: replace 23.0 w/ regenerated
    fig, ax = plt.subplots(figsize=(6.5, 4))
    bars = ax.bar(range(4), vals[:4], color="0.55", width=0.6)
    ax.bar(4, vals[4], color="white", edgecolor="0.3", hatch="////", width=0.6)
    ax.axhspan(37, 41, color="0.9", zorder=0)  # honest band
    for i, v in enumerate(vals):
        ax.text(i, v+1, f"{v:.1f}", ha="center", fontsize=10)
    ax.set_xticks(range(5)); ax.set_xticklabels(labels)
    ax.set_ylabel("Family-level recovery (%)"); ax.set_ylim(0, 78)
    fig.tight_layout(); fig.savefig("fig1_leakage.pdf")

# ---- Fig 2: conformal calibration ----------------------------------------
def fig2():
    a   = [0.01, 0.05, 0.10, 0.20]
    ff  = [0.0070, 0.0426, 0.1056, 0.1972]
    det = [0.0648, 0.1922, 0.3326, 0.5092]
    fig, ax = plt.subplots(figsize=(6.5, 4))
    ax.plot([0,0.2],[0,0.2], "--", color="0.6", label=r"Ideal ($y=x$)")
    ax.plot(a, ff, "o-", color="0.2", label="Empirical false-flag")
    ax.plot(a, det, "s-", color="tab:blue", label="Novel detected")
    ax.set_xlabel(r"Nominal $\alpha$"); ax.set_ylabel("Rate")
    ax.set_xlim(0,0.22); ax.set_ylim(0,0.55); ax.legend(loc="upper left")
    ax.grid(alpha=0.3); fig.tight_layout(); fig.savefig("fig2_calibration.pdf")

# ---- Fig 3: known/novel/outgroup gradient --------------------------------
def fig3():
    labels = ["Known\n(LSO)", "Novel\n(LGO)", "Outgroup\n(SSU)"]
    vals   = [4.3, 19.2, 99.9]
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.bar(range(3), vals, color=["0.7","0.5","0.3"], width=0.6)
    for i, v in enumerate(vals): ax.text(i, v+1.5, f"{v}", ha="center")
    ax.set_xticks(range(3)); ax.set_xticklabels(labels)
    ax.set_ylabel(r"Flagged novel at $\alpha=0.05$ (%)"); ax.set_ylim(0,108)
    fig.tight_layout(); fig.savefig("fig3_gradient.pdf")

# ---- Fig 4: soil dark-matter fractions -----------------------------------
def fig4():
    labels = [r"$\alpha=0.05$", r"$\alpha=0.10$", r"$\alpha=0.20$"]
    vals   = [20.7, 26.1, 38.8]
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.bar(range(3), vals, color="0.55", width=0.6)
    for i, v in enumerate(vals): ax.text(i, v+0.7, f"{v}", ha="center")
    ax.set_xticks(range(3)); ax.set_xticklabels(labels)
    ax.set_ylabel("ASVs flagged novel (%)"); ax.set_ylim(0,44)
    fig.tight_layout(); fig.savefig("fig4_darkmatter.pdf")

def fig5():
    fig, ax = plt.subplots(figsize=(5.5, 4))
    ax.bar([0,1], [396, 704], color=["tab:green","0.6"], width=0.6)
    for i,v in zip([0,1],[396,704]): ax.text(i, v+15, str(v), ha="center")
    ax.set_xticks([0,1])
    ax.set_xticklabels(["Recurs globally\n(>=97% id, >=90% cov)","No global\nmatch"])
    ax.set_ylabel("Novel ASVs (of 1,100)"); ax.set_ylim(0,780)
    ax.set_title("36% of flagged dark matter recurs in global soils", fontsize=10)
    fig.tight_layout(); fig.savefig("fig5_recurrence.pdf")

if __name__ == "__main__":
    fig1(); fig2(); fig3(); fig4(); fig5()
    print("wrote fig1_leakage.pdf, fig2_calibration.pdf, "
          "fig3_gradient.pdf, fig4_darkmatter.pdf")