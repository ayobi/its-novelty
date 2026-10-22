#!/usr/bin/env python3
"""PROTAX-fungi against the conformal absence call, on the same held-out queries.

Five measurements, all at genus rank, known = held-out-species (LSO) queries whose genus
is in the reference, novel = held-out-genus (LGO) queries whose genus is not:

  1. discrimination  AUROC of PROTAX's absence probability and of best-hit identity
  2. PROTAX as used  false-absence rate on known queries and detection on held-out-genus
                     queries when absence is called at PROTAX probability >= t
  3. conformal       false-absence rate and detection of the conformal call on each score,
                     averaged over random calibration splits, with the cutoff it implies
  4. mixtures        test sets with a chosen share pi of held-out-genus queries; realised
                     share of wrong absence calls (FDR) and detection for PROTAX >= t,
                     PROTAX's own expected-FDR selection, the per-query conformal call, and
                     Benjamini-Hochberg on conformal p-values (Bates et al. 2023)
  5. per-query table every score, for reuse

Queries with no identity hit at the search floor are kept at --nohit-identity (default 0,
the most novel value) in calibration and evaluation alike. conformal_novelty.py drops
them from both; --hits-only reproduces that behaviour for comparison with Table 1.

  python compare_protax_conformal.py \
      --known-queries out_lso_its2/queries.fasta --novel-queries out_lgo_its2/queries.fasta \
      --known-hits lso_its2.hits --novel-hits lgo_its2.hits \
      --known-protax protax_lso/scores_genus.tsv --novel-protax protax_lgo/scores_genus.tsv \
      --out-prefix results/protax_vs_conformal
"""
import argparse, os, sys
import numpy as np

# ---------------------------------------------------------------- io

def fasta_ids(path):
    with open(path) as fh:
        return [l[1:].split()[0] for l in fh if l.startswith(">")]

def best_identity(path):
    best = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            try:
                s = float(f[2])
            except ValueError:
                continue
            if f[0] not in best or s > best[f[0]]:
                best[f[0]] = s
    return best

def protax_table(path):
    out = {}
    with open(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        i, j = head.index("query_id"), head.index("p_absent")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            out[f[i]] = float(f[j])
    return out

# ---------------------------------------------------------------- statistics

def pvalues(cal, test):
    """Conformal p-values; higher score = more novel. (1 + #cal >= test) / (n + 1)."""
    c = np.sort(cal)
    return (1 + len(c) - np.searchsorted(c, test, side="left")) / (len(c) + 1)

def auroc(labels, scores):
    y = np.asarray(labels, bool); s = np.asarray(scores, float)
    order = np.argsort(s, kind="mergesort")
    _, inv, cnt = np.unique(s[order], return_inverse=True, return_counts=True)
    end = np.cumsum(cnt); avg = (end - cnt + 1 + end) / 2.0
    r = np.empty(len(s)); r[order] = avg[inv]
    npos, nneg = y.sum(), (~y).sum()
    return float("nan") if npos == 0 or nneg == 0 else (r[y].sum() - npos * (npos + 1) / 2) / (npos * nneg)

def bh(p, q):
    m = len(p)
    if m == 0:
        return np.zeros(0, bool)
    ps = np.sort(p)
    ok = np.nonzero(ps <= q * np.arange(1, m + 1) / m)[0]
    return p <= ps[ok.max()] if len(ok) else np.zeros(m, bool)

def expected_fdr_select(pabs, q):
    """PROTAX's own shortlist: most-probable-absent first, while mean(1 - p) <= q."""
    order = np.argsort(-pabs, kind="mergesort")
    run = np.cumsum(1 - pabs[order]) / np.arange(1, len(pabs) + 1)
    ok = np.nonzero(run <= q)[0]
    sel = np.zeros(len(pabs), bool)
    if len(ok):
        sel[order[:ok.max() + 1]] = True
    return sel

def cutoff_at(cal, alpha):
    """Score above which a test query is flagged at level alpha (flag iff score > cutoff)."""
    j = int(np.floor(alpha * (len(cal) + 1))) - 1
    if j < 0:
        return float("inf")
    return np.sort(cal)[::-1][j]

# ---------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for side in ("known", "novel"):
        ap.add_argument(f"--{side}-queries", required=True, help=f"{side} queries.fasta")
        ap.add_argument(f"--{side}-hits", required=True, help="query, target, % identity")
        ap.add_argument(f"--{side}-protax", required=True, help="protax_scores.py output")
    ap.add_argument("--alphas", default="0.01,0.05,0.10,0.20")
    ap.add_argument("--thresholds", default="0.5,0.9,0.95,0.99")
    ap.add_argument("--fdr", default="0.05,0.10")
    ap.add_argument("--mixes", default="0.01,0.02,0.05,0.10,0.20,0.50")
    ap.add_argument("--mix-size", type=int, default=1000)
    ap.add_argument("--draws", type=int, default=500)
    ap.add_argument("--splits", type=int, default=200)
    ap.add_argument("--calib-frac", type=float, default=0.5)
    ap.add_argument("--nohit-identity", type=float, default=0.0)
    ap.add_argument("--hits-only", action="store_true", help="drop no-hit queries (as Table 1 did)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out-prefix", required=True)
    a = ap.parse_args()
    alphas = [float(x) for x in a.alphas.split(",")]
    ths = [float(x) for x in a.thresholds.split(",")]
    qs = [float(x) for x in a.fdr.split(",")]
    mixes = [float(x) for x in a.mixes.split(",")]
    os.makedirs(os.path.dirname(a.out_prefix) or ".", exist_ok=True)
    lines = []
    def say(s=""):
        print(s); lines.append(s)

    # ---- assemble the two populations
    pop = {}
    for side in ("known", "novel"):
        ids = fasta_ids(getattr(a, f"{side}_queries"))
        hits = best_identity(getattr(a, f"{side}_hits"))
        prot = protax_table(getattr(a, f"{side}_protax"))
        miss = [q for q in ids if q not in prot]
        if miss:
            sys.exit(f"[compare] {len(miss)} {side} queries missing from the PROTAX table, e.g. {miss[:3]}")
        has = np.array([q in hits for q in ids])
        ident = np.array([hits.get(q, a.nohit_identity) for q in ids], float)
        pab = np.array([prot[q] for q in ids], float)
        if a.hits_only:
            ids = [q for q, h in zip(ids, has) if h]; ident, pab, has = ident[has], pab[has], has[has]
        pop[side] = dict(ids=ids, ident=ident, pabs=pab, has=has)
        say(f"[{side}] {len(ids)} queries, {int((~has).sum())} without an identity hit "
            f"(kept at identity {a.nohit_identity:g})" if not a.hits_only else
            f"[{side}] {len(ids)} queries with an identity hit")
    K, N = pop["known"], pop["novel"]
    scores = {"identity": (-K["ident"], -N["ident"]), "PROTAX": (K["pabs"], N["pabs"])}
    say()

    # ---- 1. discrimination
    say("1. Discrimination (AUROC, held-out genus vs known genus)")
    y = np.r_[np.zeros(len(K["ids"])), np.ones(len(N["ids"]))]
    for name, (sk, sn) in scores.items():
        say(f"   {name:9s} {auroc(y, np.r_[sk, sn]):.4f}")
    say()

    # ---- 2. PROTAX as people use it
    say("2. PROTAX as used: call absent at P(absent) >= t")
    say(f"   {'t':>5} {'false-absence (known)':>22} {'detection (held-out)':>21}")
    rows2 = []
    for t in ths:
        fa, de = float(np.mean(K["pabs"] >= t)), float(np.mean(N["pabs"] >= t))
        rows2.append((t, fa, de)); say(f"   {t:5.2f} {fa:22.4f} {de:21.4f}")
    say()

    # ---- 3. conformal, averaged over calibration splits
    rng = np.random.default_rng(a.seed)
    nk = len(K["ids"]); ncal = int(nk * a.calib_frac)
    res3 = {(s, al): ([], [], []) for s in scores for al in alphas}
    for _ in range(a.splits):
        perm = rng.permutation(nk); cal, tst = perm[:ncal], perm[ncal:]
        for s, (sk, sn) in scores.items():
            pt, pn = pvalues(sk[cal], sk[tst]), pvalues(sk[cal], sn)
            for al in alphas:
                ff, de, cu = res3[(s, al)]
                ff.append(np.mean(pt <= al)); de.append(np.mean(pn <= al)); cu.append(cutoff_at(sk[cal], al))
    say(f"3. Conformal call, {a.splits} random calibration splits ({ncal} calibration knowns each)")
    say(f"   {'score':9s} {'alpha':>5} {'false-absence':>14} {'detection':>16}   implied cutoff (median)")
    rows3 = []
    for s in scores:
        for al in alphas:
            ff, de, cu = (np.array(v, float) for v in res3[(s, al)])
            med = float(np.median(cu))                     # flag iff score > med (median split)
            cut = -med if s == "identity" else med
            if s == "identity":
                lab = "nothing can be flagged" if not np.isfinite(med) or cut <= a.nohit_identity \
                    else f"identity < {cut:.1f}%"
            else:
                lab = "nothing can be flagged" if not np.isfinite(med) or med >= 1 else f"P > {cut:.3f}"
            rows3.append((s, al, ff.mean(), ff.std(), de.mean(), de.std(), cut))
            say(f"   {s:9s} {al:5.2f} {ff.mean():8.4f} ±{ff.std():.4f} {de.mean():9.4f} ±{de.std():.4f}   {lab}")
    say()

    # ---- 4. mixtures
    methods = [f"PROTAX P>={t:g}" for t in ths] + [f"PROTAX expected FDR<={q:g}" for q in qs] \
            + [f"conformal identity a={q:g}" for q in qs] + [f"conformal PROTAX a={q:g}" for q in qs] \
            + [f"BH identity q={q:g}" for q in qs] + [f"BH PROTAX q={q:g}" for q in qs]
    fdp = {(pi, m): [] for pi in mixes for m in methods}
    pw = {(pi, m): [] for pi in mixes for m in methods}
    nsel = {(pi, m): [] for pi in mixes for m in methods}
    for pi in mixes:
        for _ in range(a.draws):
            perm = rng.permutation(nk); cal, pool = perm[:ncal], perm[ncal:]
            n_nov = max(1, int(round(pi * a.mix_size))); n_kn = a.mix_size - n_nov
            if n_kn > len(pool):
                n_kn = len(pool); n_nov = max(1, int(round(n_kn * pi / (1 - pi))))
            n_nov = min(n_nov, len(N["ids"]))
            tk = rng.choice(pool, n_kn, replace=False); tn = rng.choice(len(N["ids"]), n_nov, replace=False)
            is_nov = np.r_[np.zeros(n_kn, bool), np.ones(n_nov, bool)]
            pab = np.r_[K["pabs"][tk], N["pabs"][tn]]
            pv = {s: pvalues(sk[cal], np.r_[sk[tk], sn[tn]]) for s, (sk, sn) in scores.items()}
            sel = {}
            for t in ths:
                sel[f"PROTAX P>={t:g}"] = pab >= t
            for q in qs:
                sel[f"PROTAX expected FDR<={q:g}"] = expected_fdr_select(pab, q)
                sel[f"conformal identity a={q:g}"] = pv["identity"] <= q
                sel[f"conformal PROTAX a={q:g}"] = pv["PROTAX"] <= q
                sel[f"BH identity q={q:g}"] = bh(pv["identity"], q)
                sel[f"BH PROTAX q={q:g}"] = bh(pv["PROTAX"], q)
            for m, s in sel.items():
                ns = int(s.sum())
                fdp[(pi, m)].append((s & ~is_nov).sum() / ns if ns else 0.0)
                pw[(pi, m)].append((s & is_nov).sum() / n_nov)
                nsel[(pi, m)].append(ns)
    say(f"4. Mixtures: {a.draws} draws of {a.mix_size} queries per share of held-out-genus queries (pi)")
    say("   cells: realised FDR (share of absence calls that are wrong) / detection")
    w = max(len(m) for m in methods)
    say("   " + " " * w + "".join(f"{'pi=' + format(pi, 'g'):>17}" for pi in mixes))
    rows4 = []
    for m in methods:
        cells = []
        for pi in mixes:
            f, p, n = np.array(fdp[(pi, m)]), np.array(pw[(pi, m)]), np.array(nsel[(pi, m)])
            rows4.append((m, pi, f.mean(), f.std() / np.sqrt(len(f)), p.mean(), n.mean()))
            cells.append(f"{f.mean():7.3f}/{p.mean():.3f}".rjust(17))
        say("   " + m.ljust(w) + "".join(cells))
    say()

    # ---- outputs
    P = a.out_prefix
    with open(P + "_summary.txt", "w") as fh:
        fh.write("\n".join(lines) + "\n")
    with open(P + "_protax_thresholds.tsv", "w") as fh:
        fh.write("threshold\tfalse_absence_known\tdetection_heldout\n")
        fh.writelines(f"{t}\t{fa:.5f}\t{de:.5f}\n" for t, fa, de in rows2)
    with open(P + "_conformal.tsv", "w") as fh:
        fh.write("score\talpha\tfalse_absence_mean\tfalse_absence_sd\tdetection_mean\tdetection_sd\tmedian_cutoff\n")
        fh.writelines(f"{s}\t{al}\t{a1:.5f}\t{a2:.5f}\t{a3:.5f}\t{a4:.5f}\t{a5:.5g}\n" for s, al, a1, a2, a3, a4, a5 in rows3)
    with open(P + "_mixtures.tsv", "w") as fh:
        fh.write("method\tpi\tfdr_mean\tfdr_se\tdetection_mean\tcalls_mean\n")
        fh.writelines(f"{m}\t{pi}\t{f:.5f}\t{se:.5f}\t{p:.5f}\t{n:.1f}\n" for m, pi, f, se, p, n in rows4)
    with open(P + "_perquery.tsv", "w") as fh:
        fh.write("query_id\tpopulation\tidentity\tidentity_hit\tp_absent_protax\n")
        for side in ("known", "novel"):
            d = pop[side]
            for q, i, h, p in zip(d["ids"], d["ident"], d["has"], d["pabs"]):
                fh.write(f"{q}\t{side}\t{i:.2f}\t{int(h)}\t{p:.6g}\n")
    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        q0 = qs[0]
        show = [(f"PROTAX P>={t:g}", "#b2182b", ls) for t, ls in zip([x for x in ths if x >= 0.9][:2], ["-", "--"])] + [
            (f"PROTAX expected FDR<={q0:g}", "#ef8a62", ":"),
            (f"conformal identity a={q0:g}", "#67a9cf", "--"),
            (f"BH identity q={q0:g}", "#2166ac", "-"),
            (f"BH PROTAX q={q0:g}", "#1b7837", "-")]
        fig, ax = plt.subplots(1, 2, figsize=(10, 4))
        for m, col, ls in show:
            r = [x for x in rows4 if x[0] == m]
            ax[0].errorbar([x[1] for x in r], [x[2] for x in r], yerr=[1.96 * x[3] for x in r],
                           color=col, ls=ls, marker="o", ms=3, label=m)
            ax[1].plot([x[1] for x in r], [x[4] for x in r], color=col, ls=ls, marker="o", ms=3)
        ax[0].axhline(q0, color="0.5", lw=0.8, ls=":")
        for x in ax:
            x.set_xscale("log"); x.set_xlabel("share of held-out-genus queries in the sample")
            x.set_xticks(mixes); x.set_xticklabels([f"{m:g}" for m in mixes]); x.minorticks_off()
        ax[0].set_ylabel("share of absence calls that are wrong (FDR)")
        ax[1].set_ylim(-0.02, 1.02)
        ax[1].set_ylabel("detection of held-out-genus queries")
        ax[0].legend(fontsize=7, frameon=False)
        fig.tight_layout(); fig.savefig(P + "_mixtures.pdf"); fig.savefig(P + "_mixtures.png", dpi=150)
    except ImportError:
        print("[compare] matplotlib not available; skipped the figure")
    print(f"[compare] wrote {P}_summary.txt, _protax_thresholds.tsv, _conformal.tsv, "
          f"_mixtures.tsv, _perquery.tsv, _mixtures.pdf/.png")

if __name__ == "__main__":
    main()
