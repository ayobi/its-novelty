#!/usr/bin/env python3
"""Three-level calls at one rank: identified, absent from the reference, or unresolved.

Two conformal steps share one error budget alpha (split by --absence-share):

  absent      the absence test on best-hit identity at alpha1: a query is called absent when
              its best hit is lower than all but a fraction alpha1 of known queries' best hits
  otherwise   the set of reference taxa whose best identity is within d points of the
              query's best hit; one taxon = identified, two or more = unresolved. d is
              calibrated on known queries that pass the absence test, so that their true
              taxon is left out at most alpha2 of the time

The calibration queries are split in two: one half calibrates the absence test, the other
the width d. For a known query (own taxon in the reference) every error, called absent,
identified as another taxon, or unresolved among wrong candidates, then happens at most
alpha1 + alpha2 = alpha of the time. Calls on held-out queries (taxon not in the
reference) are not covered by the guarantee; they are measured.

Sets follow the query: one with a 99% hit is compared with taxa near 99%, one with an 80%
hit with taxa near 80%. (A single identity cutoff for all queries had to fall to 59% at
genus rank, putting hundreds of genera into every set.)

Inputs come from best_per_taxon.py. Genus rank: known = held-out-species (LSO) queries,
held-out = held-out-genus (LGO) queries, each against its own split's reference.

  python conformal_tiers.py --rank g \\
      --known-best lso_genus_best.tsv --known-taxonomy out_lso_its2/taxonomy.tsv \\
      --heldout-best lgo_genus_best.tsv --out-prefix results/tiers_genus

Add --apply-best new_best.tsv to label new queries, calibrated on all known queries.
"""
import argparse, bisect, math, os
from collections import defaultdict
import numpy as np

RANK_NAME = {"s": "species", "g": "genus", "f": "family", "o": "order"}

def read_best(path):
    best, order = {}, []
    with open(path) as fh:
        fh.readline()
        for line in fh:
            q, t, s = line.rstrip("\n").split("\t")
            if q not in best:
                best[q] = {}; order.append(q)
            if t != "-":
                best[q][t] = float(s)
    return order, best

def read_truth(path, rank):
    out = {}
    with open(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        i, j = head.index("seq_id"), head.index(rank)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            out[f[i]] = f[j]
    return out

FLOOR = 50.0   # vsearch --id used for best_per_taxon; taxa below it were never seen

class Q:
    """A query's per-taxon identities. Taxa with no hit are only known to sit below FLOOR, so a
    set reaching below FLOOR must include every taxon (it is vacuous: 'could be any taxon')."""
    __slots__ = ("taxa", "ids_desc", "ids_asc", "best", "own")
    def __init__(self, d, truth=None):
        items = sorted(d.items(), key=lambda x: -x[1])
        self.taxa = [t for t, _ in items]
        self.ids_desc = [s for _, s in items]
        self.ids_asc = self.ids_desc[::-1]
        self.best = self.ids_desc[0] if items else 0.0
        self.own = d.get(truth, 0.0) if truth is not None else None
    def vacuous(self, delta):
        return self.best - delta <= FLOOR
    def size(self, delta):
        """Taxa within delta points of the best hit; -1 when the set reaches below FLOOR."""
        if self.vacuous(delta):
            return -1
        return len(self.ids_asc) - bisect.bisect_left(self.ids_asc, self.best - delta)
    def gap(self):
        """Smallest width that puts the true taxon in the set."""
        return self.best - (self.own if self.own > 0 else FLOOR)
    def covered(self, delta):
        return self.gap() <= delta
    def members(self, delta, k=None):
        n = self.size(delta)
        return list(zip(self.taxa[:n], self.ids_desc[:n]))[: (k or n)]

def pvals(cal_best, test_best):
    c = np.sort(-np.asarray(cal_best, float))
    return (1 + len(c) - np.searchsorted(c, -np.asarray(test_best, float), side="left")) / (len(c) + 1)

def upper_q(scores, alpha):
    s = np.sort(np.asarray(scores, float)); n = len(s)
    k = math.ceil((n + 1) * (1 - alpha))
    return math.inf if k > n or n == 0 else s[k - 1]

def calibrate(c1, c2, a1, a2):
    """c1 calibrates the absence test; c2's queries that pass it calibrate the width d."""
    best1 = np.array([q.best for q in c1])
    p2 = pvals(best1, [q.best for q in c2])
    gaps = [q.gap() for q, p in zip(c2, p2) if p > a1]
    return best1, upper_q(gaps, a2)

def call(q, p, a1, delta):
    if p <= a1:
        return "absent", 0
    k = q.size(delta)
    return ("identified", 1) if k == 1 else ("unresolved", k)   # k == -1: could be any taxon

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rank", required=True, choices=list(RANK_NAME))
    ap.add_argument("--known-best", required=True)
    ap.add_argument("--known-taxonomy", required=True, help="taxonomy.tsv of the known queries' split")
    ap.add_argument("--heldout-best", help="best_per_taxon output for queries whose taxon is not in the reference")
    ap.add_argument("--alphas", default="0.05,0.10,0.20", help="total error budget(s)")
    ap.add_argument("--absence-share", type=float, default=0.5, help="share of alpha spent on the absence test")
    ap.add_argument("--splits", type=int, default=200)
    ap.add_argument("--calib-frac", type=float, default=0.5)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--mixes", default="0.01,0.02,0.05,0.10,0.20,0.50")
    ap.add_argument("--transfer-best", help="another known population (own taxon in its reference) to test the calls on")
    ap.add_argument("--transfer-taxonomy", help="taxonomy.tsv for --transfer-best")
    ap.add_argument("--transfer-label", default="transfer", help="name for that population in the summary")
    ap.add_argument("--apply-best", help="best_per_taxon output for new queries to label")
    ap.add_argument("--apply-alpha", type=float, default=0.10)
    ap.add_argument("--floor", type=float, default=50.0, help="the vsearch --id used, in percent")
    ap.add_argument("--out-prefix", required=True)
    a = ap.parse_args()
    global FLOOR
    FLOOR = a.floor
    alphas = [float(x) for x in a.alphas.split(",")]
    mixes = [float(x) for x in a.mixes.split(",")]
    sh = a.absence_share
    rname = RANK_NAME[a.rank]
    os.makedirs(os.path.dirname(a.out_prefix) or ".", exist_ok=True)
    lines = []
    def say(s=""):
        print(s); lines.append(s)

    kid, kbest = read_best(a.known_best)
    truth = read_truth(a.known_taxonomy, a.rank)
    miss = [q for q in kid if q not in truth]
    if miss:
        raise SystemExit(f"{len(miss)} known queries have no {rname} in {a.known_taxonomy}, e.g. {miss[:3]}")
    K = [Q(kbest[q], truth[q]) for q in kid]
    hid, H = [], []
    if a.heldout_best:
        hid, hbest = read_best(a.heldout_best)
        H = [Q(hbest[q]) for q in hid]
    say(f"Three-level calls at {rname} rank: {len(K)} known queries" + (f", {len(H)} held-out queries" if H else ""))
    say(f"  known queries with no hit to their own {rname}: {sum(q.own == 0 for q in K)}; with no hit at all: "
        f"{sum(q.best == 0 for q in K)}" + (f"; held-out queries with no hit at all: {sum(q.best == 0 for q in H)}" if H else ""))
    T = []
    if a.transfer_best:
        tid, tbest = read_best(a.transfer_best)
        ttruth = read_truth(a.transfer_taxonomy, a.rank)
        T = [Q(tbest[q], ttruth[q]) for q in tid if q in ttruth]
        say(f"  {a.transfer_label}: {len(T)} further known queries, to test calls calibrated on the first set")
    say(f"  alpha is split {sh:g} : {1 - sh:g} between the absence test and the set width")
    say()

    rng = np.random.default_rng(a.seed)
    nk = len(K); ncal = int(nk * a.calib_frac)
    acc = defaultdict(list)
    first = None
    for sp in range(a.splits):
        perm = rng.permutation(nk); cal, tst = perm[:ncal], perm[ncal:]
        c1, c2 = cal[: ncal // 2], cal[ncal // 2:]
        calbest = np.array([K[i].best for i in cal])
        pk_full = pvals(calbest, [K[i].best for i in tst])
        ph_full = pvals(calbest, [q.best for q in H]) if H else None
        for al in alphas:
            a1, a2 = sh * al, (1 - sh) * al
            best1, delta = calibrate([K[i] for i in c1], [K[i] for i in c2], a1, a2)
            acc[(al, "width")].append(delta)
            pk = pvals(best1, [K[i].best for i in tst])
            cnt = defaultdict(int); sizes = []
            for i, p in zip(tst, pk):
                q = K[i]; tier, k = call(q, p, a1, delta)
                ok = tier != "absent" and q.covered(delta)
                cnt["error"] += not ok
                if tier == "absent": cnt["absent"] += 1
                elif tier == "identified": cnt["id_ok" if ok else "id_wrong"] += 1
                else:
                    cnt["unres_ok" if ok else "unres_miss"] += 1
                    if k == -1: cnt["vacuous"] += 1
                    else: sizes.append(k)
            for x in ("error", "absent", "id_ok", "id_wrong", "unres_ok", "unres_miss", "vacuous"):
                acc[(al, "known_" + x)].append(cnt[x] / len(tst))
            acc[(al, "known_unres_size_median")].append(float(np.median(sizes)) if sizes else 0.0)
            acc[(al, "test_alone_false_absence")].append(float(np.mean(pk_full <= al)))
            if H:
                ph = pvals(best1, [q.best for q in H])
                cnt = defaultdict(int); sizes = []
                for q, p in zip(H, ph):
                    tier, k = call(q, p, a1, delta); cnt[tier] += 1
                    if tier == "unresolved":
                        if k == -1: cnt["vacuous"] += 1
                        else: sizes.append(k)
                for x in ("absent", "identified", "unresolved", "vacuous"):
                    acc[(al, "held_" + x)].append(cnt[x] / len(H))
                acc[(al, "held_unres_size_median")].append(float(np.median(sizes)) if sizes else 0.0)
                acc[(al, "test_alone_detection")].append(float(np.mean(ph_full <= al)))
            if T:
                cnt = defaultdict(int); sizes = []
                for q, p in zip(T, pvals(best1, [q.best for q in T])):
                    tier, k = call(q, p, a1, delta)
                    ok = tier != "absent" and q.covered(delta)
                    cnt["error"] += not ok
                    if tier == "absent": cnt["absent"] += 1
                    elif tier == "identified": cnt["id_ok" if ok else "id_wrong"] += 1
                    else:
                        cnt["unres_ok" if ok else "unres_miss"] += 1
                        if k == -1: cnt["vacuous"] += 1
                        else: sizes.append(k)
                for x in ("error", "absent", "id_ok", "id_wrong", "unres_ok", "unres_miss", "vacuous"):
                    acc[(al, "transfer_" + x)].append(cnt[x] / len(T))
                acc[(al, "transfer_unres_size_median")].append(float(np.median(sizes)) if sizes else 0.0)
            if sp == 0:
                first = first or {}
                first[al] = (tst, best1, delta, a1)
    def m(al, x):
        v = acc[(al, x)]
        with np.errstate(invalid="ignore"):
            return float(np.mean(v)) if v else float("nan")
    def sd(al, x):
        v = acc[(al, x)]
        with np.errstate(invalid="ignore"):
            return float(np.std(v)) if v else float("nan")

    say(f"Known queries (own {rname} in the reference), means over {a.splits} calibration splits")
    say(f"  {'alpha':>5} {'width d':>9} {'errors':>6} | {'identified':>10} {'(wrong)':>8} {'unresolved':>10} "
        f"{'(misses)':>8} {'(any)':>6} {'size':>5} | {'absent':>7}")
    for al in alphas:
        w = np.array(acc[(al, "width")]); d = float(np.median(w))
        ds = (f"{d:7.1f}" if np.isfinite(d) else "    all") + ("*" if not np.all(np.isfinite(w)) else " ")
        idn = m(al, "known_id_ok") + m(al, "known_id_wrong"); unr = m(al, "known_unres_ok") + m(al, "known_unres_miss")
        say(f"  {al:5.2f} {ds} {m(al, 'known_error'):6.3f} | {idn:10.3f} {m(al, 'known_id_wrong'):8.3f} {unr:10.3f} "
            f"{m(al, 'known_unres_miss'):8.3f} {m(al, 'known_vacuous'):6.3f} {m(al, 'known_unres_size_median'):5.0f} | "
            f"{m(al, 'known_absent'):7.3f}")
    if any(not np.all(np.isfinite(acc[(al, "width")])) for al in alphas):
        say("  * in some calibration splits no finite width met the budget; width shown is the median over splits")
    say("  errors = called absent, identified as another taxon, or unresolved without the true taxon (bounded by alpha);")
    say("  width d = identity points below the best hit a taxon may sit and still be in the set;"
        " (any) = sets reaching below the")
    say(f"  {FLOOR:g}% search floor, so the query could belong to any taxon; size = median candidates of the other"
        " unresolved sets")
    say()
    if T:
        say(f"{a.transfer_label}: calls calibrated on the known queries above, applied to these known queries")
        say(f"  {'alpha':>5} {'errors':>7} | {'identified':>10} {'(wrong)':>8} {'unresolved':>10} {'(misses)':>8} "
            f"{'(any)':>6} {'size':>5} | {'absent':>7}")
        for al in alphas:
            idn = m(al, "transfer_id_ok") + m(al, "transfer_id_wrong"); unr = m(al, "transfer_unres_ok") + m(al, "transfer_unres_miss")
            say(f"  {al:5.2f} {m(al, 'transfer_error'):7.3f} | {idn:10.3f} {m(al, 'transfer_id_wrong'):8.3f} {unr:10.3f} "
                f"{m(al, 'transfer_unres_miss'):8.3f} {m(al, 'transfer_vacuous'):6.3f} "
                f"{m(al, 'transfer_unres_size_median'):5.0f} | {m(al, 'transfer_absent'):7.3f}")
        say()
    if H:
        say(f"Held-out queries ({rname} not in the reference): an identified call is always wrong at this rank")
        say(f"  {'alpha':>5} {'absent':>7} {'identified':>11} {'unresolved':>11} {'(any)':>6} {'size':>5} | "
            f"{'absence test alone at alpha: detection':>39} {'false-absence':>14}")
        for al in alphas:
            say(f"  {al:5.2f} {m(al, 'held_absent'):7.3f} {m(al, 'held_identified'):11.3f} {m(al, 'held_unresolved'):11.3f} "
                f"{m(al, 'held_vacuous'):6.3f} "
                f"{m(al, 'held_unres_size_median'):5.0f} | {m(al, 'test_alone_detection'):39.3f} "
                f"{m(al, 'test_alone_false_absence'):14.3f}")
        say("  absence test alone = the single absence call at the full alpha, for comparison")
        say()
        say("Share of calls that are wrong, by share of held-out queries in a sample (pi)")
        say(f"  {'alpha':>5} {'call':>11} " + "".join(f"{'pi=' + format(p, 'g'):>9}" for p in mixes))
        rows_mix = []
        for al in alphas:
            AK, AN = m(al, "known_absent"), m(al, "held_absent")
            IK = m(al, "known_id_ok") + m(al, "known_id_wrong"); IWK = m(al, "known_id_wrong"); IN = m(al, "held_identified")
            ab = [(1 - p) * AK / ((1 - p) * AK + p * AN) if (AK + AN) > 0 else float("nan") for p in mixes]
            idw = [((1 - p) * IWK + p * IN) / ((1 - p) * IK + p * IN) if (IK + IN) > 0 else float("nan") for p in mixes]
            say(f"  {al:5.2f} {'absent':>11} " + "".join(f"{x:9.3f}" for x in ab))
            say(f"  {'':5} {'identified':>11} " + "".join(f"{x:9.3f}" for x in idw))
            rows_mix += [(al, "absent", p, x) for p, x in zip(mixes, ab)] + [(al, "identified", p, x) for p, x in zip(mixes, idw)]
        say()

    P = a.out_prefix
    with open(P + "_summary.txt", "w") as fh:
        fh.write("\n".join(lines) + "\n")
    keys = sorted({k for (_, k) in acc})
    with open(P + "_tiers.tsv", "w") as fh:
        fh.write("alpha\tmeasure\tmean\tsd\n")
        for al in alphas:
            for k in keys:
                if acc[(al, k)]:
                    fh.write(f"{al}\t{k}\t{m(al, k):.5f}\t{sd(al, k):.5f}\n")
    if H:
        with open(P + "_mixtures.tsv", "w") as fh:
            fh.write("alpha\tcall\tpi\tshare_wrong\n")
            fh.writelines(f"{al}\t{c}\t{p}\t{x:.5f}\n" for al, c, p, x in rows_mix)
    al0 = 0.10 if 0.10 in first else alphas[0]
    tst, best1, delta, a1 = first[al0]
    with open(P + "_perquery.tsv", "w") as fh:
        fh.write(f"query_id\tpopulation\ttrue_{rname}\ttier_a{al0:g}\tset_size\ttruth_in_set\tbest_identity\tcandidates\n")
        pk = pvals(best1, [K[i].best for i in tst])
        for i, p in zip(tst, pk):
            q = K[i]; tier, k = call(q, p, a1, delta)
            cand = ("any" if k == -1 else ",".join(f"{t}:{s:.1f}" for t, s in q.members(delta, 8))) if tier != "absent" else ""
            fh.write(f"{kid[i]}\tknown\t{truth[kid[i]]}\t{tier}\t{k}\t{int(tier != 'absent' and q.covered(delta))}\t{q.best:.1f}\t{cand}\n")
        if H:
            for qid, q, p in zip(hid, H, pvals(best1, [q.best for q in H])):
                tier, k = call(q, p, a1, delta)
                cand = ("any" if k == -1 else ",".join(f"{t}:{s:.1f}" for t, s in q.members(delta, 8))) if tier != "absent" else ""
                fh.write(f"{qid}\theld-out\t\t{tier}\t{k}\t0\t{q.best:.1f}\t{cand}\n")

    if a.apply_best:
        aid, abest = read_best(a.apply_best)
        al = a.apply_alpha; a1, a2 = sh * al, (1 - sh) * al
        perm = np.random.default_rng(a.seed).permutation(nk)
        best1, delta = calibrate([K[i] for i in perm[: nk // 2]], [K[i] for i in perm[nk // 2:]], a1, a2)
        A = [Q(abest[x]) for x in aid]
        pa = pvals(best1, [q.best for q in A])
        pall = pvals([q.best for q in K], [q.best for q in A])   # single absence call, all known queries
        counts = defaultdict(int)
        with open(P + f"_apply_a{al:g}.tsv", "w") as fh:
            fh.write("query_id\ttier\tset_size\tbest_identity\tabsence_p_tier\tabsence_p_single\tcandidates\n")
            for x, q, p, p2 in zip(aid, A, pa, pall):
                tier, k = call(q, p, a1, delta); counts[tier] += 1
                cand = ("any" if k == -1 else ",".join(f"{t}:{s:.1f}" for t, s in q.members(delta, 10))) if tier != "absent" else ""
                fh.write(f"{x}\t{tier}\t{k}\t{q.best:.2f}\t{p:.5g}\t{p2:.5g}\t{cand}\n")
        n = len(aid)
        say(f"Applied to {n} new queries at alpha={al:g} (width d = {delta:.1f}, calibrated on all {nk} known queries): "
            + ", ".join(f"{t} {counts[t]} ({counts[t] / n:.1%})" for t in ("identified", "unresolved", "absent")))
        say("  single absence call, calibrated on all known queries: " + ", ".join(
            f"alpha {x:g}: {int(np.sum(pall <= x))} ({np.mean(pall <= x):.1%})" for x in (0.05, 0.10, 0.20)))
        with open(P + "_summary.txt", "a") as fh:
            fh.write(lines[-2] + "\n" + lines[-1] + "\n")

    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        cols = {"identified, right": "#009E73", "identified, wrong": "#E69F00",
                "unresolved, includes truth": "#BBBBBB", "unresolved, misses truth": "#666666",
                "unresolved, could be any taxon": "#E5E5E5", "absent": "#0072B2"}
        fig, ax = plt.subplots(1, 2 if H else 1, figsize=(6.85 if H else 3.6, 2.8), squeeze=False)
        x = np.arange(len(alphas))
        val = lambda k: np.array([m(al, k) for al in alphas])
        parts = [("identified, right", val("known_id_ok")), ("identified, wrong", val("known_id_wrong")),
                 ("unresolved, includes truth", val("known_unres_ok") - val("known_vacuous")),
                 ("unresolved, misses truth", val("known_unres_miss")),
                 ("unresolved, could be any taxon", val("known_vacuous")), ("absent", val("known_absent"))]
        bottom = np.zeros(len(alphas))
        for lab, v in parts:
            ax[0][0].bar(x, v, bottom=bottom, color=cols[lab], label=lab, width=0.7, edgecolor="0.6", linewidth=0.3); bottom += v
        for i, al in enumerate(alphas):
            ax[0][0].text(i, 1.01, f"{m(al, 'known_unres_size_median'):.0f}", ha="center", va="bottom", fontsize=6.5, color="0.35")
        ax[0][0].set_title(f"A  {rname} in the reference", loc="left", fontsize=9, fontweight="bold", pad=12)
        if H:
            bottom = np.zeros(len(alphas))
            for lab, v in (("identified, wrong", val("held_identified")),
                           ("unresolved, misses truth", val("held_unresolved") - val("held_vacuous")),
                           ("unresolved, could be any taxon", val("held_vacuous")), ("absent", val("held_absent"))):
                ax[0][1].bar(x, v, bottom=bottom, color=cols[lab], width=0.7, edgecolor="0.6", linewidth=0.3); bottom += v
            for i, al in enumerate(alphas):
                ax[0][1].text(i, 1.01, f"{m(al, 'held_unres_size_median'):.0f}", ha="center", va="bottom", fontsize=6.5, color="0.35")
            ax[0][1].set_title(f"B  {rname} absent from the reference", loc="left", fontsize=9, fontweight="bold", pad=12)
        for axx in ax[0]:
            axx.set_xticks(x); axx.set_xticklabels([f"{al:g}" for al in alphas]); axx.set_xlabel("total error budget \u03b1")
            axx.set_ylim(0, 1); axx.set_ylabel("share of queries"); axx.tick_params(labelsize=7)
        fig.legend(*ax[0][0].get_legend_handles_labels(), loc="lower center", ncol=3, fontsize=7, frameon=False)
        fig.tight_layout(rect=(0, 0.14, 1, 1))
        fig.savefig(P + "_tiers.pdf"); fig.savefig(P + "_tiers.png", dpi=300)
    except ImportError:
        pass
    print(f"[tiers] wrote {P}_summary.txt, _tiers.tsv, _perquery.tsv" + (", _mixtures.tsv" if H else "")
          + (f", _apply_a{a.apply_alpha:g}.tsv" if a.apply_best else "") + ", _tiers.pdf/.png")

if __name__ == "__main__":
    main()
