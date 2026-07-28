#!/usr/bin/env python3
"""Calibrated novelty flag via inductive conformal anomaly detection.
Calibrate on held-out KNOWN queries; p<=alpha flags novel while controlling the
false-novelty rate on knowns at alpha. Score-agnostic (identity or cosine)."""
import argparse, random
import numpy as np

def best_scores(path):
    best = {}
    with open(path) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 3: continue
            try: s = float(p[2])
            except ValueError: continue
            q = p[0]
            if q not in best or s > best[q]: best[q] = s
    return best

def pvalues(calib_novelty, test_novelty):
    c = np.sort(np.asarray(calib_novelty, float)); n = len(c)
    ge = n - np.searchsorted(c, np.asarray(test_novelty, float), side="left")
    return (1 + ge) / (n + 1)

def auroc(labels, scores):
    labels = np.asarray(labels, float); scores = np.asarray(scores, float)
    order = np.argsort(scores, kind="mergesort"); ss = scores[order]
    r = np.empty(len(scores)); i = 0
    while i < len(scores):
        j = i
        while j + 1 < len(scores) and ss[j + 1] == ss[i]: j += 1
        r[order[i:j + 1]] = (i + j) / 2.0 + 1; i = j + 1
    npos = labels.sum(); nneg = len(labels) - npos
    return float("nan") if npos == 0 or nneg == 0 else \
        (r[labels == 1].sum() - npos * (npos + 1) / 2) / (npos * nneg)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--known", required=True); ap.add_argument("--novel", required=True)
    ap.add_argument("--label", default="signal"); ap.add_argument("--calib-frac", type=float, default=0.5)
    ap.add_argument("--alphas", default="0.01,0.05,0.10,0.20"); ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    known = best_scores(a.known); novel = best_scores(a.novel)
    kq = list(known); random.Random(a.seed).shuffle(kq)
    ncal = int(len(kq) * a.calib_frac); calib_ids, testk_ids = kq[:ncal], kq[ncal:]; nq = list(novel)
    calib_s = np.array([-known[q] for q in calib_ids])
    testk_s = np.array([-known[q] for q in testk_ids]); novel_s = np.array([-novel[q] for q in nq])
    p_testk = pvalues(calib_s, testk_s); p_novel = pvalues(calib_s, novel_s)
    au = auroc([0]*len(testk_s)+[1]*len(novel_s), np.concatenate([testk_s, novel_s]))
    calib_sim = np.array([known[q] for q in calib_ids])
    print(f"=== conformal novelty [{a.label}] ===")
    print(f"  calibration knowns : {len(calib_ids)}")
    print(f"  test knowns        : {len(testk_ids)}    test novels: {len(nq)}")
    print(f"  AUROC (test split) : {au:.4f}")
    print(f"  {'alpha':>6} {'sim_cutoff':>11} {'known_falseflag':>16} {'novel_detected':>15}")
    for alpha in [float(x) for x in a.alphas.split(",")]:
        cutoff = np.quantile(calib_sim, alpha)
        ff = float(np.mean(p_testk <= alpha)); det = float(np.mean(p_novel <= alpha))
        print(f"  {alpha:6.2f} {cutoff:11.3f} {ff:16.4f} {det:15.4f}")
    if a.out:
        with open(a.out, "w") as o:
            o.write("query_id\tpopulation\tbest_sim\tp_value\n")
            for q, p in zip(testk_ids, p_testk): o.write(f"{q}\tknown\t{known[q]:.4f}\t{p:.4f}\n")
            for q, p in zip(nq, p_novel): o.write(f"{q}\tnovel\t{novel[q]:.4f}\t{p:.4f}\n")
        print(f"  wrote per-query p-values to {a.out}")

if __name__ == "__main__": main()
