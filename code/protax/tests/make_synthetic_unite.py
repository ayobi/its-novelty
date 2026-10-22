#!/usr/bin/env python3
"""Synthetic UNITE-style ITS2 data for testing the PROTAX harness end to end.

Writes, into --out:
  unite_synth.fasta     UNITE-style headers, sequences evolved down a 6-rank tree
  official_taxonomy     a PROTAX-format taxonomy (id pid level name prior) that mimics
                        PROTAX-fungi's: its per-species prior scale, Index Fungorum-style
                        numeric suffixes on genus and
                        species names, extra named-but-unsequenced genera and species,
                        a nonFungi branch, and two sequenced genera deliberately missing
                        (so the builder has to create nodes for them)
  soil_asvs.fasta       a few query sequences, some from unsampled genera

Plumbing test only: the numbers it produces mean nothing biologically.
"""
import argparse, os, random

BASES = "ACGT"

def mutate(seq, rate, rng):
    out = []
    for b in seq:
        r = rng.random()
        if r < rate * 0.85:
            out.append(rng.choice([x for x in BASES if x != b]))
        elif r < rate * 0.92:
            continue                                   # deletion
        elif r < rate:
            out.append(b); out.append(rng.choice(BASES))  # insertion
        else:
            out.append(b)
    return "".join(out)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    rng = random.Random(a.seed)
    os.makedirs(a.out, exist_ok=True)

    # per-rank divergence: phylum, class, order, family, genus, species, individual
    rates = [0.20, 0.10, 0.07, 0.06, 0.05, 0.025, 0.006]
    root = "".join(rng.choice(BASES) for _ in range(260))
    recs, acc = [], 1
    tax_lines = []          # (level, path tuple, sequenced species count, unsequenced species count)
    gid = 100000
    for p in range(2):
        ps = mutate(root, rates[0], rng); P = f"Phy{p}mycota"
        for c in range(2):
            cs = mutate(ps, rates[1], rng); C = f"{P[:4]}Cls{c}mycetes"
            for o in range(2):
                os_ = mutate(cs, rates[2], rng); O = f"{C[:6]}Ord{o}ales"
                for f in range(3):
                    fs = mutate(os_, rates[3], rng); F = f"{O[:7]}Fam{f}aceae"
                    for g in range(4):
                        gs = mutate(fs, rates[4], rng); G = f"Gen{p}{c}{o}{f}{g}"
                        n_sp = 3
                        for s in range(n_sp):
                            ss = mutate(gs, rates[5], rng); S = f"{G}_epith{s}"
                            for _ in range(3):
                                seq = mutate(ss, rates[6], rng)
                                header = (f"{S}|UDB{acc:07d}|SH{acc+900000}.10FU|reps|"
                                          f"k__Fungi;p__{P};c__{C};o__{O};f__{F};g__{G};s__{S}")
                                recs.append((header, seq)); acc += 1
                        # one genus-only (placeholder species) sequence per genus
                        seq = mutate(mutate(gs, rates[5], rng), rates[6], rng)
                        recs.append((f"{G}_sp|UDB{acc:07d}|SH{acc+900000}.10FU|reps|"
                                     f"k__Fungi;p__{P};c__{C};o__{O};f__{F};g__{G};s__{G}_sp", seq)); acc += 1
    with open(os.path.join(a.out, "unite_synth.fasta"), "w") as fh:
        for h, s in recs:
            fh.write(f">{h}\n{s}\n")

    # ---- PROTAX-format "official" taxonomy ----
    lineages = {}
    for h, _ in recs:
        lin = h.split("|")[-1].split(";")
        r = {t.split("__")[0]: t.split("__")[1] for t in lin}
        lineages.setdefault(r["g"], r)
    missing = {"Gen00000", "Gen10102"}          # sequenced in "UNITE", absent from the taxonomy
    nodes = [(0, 0, 0, "root", 1.0)]
    nid = 1
    nodes.append((nid, 0, 1, "Fungi", 0.99)); fungi = nid; nid += 1
    nodes.append((nid, 0, 1, "nonFungi", 0.01)); nonf = nid; nid += 1
    nodes.append((nid, nonf, 2, "nonFungi,Plantae", 0.01)); nid += 1
    path_id = {(): fungi}
    sp_rows = []
    def ensure(path, level):
        nonlocal nid
        if path in path_id: return path_id[path]
        parent = ensure(path[:-1], level - 1)
        path_id[path] = nid
        nodes.append([nid, parent, level, "Fungi," + ",".join(path), 0.0]); nid += 1
        return path_id[path]
    for g, r in sorted(lineages.items()):
        if g in missing: continue
        base = (r["p"], r["c"], r["o"], r["f"])
        gpath = base + (f"{g}_{gid}",); gid += 1
        ensure(gpath, 6)
        for s in range(3):
            sp_rows.append((gpath + (f"{g}_epith{s}_{gid}",))); gid += 1
        for s in range(rng.randint(0, 3)):           # named, never sequenced
            sp_rows.append((gpath + (f"{g}_ignotus{s}_{gid}",))); gid += 1
    fams = sorted({(r["p"], r["c"], r["o"], r["f"]) for r in lineages.values()})
    for fi, fam in enumerate(fams):                   # named genera with no sequences at all
        for k in range(2):
            gpath = fam + (f"Nudogen{fi}x{k}_{gid}",); gid += 1
            ensure(gpath, 6)
            for s in range(rng.randint(1, 4)):
                sp_rows.append(gpath + (f"Nudogen{fi}x{k}_ignotus{s}_{gid}",)); gid += 1
    for sp in sp_rows:
        ensure(sp, 7)
    sp_prior = 7.5294e-6     # PROTAX-fungi's per-species prior: this tree is a corner of a big one
    # absolute priors proportional to species counts, as in PROTAX-fungi's taxonomy
    count = {}
    for sp in sp_rows:
        for k in range(1, len(sp) + 1):
            count[sp[:k]] = count.get(sp[:k], 0) + 1
    out = []
    for n in nodes:
        n = list(n)
        if n[2] >= 2 and n[3].startswith("Fungi,"):
            path = tuple(n[3].split(",")[1:])
            n[4] = sp_prior * count.get(path, 1)
        out.append(n)
    with open(os.path.join(a.out, "official_taxonomy"), "w") as fh:
        for n in out:
            fh.write(f"{n[0]}\t{n[1]}\t{n[2]}\t{n[3]}\t{n[4]:.10f}\n")

    # ---- a handful of "soil" ASVs: some near sampled genera, some far ----
    with open(os.path.join(a.out, "soil_asvs.fasta"), "w") as fh:
        for i in range(30):
            h, s = recs[rng.randrange(len(recs))]
            s = mutate(s, 0.02 if i % 3 else 0.25, rng)
            fh.write(f">ASV{i+1}\n{s}\n")
    print(f"[synth] {len(recs)} sequences, {len(out)} taxonomy nodes, sp_prior={sp_prior:.3g}")

if __name__ == "__main__":
    main()
