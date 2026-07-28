#!/usr/bin/env python3
"""
Phase-0 step 1a: prepare a held-out benchmark split from UNITE ITS2.

No third-party dependencies (stdlib only). Produces, in the output dir:
  refs.fasta      -- reference sequences (build the TaxoTagger/BLAST/VSEARCH DB from THIS only)
  queries.fasta   -- held-out query sequences (the "novel" fungi)
  taxonomy.tsv    -- seq_id -> p,c,o,f,g,s  (for every input seq; used by later metrics)
  split.tsv       -- seq_id -> role(ref|query), true_genus, true_family (the answer key)

Split modes:
  lgo  leave-one-genus-out : hold out one genus per family that has >=2 genera
                             (query's true FAMILY stays in refs)  -> DEEP novelty
  lso  leave-one-species-out: hold out one species per genus that has >=2 species
                             (query's true GENUS stays in refs)   -> SHALLOW novelty

Run the self-test first:   python prep_split.py --selftest
Then on real data:         python prep_split.py -i unite_its2.fasta -o out_lgo --mode lgo --seed 0
"""

import argparse
import random
import sys
from collections import defaultdict
import re
ACC_RE = re.compile(r'^(UDB\d+|[A-Z]{1,3}\d{5,8}(\.\d+)?)$')
SH_RE  = re.compile(r'^SH\d+\.\d+\w*$')

RANKS = ["p", "c", "o", "f", "g", "s"]
# Placeholder / unreliable annotations to drop from the backbone.
BAD = {"", "unidentified", "incertae_sedis", "environmental", "unclassified", "sp", "sp."}


# ---------- parsing ----------

def parse_fasta(path):
    """Yield (header_without_gt, sequence). Minimal, memory-light FASTA reader."""
    header, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq)
                header, seq = line[1:], []
            else:
                seq.append(line)
    if header is not None:
        yield header, "".join(seq)


def parse_unite_header(header):
    parts = [p.strip() for p in header.split("|")]
    ranks = {r: "" for r in RANKS}
    lineage_field = next((p for p in parts if "__" in p and ";" in p), "")
    for tok in lineage_field.split(";"):
        tok = tok.strip()
        if "__" in tok:
            pref, _, val = tok.partition("__")
            if pref in ranks:
                ranks[pref] = val.strip()
    seq_id = next((p for p in parts if ACC_RE.match(p)), None)
    if seq_id is None:
        seq_id = next((p for p in parts if SH_RE.match(p)), None)
    if seq_id is None:
        seq_id = next((p for p in parts if p and p != lineage_field), header.split()[0])
    return seq_id, ranks


def norm(x):
    return (x or "").strip().lower().replace(" ", "_")


def is_placeholder_species(sp):
    s = norm(sp)
    return s in BAD or s.endswith("_sp.") or s.endswith("_sp") or s == ""


# ---------- backbone cleaning ----------

def clean_backbone(records, need_species):
    """
    records: dict seq_id -> ranks(dict). Returns a filtered dict.
    Requires real genus AND family; if need_species, also a real binomial.
    """
    out = {}
    for sid, r in records.items():
        if "incertae" in norm(r["g"]) or "incertae" in norm(r["f"]) or norm(r["g"]) in BAD or norm(r["f"]) in BAD:
            continue
        if need_species and is_placeholder_species(r["s"]):
            continue
        out[sid] = r
    return out


# ---------- splits ----------

def leave_one_genus_out(records, seed=0):
    """One held-out genus per family with >=2 genera. Held-out genus' FAMILY stays in refs."""
    rng = random.Random(seed)
    fam_to_genera = defaultdict(set)
    for r in records.values():
        fam_to_genera[r["f"]].add(r["g"])
    holdout_genera = set()
    for fam, genera in fam_to_genera.items():
        if len(genera) >= 2:
            holdout_genera.add(rng.choice(sorted(genera)))
    q = [sid for sid, r in records.items() if r["g"] in holdout_genera]
    ref = [sid for sid, r in records.items() if r["g"] not in holdout_genera]
    return set(q), set(ref)


def leave_one_species_out(records, seed=0):
    """One held-out species per genus with >=2 species. Held-out species' GENUS stays in refs."""
    rng = random.Random(seed)
    genus_to_species = defaultdict(set)
    for r in records.values():
        genus_to_species[r["g"]].add(r["s"])
    holdout_species = set()  # keyed by (genus, species) to avoid cross-genus name collisions
    for genus, species in genus_to_species.items():
        if len(species) >= 2:
            holdout_species.add((genus, rng.choice(sorted(species))))
    q = [sid for sid, r in records.items() if (r["g"], r["s"]) in holdout_species]
    ref = [sid for sid, r in records.items() if (r["g"], r["s"]) not in holdout_species]
    return set(q), set(ref)


# ---------- correctness checks ----------

def assert_split_valid(records, q_ids, ref_ids, mode):
    """The whole point of the sibling conditions: every query's true rank must remain in refs."""
    level = "f" if mode == "lgo" else "g"
    ref_levels = {records[s][level] for s in ref_ids}
    missing = sorted({records[s][level] for s in q_ids} - ref_levels)
    assert not missing, (
        f"INVALID SPLIT ({mode}): {len(missing)} query {('family' if level=='f' else 'genus')}(s) "
        f"absent from refs, e.g. {missing[:5]}. Sibling condition failed."
    )
    assert not (q_ids & ref_ids), "A sequence is in both refs and queries."


# ---------- IO ----------

def write_outputs(outdir, records, seqs, q_ids, ref_ids):
    import os
    os.makedirs(outdir, exist_ok=True)

    def write_fasta(path, ids):
        with open(path, "w") as fh:
            for sid in ids:
                fh.write(f">{sid}\n{seqs[sid]}\n")

    write_fasta(os.path.join(outdir, "refs.fasta"), sorted(ref_ids))
    write_fasta(os.path.join(outdir, "queries.fasta"), sorted(q_ids))
    with open(os.path.join(outdir, "taxonomy.tsv"), "w") as fh:
        fh.write("seq_id\t" + "\t".join(RANKS) + "\n")
        for sid, r in records.items():
            fh.write(sid + "\t" + "\t".join(r[x] for x in RANKS) + "\n")
    with open(os.path.join(outdir, "split.tsv"), "w") as fh:
        fh.write("seq_id\trole\ttrue_genus\ttrue_family\n")
        for sid in sorted(records):
            role = "query" if sid in q_ids else "ref"
            fh.write(f"{sid}\t{role}\t{records[sid]['g']}\t{records[sid]['f']}\n")


def summarize(records, q_ids, ref_ids, mode):
    ref_g = {records[s]["g"] for s in ref_ids}
    ref_f = {records[s]["f"] for s in ref_ids}
    print(f"  mode                : {mode}")
    print(f"  backbone sequences  : {len(records)}")
    print(f"  reference sequences : {len(ref_ids)}  ({len(ref_g)} genera, {len(ref_f)} families)")
    print(f"  query sequences     : {len(q_ids)}  "
          f"({'held-out genera' if mode=='lgo' else 'held-out species'})")
    novelty = "family present, genus ABSENT (deep)" if mode == "lgo" else "genus present (shallow)"
    print(f"  novelty regime      : queries have {novelty}")


# ---------- main ----------

def run_real(args):
    seqs, recs = {}, {}
    for header, seq in parse_fasta(args.input):
        sid, ranks = parse_unite_header(header)
        seqs[sid] = seq
        recs[sid] = ranks
    print(f"[read] {len(recs)} sequences from {args.input}")
    need_species = (args.mode == "lso")
    recs = clean_backbone(recs, need_species=need_species)
    print(f"[clean] {len(recs)} remain after backbone filtering "
          f"(need_species={need_species})")
    split = leave_one_genus_out if args.mode == "lgo" else leave_one_species_out
    q_ids, ref_ids = split(recs, seed=args.seed)
    assert_split_valid(recs, q_ids, ref_ids, args.mode)
    print("[check] split valid: every query's true rank is present in refs  ✅")
    write_outputs(args.out, recs, seqs, q_ids, ref_ids)
    print(f"[write] refs.fasta / queries.fasta / taxonomy.tsv / split.tsv -> {args.out}/")
    summarize(recs, q_ids, ref_ids, args.mode)


def run_selftest():
    # Synthetic taxonomy exercising every branch of the sibling conditions.
    # Family F1: G1{sp1,sp2}, G2{sp3}          -> F1 eligible for LGO (>=2 genera)
    # Family F2: G3{sp4,sp5,sp6}                -> F2 NOT eligible for LGO (1 genus)
    # Family F3: G4{sp7}                        -> F3 NOT eligible for LGO (1 genus)
    # LSO-eligible genera (>=2 species): G1, G3.  Ineligible: G2, G4.
    # Plus placeholder rows that MUST be dropped.
    def rec(p, c, o, f, g, s):
        return {"p": p, "c": c, "o": o, "f": f, "g": g, "s": s}

    recs = {}
    spec = {
        "F1": {"G1": ["G1_sp1", "G1_sp2"], "G2": ["G2_sp3"]},
        "F2": {"G3": ["G3_sp4", "G3_sp5", "G3_sp6"]},
        "F3": {"G4": ["G4_sp7"]},
    }
    sid = 0
    for f, genera in spec.items():
        for g, species in genera.items():
            for s in species:
                for _ in range(3):  # 3 sequences per species
                    recs[f"seq{sid}"] = rec("P", "C", "O", f, g, s)
                    sid += 1
    # placeholders (should be filtered out entirely)
    recs["bad1"] = rec("P", "C", "O", "F1", "unidentified", "unidentified")
    recs["bad2"] = rec("P", "C", "O", "F2", "G3", "G3_sp.")
    n_real_genus = 4

    ok = True

    # --- LGO ---
    lgo_recs = clean_backbone(recs, need_species=False)
    assert "bad1" not in lgo_recs, "placeholder genus not dropped"
    q, ref = leave_one_genus_out(lgo_recs, seed=1)
    assert_split_valid(lgo_recs, q, ref, "lgo")
    held_genera = {lgo_recs[s]["g"] for s in q}
    # only a genus from F1 may be held out (only F1 has >=2 genera)
    if not held_genera or not held_genera.issubset({"G1", "G2"}):
        print(f"  [FAIL] LGO held out {held_genera}, expected exactly one of {{G1,G2}}"); ok = False
    else:
        print(f"  [ok]   LGO held out {held_genera} (only F1 is eligible)  ✅")
    ref_fams = {lgo_recs[s]["f"] for s in ref}
    if {lgo_recs[s]["f"] for s in q}.issubset(ref_fams):
        print("  [ok]   LGO: every query family still present in refs  ✅")
    else:
        print("  [FAIL] LGO emptied a family"); ok = False

    # --- LSO ---
    lso_recs = clean_backbone(recs, need_species=True)
    if "bad2" in lso_recs:
        print("  [FAIL] placeholder species 'G3_sp.' not dropped"); ok = False
    else:
        print("  [ok]   LSO: placeholder species dropped  ✅")
    q, ref = leave_one_species_out(lso_recs, seed=1)
    assert_split_valid(lso_recs, q, ref, "lso")
    held_genera_of_species = {lso_recs[s]["g"] for s in q}
    # held-out species may only come from genera with >=2 species: G1, G3
    if held_genera_of_species and held_genera_of_species.issubset({"G1", "G3"}):
        print(f"  [ok]   LSO held species from {held_genera_of_species} (only G1,G3 eligible)  ✅")
    else:
        print(f"  [FAIL] LSO held species from {held_genera_of_species}, expected ⊆ {{G1,G3}}"); ok = False
    if {lso_recs[s]["g"] for s in q}.issubset({lso_recs[s]["g"] for s in ref}):
        print("  [ok]   LSO: every query genus still present in refs  ✅")
    else:
        print("  [FAIL] LSO emptied a genus"); ok = False

    print("\nSELFTEST:", "PASS ✅" if ok else "FAIL ❌")
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser(description="Phase-0 step 1a: UNITE ITS2 held-out split.")
    ap.add_argument("-i", "--input", help="UNITE ITS2 FASTA (headers with x__ lineage)")
    ap.add_argument("-o", "--out", help="output directory")
    ap.add_argument("--mode", choices=["lgo", "lso"], default="lgo")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--selftest", action="store_true", help="run synthetic correctness test and exit")
    args = ap.parse_args()

    if args.selftest:
        sys.exit(run_selftest())
    if not (args.input and args.out):
        ap.error("provide -i/--input and -o/--out (or --selftest)")
    run_real(args)


if __name__ == "__main__":
    main()
