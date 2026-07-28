#!/usr/bin/env python3
"""Path B step 1: mine GlobalFungi sample metadata to shortlist independent ITS2
studies for the cross-study reproducibility test."""
import argparse, gzip, csv
from collections import defaultdict

def openmaybe(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--region", default="ITS2")
    ap.add_argument("--platform", default="")
    ap.add_argument("--habitat", default="")
    ap.add_argument("--min-samples", type=int, default=20)
    ap.add_argument("--out", default="study_shortlist.tsv")
    a = ap.parse_args()
    rl, pl, hl = a.region.lower(), a.platform.lower(), a.habitat.lower()
    S = defaultdict(lambda: {"n":0,"primers":defaultdict(int),"platform":set(),"region":set(),
                             "habitat":defaultdict(int),"continent":set(),"year":"","doi":"","title":""})
    n = 0
    with openmaybe(a.metadata) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            n += 1
            region = (row.get("barcoding_region") or "").strip()
            plat   = (row.get("sequencing_platform") or "").strip()
            styp   = (row.get("sample_type") or "").strip()
            env    = (row.get("environment_type") or "").strip()
            if rl and rl not in region.lower(): continue
            if pl and pl not in plat.lower(): continue
            if hl and hl not in (styp + " " + env).lower(): continue
            s = S[row.get("paper_ID") or "NA"]
            s["n"] += 1
            s["primers"][(row.get("PCR_primers") or "NA").strip()] += 1
            s["platform"].add(plat); s["region"].add(region)
            s["habitat"][styp or env] += 1
            s["continent"].add((row.get("continent") or "").strip())
            s["year"] = (row.get("paper_year") or "").strip()
            s["doi"]  = (row.get("paper_doi") or "").strip()
            s["title"] = (row.get("paper_title") or "").strip()
    kept = {p:s for p,s in S.items() if s["n"] >= a.min_samples}
    def topprimer(s): return max(s["primers"].items(), key=lambda kv: kv[1])[0]
    cohorts = defaultdict(list)
    for pid, s in kept.items(): cohorts[topprimer(s)].append(pid)
    print(f"[pick] rows read: {n}")
    print(f"[pick] studies matching region='{a.region}' platform~'{a.platform}' "
          f"habitat~'{a.habitat}' with >= {a.min_samples} samples: {len(kept)}")
    print(f"\n[pick] independent-study cohorts sharing a primer pair:")
    for primer, pids in sorted(cohorts.items(), key=lambda kv: -len(kv[1])):
        nsamp = sum(kept[p]["n"] for p in pids)
        print(f"   {len(pids):3d} studies  {nsamp:6d} samples   primer={primer}")
    with open(a.out, "w") as o:
        o.write("paper_ID\tyear\tn_samples\ttop_primer\tplatform\tcontinent\thabitat\tdoi\ttitle\n")
        for pid, s in sorted(kept.items(), key=lambda kv: (topprimer(kv[1]), -kv[1]["n"])):
            hab = max(s["habitat"].items(), key=lambda kv: kv[1])[0] if s["habitat"] else ""
            o.write("\t".join([pid, s["year"], str(s["n"]), topprimer(s),
                               "/".join(sorted(x for x in s["platform"] if x)),
                               "/".join(sorted(x for x in s["continent"] if x)),
                               hab, s["doi"], s["title"][:80]]) + "\n")
    print(f"\n[pick] wrote per-study shortlist -> {a.out} (shared-primer studies cluster together)")

if __name__ == "__main__": main()
