#!/usr/bin/env python3
"""Build a PROTAX-fungi model directory for one reference set.

PROTAX-fungi (Abarenkov et al. 2018, github.com/psomervuo/protaxfungi) ships its Perl
classifier and per-rank parameters (params.its2.level2..7) but no reference data; the
user supplies reference sequences and a taxonomy. This script supplies them so that
PROTAX sees exactly the reference the identity scores were computed against:

  reference sequences  the split's refs.fasta (or every record of a UNITE FASTA)
  taxonomy and priors  PROTAX-fungi's own taxonomy (--protax-taxonomy), so a held-out
                       genus stays in the tree as a named genus without reference
                       sequences, next to every other named-but-unsequenced genus, on
                       the prior scale the parameters were fitted under

Split names are matched to taxonomy nodes by name (Index Fungorum suffixes ignored,
genus first, ties broken by agreement of family, order, class and phylum); names the
taxonomy lacks get new nodes. Writes what run_protax.sh needs: taxonomy, taxonomy.ascii7,
tax1..7, seqid2tax, ref.tax2..7, rseqs2..7, <region>.fa, sintax<region>train.fa,
params.<region>.level2..7, startnode, build_report.txt.

Without --protax-taxonomy the tree is built from the input's own names. That runs, but
the only named genera without sequences are then the held-out ones, which pushes PROTAX
towards "known". Use it for testing, not for the comparison in the paper.

  python build_protax_model.py --split-dir out_lgo_its2 \
      --protax-taxonomy db_files/tax7 --params protaxfungi/model --out model_lgo_its2
"""
import argparse, os, re, shutil, sys
from collections import Counter, defaultdict

RANKS = ["p", "c", "o", "f", "g", "s"]            # taxonomy.tsv columns = PROTAX levels 2..7
SINTAX = {1: "d", 2: "p", 3: "c", 4: "o", 5: "f", 6: "g", 7: "s"}
BAD = {"", "unidentified", "incertae_sedis", "environmental", "unclassified", "sp", "sp."}
UNSAFE = re.compile(r"[^A-Za-z0-9_.\-]")
NUMSUFFIX = re.compile(r"_\d+$")                  # Index Fungorum ids on genus/species names
DEFAULT_SP_PRIOR = 7.5294e-6                      # per-species prior in PROTAX-fungi's taxonomy

def safe(x):
    x = UNSAFE.sub("_", (x or "").strip())
    return x or "NA"

def key(base):
    return NUMSUFFIX.sub("", base).lower()

def norm(x):
    return (x or "").strip().lower().replace(" ", "_")

def real_rank(x):
    v = norm(x)
    return v not in BAD and "incertae" not in v and not v.endswith("_sp") and not v.endswith("_sp.")

def real_species(x):
    v = norm(x)
    return v not in BAD and not v.endswith("_sp") and not v.endswith("_sp.")

# ---------------------------------------------------------------- input

def parse_fasta(path):
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
                seq.append(line.strip())
    if header is not None:
        yield header, "".join(seq)

ACC_RE = re.compile(r"^(UDB\d+|[A-Z]{1,3}\d{5,8}(\.\d+)?)$")
SH_RE = re.compile(r"^SH\d+\.\d+\w*$")

def parse_unite_header(header):
    """Same id and lineage rules as code/prep_split.py, so ids line up with the splits."""
    parts = [p.strip() for p in header.split("|")]
    ranks = {r: "" for r in RANKS}
    lineage = next((p for p in parts if "__" in p and ";" in p), "")
    for tok in lineage.split(";"):
        pref, _, val = tok.strip().partition("__")
        if pref in ranks:
            ranks[pref] = val.strip()
    sid = next((p for p in parts if ACC_RE.match(p)), None) \
        or next((p for p in parts if SH_RE.match(p)), None) \
        or next((p for p in parts if p and p != lineage), header.split()[0])
    return sid, ranks

def read_taxonomy_tsv(path):
    tax = {}
    with open(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        ix = {h: i for i, h in enumerate(head)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < len(head):
                continue
            tax[f[ix["seq_id"]]] = {r: f[ix[r]] for r in RANKS}
    return tax

def read_split(path):
    roles = {}
    with open(path) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2:
                roles[f[0]] = f[1]
    return roles

# ---------------------------------------------------------------- tree

class Tree:
    def __init__(self):
        self.pid, self.level, self.base, self.prior = {}, {}, {}, {}
        self.children = defaultdict(list)
        self.created = set()
        self.idx = defaultdict(lambda: defaultdict(list))   # level -> key(base) -> [nid]
        self.next_id = 0

    def add(self, nid, pid, level, base, prior, created=False):
        self.pid[nid], self.level[nid], self.base[nid], self.prior[nid] = pid, level, base, prior
        if nid != pid:
            self.children[pid].append(nid)
        if level > 0:
            self.idx[level][key(base)].append(nid)
        if created:
            self.created.add(nid)
        self.next_id = max(self.next_id, nid + 1)
        return nid

    def new(self, pid, level, base, suffix=False):
        nid = self.next_id
        return self.add(nid, pid, level, f"{base}_{nid}" if suffix else base, 0.0, created=True)

    def anc(self, nid, level):
        while self.level[nid] > level:
            nid = self.pid[nid]
        return nid

def load_protax_taxonomy(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            f = line.split()
            if not f:
                continue
            if len(f) < 5:
                sys.exit(f"[build] {path}: expected 5 columns (id pid level name prior). "
                         "Pass PROTAX-fungi's db_files/tax7 or taxonomy file.")
            rows.append((int(f[0]), int(f[1]), int(f[2]), f[3], float(f[4])))
    t = Tree()
    t.add(0, 0, 0, "root", 1.0)
    rows.sort(key=lambda r: (r[2], r[0]))
    for nid, pid, lev, name, prior in rows:
        if lev == 0:
            continue
        if name == "nonFungi" or name.startswith("nonFungi,"):
            continue                                # cptaxonomy.pl drops this branch too
        if pid not in t.level:
            continue                                # orphan under a dropped branch
        base = safe(name.split(",")[-1])
        if lev == 1 and base == "Fungi":
            prior = 1.0
        t.add(nid, pid, lev, base, prior)
    fungi = [n for n in t.idx[1].get("fungi", [])]
    if not fungi:
        sys.exit(f"[build] {path}: no level-1 node named Fungi")
    return t, fungi[0]

def new_tree():
    t = Tree()
    t.add(0, 0, 0, "root", 1.0)
    return t, t.add(1, 0, 1, "Fungi", 1.0)

# ---------------------------------------------------------------- name -> node mapping

class Mapper:
    def __init__(self, tree, fungi):
        self.t, self.fungi = tree, fungi
        self.rank_cache, self.genus_cache = {}, {}
        self.stats = Counter()

    def _score(self, nid, r, top):
        """How many of the lineage's names above `top` agree with the node's ancestors."""
        s = 0
        for lev in range(top, 1, -1):
            if key(self.t.base[self.t.anc(nid, lev)]) == key(safe(r[RANKS[lev - 2]])):
                s += 1
        return s

    def _pick(self, cands, r, level):
        if len(cands) == 1:
            return cands[0]
        best = max(self._score(c, r, level - 1) for c in cands)
        self.stats[f"ambiguous_level{level}"] += 1
        return min(c for c in cands if self._score(c, r, level - 1) == best)

    def rank_node(self, r, level):
        """Node for the lineage's name at level 2..5, matched or created."""
        ck = tuple(r[x] for x in RANKS[:level - 1])
        if ck in self.rank_cache:
            return self.rank_cache[ck]
        name = safe(r[RANKS[level - 2]])
        cands = self.t.idx[level].get(key(name), [])
        if cands:
            nid = self._pick(cands, r, level)
        else:
            parent = self.fungi if level == 2 else self.rank_node(r, level - 1)
            nid = self.t.new(parent, level, name)
            self.stats[f"created_level{level}"] += 1
        self.rank_cache[ck] = nid
        return nid

    def genus_node(self, r):
        ck = (r["p"], r["c"], r["o"], r["f"], r["g"])
        if ck in self.genus_cache:
            return self.genus_cache[ck]
        cands = self.t.idx[6].get(key(safe(r["g"])), [])
        if cands:
            nid = self._pick(cands, r, 6)
            self.stats["genus_matched"] += 1
        else:
            nid = self.t.new(self.rank_node(r, 5), 6, safe(r["g"]), suffix=True)
            self.stats["genus_created"] += 1
        self.genus_cache[ck] = nid
        return nid

    def species_node(self, gid, s):
        k = key(safe(s))
        for c in self.t.children[gid]:
            if self.t.level[c] == 7 and key(self.t.base[c]) == k:
                return c
        self.stats["species_created"] += 1
        return self.t.new(gid, 7, safe(s), suffix=True)

    def place(self, r):
        """(node, level) for a record: deepest real rank, species when it is a binomial."""
        if real_rank(r["g"]):
            g = self.genus_node(r)
            if real_species(r["s"]):
                return self.species_node(g, r["s"]), 7
            return g, 6
        deepest = 1
        for lev in range(2, 6):
            if real_rank(r[RANKS[lev - 2]]):
                deepest = lev
            else:
                break
        return (self.fungi, 1) if deepest == 1 else (self.rank_node(r, deepest), deepest)

# ---------------------------------------------------------------- priors, names, output

def finish_tree(t, sp_prior):
    # created nodes: species get sp_prior, others the sum of their children (min sp_prior)
    for nid in sorted(t.created, key=lambda n: -t.level[n]):
        if t.level[nid] == 7:
            t.prior[nid] = sp_prior
        else:
            t.prior[nid] = max(sp_prior, sum(t.prior[c] for c in t.children[nid]))
    # unique basenames per level (sintax2sa.pl maps rank names back to nodes by basename)
    for lev in range(1, 8):
        seen = defaultdict(list)
        for nid, l in t.level.items():
            if l == lev:
                seen[t.base[nid]].append(nid)
        for base, ids in seen.items():
            for nid in sorted(ids)[1:]:
                t.base[nid] = f"{base}_n{nid}"
    names = {0: "root"}
    for nid in sorted(t.level, key=lambda n: (t.level[n], n)):
        if t.level[nid] == 1:
            names[nid] = t.base[nid]
        elif t.level[nid] > 1:
            names[nid] = names[t.pid[nid]] + "," + t.base[nid]
    return names

def write_model(out, t, names, refs, placed, region, params_dir, fungi):
    os.makedirs(out, exist_ok=True)
    order = sorted(t.level, key=lambda n: (t.level[n], n))
    line = lambda n: f"{n}\t{t.pid[n]}\t{t.level[n]}\t{names[n]}\t{t.prior[n]:.10g}\n"
    for fname in ("taxonomy", "taxonomy.ascii7"):
        with open(os.path.join(out, fname), "w") as fh:
            fh.writelines(line(n) for n in order)
    for lev in range(1, 8):
        with open(os.path.join(out, f"tax{lev}"), "w") as fh:
            fh.writelines(line(n) for n in order if t.level[n] <= lev)
    with open(os.path.join(out, "seqid2tax"), "w") as fh:
        for sid, (node, lev) in placed.items():
            fh.write(f"{sid}\t{lev}\t{names[node]}\n")
    rseqs = {lev: defaultdict(list) for lev in range(2, 8)}
    for lev in range(2, 8):
        with open(os.path.join(out, f"ref.tax{lev}"), "w") as fh:
            for sid, (node, slev) in placed.items():
                if slev >= lev:
                    a = t.anc(node, lev)
                    fh.write(f"{sid}\t{names[a]}\n")
                    rseqs[lev][a].append(sid)
        with open(os.path.join(out, f"rseqs{lev}"), "w") as fh:
            for a in sorted(rseqs[lev]):
                fh.write(f"{a}\t{','.join(rseqs[lev][a])}\n")
    with open(os.path.join(out, f"{region}.fa"), "w") as fa, \
         open(os.path.join(out, f"sintax{region}train.fa"), "w") as st:
        for sid, (node, slev) in placed.items():
            fa.write(f">{sid}\n{refs[sid]}\n")
            tax = ",".join(f"{SINTAX[l]}:{t.base[t.anc(node, l)]}" for l in range(1, slev + 1))
            st.write(f">{sid};tax={tax};\n{refs[sid]}\n")
    for lev in range(2, 8):
        src = os.path.join(params_dir, f"params.{region}.level{lev}")
        if not os.path.exists(src):
            sys.exit(f"[build] missing {src}; point --params at protaxfungi/model")
        shutil.copy(src, out)
    with open(os.path.join(out, "startnode"), "w") as fh:
        fh.write(f"{fungi}\n")
    return rseqs

# ---------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--split-dir", help="prep_split.py output: refs.fasta, taxonomy.tsv, split.tsv")
    src.add_argument("--unite-fasta", help="every record is a reference (e.g. the full ITS2 pool)")
    ap.add_argument("--protax-taxonomy", help="PROTAX-fungi taxonomy with priors (db_files/tax7)")
    ap.add_argument("--params", required=True, help="dir with params.<region>.level2..7")
    ap.add_argument("--region", default="its2", choices=["its1", "its2", "itsfull"])
    ap.add_argument("--species-prior", type=float, default=None,
                    help="prior for created species (default: modal species prior of the "
                         f"PROTAX taxonomy, else {DEFAULT_SP_PRIOR:g})")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    # references and the names they carry
    if a.split_dir:
        refs = {h.split()[0]: s.upper() for h, s in parse_fasta(os.path.join(a.split_dir, "refs.fasta"))}
        tax = read_taxonomy_tsv(os.path.join(a.split_dir, "taxonomy.tsv"))
        roles = read_split(os.path.join(a.split_dir, "split.tsv"))
    else:
        refs, tax, roles = {}, {}, {}
        for h, s in parse_fasta(a.unite_fasta):
            sid, r = parse_unite_header(h)
            refs[sid], tax[sid] = s.upper(), r
    missing = [s for s in refs if s not in tax]
    if missing:
        sys.exit(f"[build] {len(missing)} reference ids have no taxonomy row, e.g. {missing[:3]}")

    if a.protax_taxonomy:
        t, fungi = load_protax_taxonomy(a.protax_taxonomy)
        sp7 = Counter(round(p, 14) for n, p in t.prior.items() if t.level[n] == 7)
        sp_prior = a.species_prior or (sp7.most_common(1)[0][0] if sp7 else DEFAULT_SP_PRIOR)
        n_official = len(t.level)
    else:
        print("[build] WARNING: no --protax-taxonomy; building the tree from the input's own "
              "names. Fine for testing, biased towards 'known' for the comparison.")
        t, fungi = new_tree()
        sp_prior = a.species_prior or DEFAULT_SP_PRIOR
        n_official = 0
    official_genera = {n for n in t.level if t.level[n] == 6}

    m = Mapper(t, fungi)
    placed_all = {sid: m.place(r) for sid, r in tax.items()}     # queries too: their names
    placed = {sid: placed_all[sid] for sid in refs}              # must exist as nodes
    names = finish_tree(t, sp_prior)
    rseqs = write_model(a.out, t, names, refs, placed, a.region, a.params, fungi)

    # ---- checks and report
    with_refs6 = set(rseqs[6])
    rep = [f"source            : {a.split_dir or a.unite_fasta}",
           f"taxonomy          : {a.protax_taxonomy or '(built from input names)'}",
           f"nodes             : {len(t.level)} ({n_official} from PROTAX taxonomy, "
           f"{len(t.created)} created)",
           f"species prior     : {sp_prior:.6g}",
           f"reference seqs    : {len(placed)} (level 7: "
           f"{sum(1 for v in placed.values() if v[1] == 7)}, level 6: "
           f"{sum(1 for v in placed.values() if v[1] == 6)}, above genus: "
           f"{sum(1 for v in placed.values() if v[1] < 6)})",
           f"genera with refs  : {len(with_refs6)} of {sum(1 for n in t.level if t.level[n] == 6)}"]
    rep += [f"mapping           : {k} = {v}" for k, v in sorted(m.stats.items())]
    status = 0
    if roles:
        q = [s for s, role in roles.items() if role == "query" and s in placed_all]
        qg = {s: t.anc(placed_all[s][0], 6) for s in q if placed_all[s][1] >= 6}
        n_ref = sum(1 for g in qg.values() if g in with_refs6)
        held = {g for g in qg.values() if g not in with_refs6}
        fam_noref = [sum(1 for c in t.children[t.pid[g]] if t.level[c] == 6 and c not in with_refs6)
                     for g in held]
        rep += [f"queries           : {len(q)}",
                f"  genus has refs  : {n_ref}  (expected all for LSO, none for LGO)",
                f"  genus has none  : {len(qg) - n_ref}",
                f"held-out genera   : {len(held)} ({sum(1 for g in held if g in official_genera)} "
                f"named in the PROTAX taxonomy, the rest created)"]
        if fam_noref:
            fam_noref.sort()
            rep.append(f"  named genera without refs in their family (incl. itself): "
                       f"median {fam_noref[len(fam_noref) // 2]}, max {fam_noref[-1]}")
        if 0 < n_ref < len(qg):
            rep.append("  WARNING: mixed: some query genera have references and some do not. "
                       "Expected for neither LGO nor LSO; check name collisions above.")
            status = 2
    with open(os.path.join(a.out, "build_report.txt"), "w") as fh:
        fh.write("\n".join(rep) + "\n")
    print("\n".join(f"[build] {x}" for x in rep))
    print(f"[build] wrote model to {a.out}/")
    sys.exit(status)

if __name__ == "__main__":
    main()
