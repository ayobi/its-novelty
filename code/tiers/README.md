# Three-level calls: identified, absent, or unresolved

Two conformal steps share one error budget α (half each by default, `--absence-share`):

- **absent**: the absence test on best-hit identity, at α/2.
- **otherwise**: the set of reference taxa within *d* identity points of the query's best hit.
  One taxon means identified, two or more unresolved (the set names the candidates). *d* is
  calibrated on known queries that pass the absence test, so their true taxon is missed at
  most α/2 of the time.

For known queries every error (called absent, identified as another taxon, unresolved without
the true taxon) then stays within α. Sets follow each query's best hit; a single identity
cutoff for all queries had to fall to 59% at genus rank and put hundreds of genera in every set.

| rank    | known queries (calibration)           | held-out queries          |
|---------|---------------------------------------|---------------------------|
| species | `--mode loo` (species in reference)   | `--mode lso`              |
| genus   | `--mode lso` (genus in reference)     | `--mode lgo`              |
| family  | `--mode lgo` (family in reference)    | `--mode lfo`              |

`prep_split.py` gained `lfo` and `loo`; `lgo` and `lso` are unchanged (outputs byte-identical).

## Per split and rank

```
vsearch --usearch_global out_X/queries.fasta --db out_X/refs.fasta --id 0.5 \
    --maxaccepts 0 --maxrejects 0 --threads 32 --quiet \
    --userfields query+target+id+qcov+tcov+ql+tl --userout /dev/stdout \
| python3 code/tiers/best_per_taxon.py --taxonomy out_X/taxonomy.tsv --rank g \
    --queries out_X/queries.fasta --out X_genus_best.tsv
```

The search must be exhaustive (`--maxaccepts 0 --maxrejects 0`): a capped search missed
about 9% of the genera a query matched at 60% identity or more. Hits must also cover at
least 80% of the shorter sequence (`--min-cover`), because vsearch's identity ignores
terminal gaps and an exhaustive search otherwise picks up short end overlaps.

## Calls and evaluation

```
python3 code/tiers/conformal_tiers.py --rank g \
    --known-best lso_genus_best.tsv --known-taxonomy out_lso_its2/taxonomy.tsv \
    --heldout-best lgo_genus_best.tsv --out-prefix results/tiers_genus
```

Writes `_summary.txt` (tier shares for known and held-out queries at each α, the single
absence test at the full α for comparison, and the share of wrong calls by how common absent taxa
are), `_tiers.tsv`, `_mixtures.tsv`, `_perquery.tsv` (one calibration split, with candidate
lists) and `_tiers.pdf/.png`. Add `--apply-best new_best.tsv` to label new queries such as
soil ASVs, calibrated on all known queries.

Self-test (synthetic data, all three ranks): `bash code/tiers/tests/selftest_tiers.sh`
