# A calibrated novelty flag for fungal ITS metabarcoding

Harness, verdict tables and figure code for the conformal novelty flag manuscript.

Each ITS2 query receives a p-value with a distribution-free bound on the rate at which
known sequences are wrongly declared novel. The bound holds across four scores of
different construction, two amplicon regions and a sevenfold change in reference size;
detection over the same settings varies more than tenfold, which is the point of
reporting both.

## Layout

```
code/       the harness, in pipeline order below
data/       verdict tables, p-values, ASV set, named-arm recurrence hits
taxa/       MycoAI label space, for the cross-score comparison strata
figures/    fig2 calibration, fig3 gradient, fig4 dark-matter fraction
logs/       cutadapt, DADA2 and vsearch run logs, as provenance
```

## What "flagged" means, precisely

An ASV is flagged novel when its conformal p-value is at most alpha. The identity cutoff
the procedure derives, 64.8% at alpha = 0.05 on this reference, is reported in the paper
so the requested error rate can be read on a familiar scale, and is used as a threshold
nowhere. `code/make_flagged_fasta.py` regenerates the flagged set from the p-value column
so the two definitions cannot drift apart.

## Not included, and why

- **The GlobalFungi variant library** (131M and 50M records) and the search output against
  it (~1.3 GB). Regenerable by `code/gf_recurrence.sh` from public data. The named arm's
  hits are small and are included.
- **UNITE and SILVA releases.** Cited, not redistributed. UNITE general FASTA release for
  Fungi version 19.02.2025, doi:10.15156/BIO/3301229.
- **Raw reads and the DADA2 intermediate objects.** NCBI BioProject PRJNA704912;
  `code/dada2_its2.R` reproduces the sequence table.
- **The manuscript source.** Added at submission.

## Pipeline

```
# 1. reference, splits, calibration
python code/prep_split.py --input unite.ITS2.fasta --out out_lso_its2 --mode lso
python code/score_recovery.py ...
python code/conformal_novelty.py --known lso.hits --novel lgo.hits \
    --alphas 0.01,0.05,0.10,0.20

# 2. soil application
Rscript code/dada2_its2.R                 # sequence table from PRJNA704912
Rscript code/dada2_half.R                 # split-half determinism check
python  code/score_asvs_conformal.py      # per-ASV p-values and verdicts
Rscript code/asv_abundance.R              # prevalence and read totals per ASV
python  code/make_flagged_fasta.py        # flagged set, from the p-value column

# 3. cross-study recurrence, both arms, against all 5,521 ASVs
bash code/gf_recurrence.sh                # ARM=named first: seconds, checks the call

# 4. species-hypothesis prevalence comparison
python code/frey_sh_verdicts.py
python code/gf_prevalence.py

# 5. figures
python code/make_figures.py
```

`code/pick_studies.py` selected the application dataset and is included for completeness.

## Three things a reader may trip over

**Three routes to ITS2, not one.** The reference and the benchmark queries drawn from it
are full-length records and were delimited with ITSx. The soil ASVs and the GlobalFungi
variants could not be, and were trimmed with cutadapt. All three end at ITS2; they are
not the same operation, and residual boundary disagreement is absorbed by scoring
coverage over the shorter sequence rather than requiring mutual full-length identity.

**The two recurrence arms are not symmetric.** The unnamed arm searches ~50.3M individual
variants, the named arm 26,096 species-hypothesis representatives. A higher rate in the
named arm does not mean a closer relationship. Each arm is comparable within itself,
across ASV groups, which is what the table uses it for.

**`--maxaccepts 20 --maxrejects 0` is deliberate.** With vsearch's default single accepted
hit, a global variant matching both a flagged and a non-flagged ASV is reported against
only one of them, which would bias the comparison between those groups mechanically.

## Environment

vsearch v2.31.0, cutadapt v4.4, DADA2 v1.38.0 on R 4.5.3.

## License

MIT, see LICENSE.
