# Identified, absent or unresolved

Code, figures and soil data for *Identified, absent or unresolved: error-controlled
identification of fungal ITS2 sequences against a reference database* (O'Brien, Marín and
Parada).

Conformal calibration turns a similarity or classifier score into calls with an error rate
the user chooses:

- the **absence call** declares a sequence's genus absent from the reference, while
  sequences of genera in the reference are wrongly called absent at most α of the time;
- the **three-level calls** report a sequence as identified, absent or unresolved among
  named candidates, at species, genus and family rank, and keep the chance that a sequence
  of a represented taxon gets a wrong call of any kind within one budget α.

"Absent" always means absent from the reference used, at the rank considered. Whether such
a sequence belongs to an undescribed taxon is a question for phylogenetic analysis and
nomenclature, not for these calls.

The first version of this work, *A calibrated novelty flag for fungal ITS metabarcoding*
(bioRxiv, doi:10.64898/2026.07.29.741524), is tagged `v0.1.0`. This version is `v0.2.0`.

## Layout

    code/                     splits, searches, soil application, Figures 1 and 4
    code/protax/              PROTAX-fungi runs and the comparison (see its README)
    code/tiers/               three-level calls and genus pairs (see its README)
    code/refresh_numbers.sh   regenerates every number in the paper
    data/                     soil ASVs: ITS2 from ITSxRust, ASVs without an ITS2 call,
                              read counts; tables from the first version
    figures/                  Figures 1 to 4 of the paper
    logs/, taxa/              first-version run logs and MycoAI label space

Scripts from the first version stay in `code/`; their pipeline is documented at `v0.1.0`.

## What "absent" means, precisely

A query's genus is called absent when its conformal p-value is at most α. The p-value
ranks the query's best-hit identity among calibration queries whose genus is in the
reference but whose own species was removed. The identity cutoff this implies, 69.1% at
α = 0.05 (median over 200 calibration splits), is reported in the paper to put the error
rate on a familiar scale; the calls themselves come from the p-values.

## Regenerating the numbers

    bash code/refresh_numbers.sh

runs the searches, the PROTAX-fungi comparison, the three-level calls, the genus pairs and
the soil application, and writes the results, with Figures 2 and 3, to `refresh/results`
under the current folder. It takes about half an hour on 32 threads. Searches are cached,
so a rerun only repeats the analyses. The held-out-species and held-out-genus splits are
required; every later step is skipped, with a message, when its inputs are missing.
Inputs are set by environment variables; apart from `SOIL`, the defaults are the author's
run folders.

| variable | points to | made by |
|---|---|---|
| `SPLITS` | `out_lso_its2/`, `out_lgo_its2/` | `code/prep_split.py --mode lso`, `--mode lgo` |
| `NEWSPLITS` | `out_loo_its2/`, `out_lfo_its2/` | `code/prep_split.py --mode loo`, `--mode lfo` |
| `PROTAX` | `protax_lso/`, `protax_lgo/`, each with `scores_genus.tsv` | `code/protax/` (see its README) |
| `SOIL` | ITS2 of the soil ASVs, default `data/asvs_its2_itsx.fasta` | ITSxRust v0.3.0, `itsxrust extract --region its2` |
| `GFDIR` | GlobalFungi search output (`gf_*_vs_all.hits`) and `dada2_out/` | `code/gf_recurrence.sh`, `code/dada2_its2.R` |
| `THREADS`, `OUT` | threads (32) and output folder (`refresh`) | |

Each split folder holds `queries.fasta`, `refs.fasta` and `taxonomy.tsv`, built with a
fixed seed from the ITS2 of the UNITE general release 19.02.2025 (extracted with ITSx):

    python code/prep_split.py --input unite.ITS2.fasta --out out_lso_its2 --mode lso

and likewise with `--mode lgo`, `loo` and `lfo`. Table S2 comes from the same script run
on splits built from the 20,144 ITS2 sequences that UNITE designates as reference
sequences.

Figures 1 and 4 are drawn from the results:

    python code/make_overview_figure.py --perquery refresh/results/protax_vs_conformal_perquery.tsv \
        --conformal refresh/results/protax_vs_conformal_conformal.tsv --out fig_overview
    python code/make_soil_figure.py --kept data/asvs_its2_itsx.fasta --skipped data/asvs_itsx_skipped.fasta \
        --abundance data/asv_abundance.tsv --calls refresh/results/soil_itsx_a0.10_apply_a0.1.tsv \
        --table3 refresh/results/table3.txt --out fig_soil

## Things a reader may trip over

**The search is exhaustive.** vsearch runs with `--maxaccepts 0 --maxrejects 0` and keeps
an alignment only if it covers at least 80% of the shorter sequence; queries without a
qualifying hit are kept with identity 0. Stopping at the first hit above the floor, as the
first version did, returned a worse match than the best one for 34% of held-out-species
queries and 60% of held-out-genus queries.

**For the calls, ITS2 is cut the same way on both sides.** The reference was delimited
with ITSx and the soil ASVs with ITSxRust using ITSx's fungal profiles, which also sets
aside amplicons with no fungal 5.8S end, plant ITS2 among them.

**The guarantee is per query, not per list.** The share of wrong calls in a list of
absence calls depends on how common absent taxa are in the sample (Figure 2B).

**The two recurrence arms are not symmetric.** The unnamed arm searches about 50.3 million
individual GlobalFungi variants, the named arm 26,096 species-hypothesis representatives.
Each arm is comparable within itself, across ASV groups, which is how Table 3 uses it.
The soil ASVs and GlobalFungi delimit ITS2 differently, so coverage is measured on the
GlobalFungi sequence.

**`--maxaccepts 20 --maxrejects 0` in the GlobalFungi search is deliberate.** With
vsearch's default of one accepted hit, a variant matching two ASVs would be reported
against only one of them, which would bias the comparison between groups.

## Not included, and why

- **The GlobalFungi variant library** (figshare, doi:10.6084/m9.figshare.c.4915392,
  `GlobalFungi_ITS_variants.zip` as deposited on 5 April 2020) and the search output
  against it (1.3 GB), regenerated by `code/gf_recurrence.sh`. The named arm's hits are
  small and are in `data/`.
- **UNITE, PROTAX-fungi and UNITE's PlutoF files.** Cited, not redistributed. UNITE
  general FASTA release for Fungi version 19.02.2025, doi:10.15156/BIO/3301229;
  `code/protax/README.md` says where to get the other two.
- **Raw reads and DADA2 intermediates.** NCBI BioProject PRJNA704912;
  `code/dada2_its2.R` reproduces the sequence table.
- **The manuscript.**

## Environment

vsearch v2.31.0, ITSx (reference), ITSxRust v0.3.0 (soil ASVs), PROTAX-fungi v1 with
Perl 5, Python 3 with numpy and matplotlib, DADA2 v1.38.0 on R 4.5.3, cutadapt v4.4
(GlobalFungi recurrence).

## License

MIT, see LICENSE.
