# PROTAX-fungi against the conformal absence call

The editor's condition for resubmission: show how the method relates to PROTAX-fungi and
whether it adds anything. This harness runs PROTAX-fungi on the same held-out queries and
the same reference sequences as the identity scores, then measures both side by side.

## Setup (once)

1. **PROTAX-fungi v1** (Perl scripts and per-rank parameters):
   `curl -LO https://raw.githubusercontent.com/psomervuo/protaxfungi/master/protaxfungi.tgz && tar xzf protaxfungi.tgz`
2. **PROTAX-fungi's taxonomy, with priors**, as UNITE deploys it: run `run_setup.sh` from
   [TU-NHM/protax_fungi_plutof_pub](https://github.com/TU-NHM/protax_fungi_plutof_pub),
   which downloads `db_files/` from PlutoF. The builder needs `db_files/tax7` (or
   `db_files/taxonomy`); both carry the fifth, prior column. If `db_files/` also holds
   `params.its2.level2..7`, pass `--params db_files` so the parameters are the ones PlutoF
   runs; otherwise `--params protaxfungi/model`. Record which in the methods.
3. vsearch (2.31.0 here), Perl 5, Python 3 with numpy; matplotlib for the figure.

Check the wiring before the long run (about a minute, synthetic data):

```
PROTAX_SCRIPTS=protaxfungi/protaxscripts PROTAX_PARAMS=protaxfungi/model \
  bash code/protax/tests/selftest.sh /tmp/protax_selftest
```

## Run

```
P=code/protax
export PROTAX_SCRIPTS=protaxfungi/protaxscripts THREADS=16     # CHUNKS defaults to THREADS
for mode in lgo lso; do
  python $P/build_protax_model.py --split-dir out_${mode}_its2 \
      --protax-taxonomy db_files/tax7 --params db_files --out model_${mode}_its2
  bash $P/run_protax.sh model_${mode}_its2 out_${mode}_its2/queries.fasta protax_${mode}
  python $P/protax_scores.py --protax-dir protax_${mode} --level 6 \
      --taxonomy model_${mode}_its2/taxonomy --out protax_${mode}/scores_genus.tsv
done
python $P/compare_protax_conformal.py \
  --known-queries out_lso_its2/queries.fasta --novel-queries out_lgo_its2/queries.fasta \
  --known-hits lso_its2.hits --novel-hits lgo_its2.hits \
  --known-protax protax_lso/scores_genus.tsv --novel-protax protax_lgo/scores_genus.tsv \
  --out-prefix results/protax_vs_conformal
```

Add `--hits-only` (and a different `--out-prefix`) to reproduce Table 1's query set.

Soil ASVs: build the model from the same reference the ASV identity hits were computed
against (`--unite-fasta` for a full pool, `--split-dir` for a split's refs), then

```
bash $P/run_protax.sh model_asvs dada2_out/asvs_its2.fasta protax_asvs
python $P/protax_scores.py --protax-dir protax_asvs --out protax_asvs/scores_genus.tsv
python $P/crosstab_soil.py --flags data/frey_asv_pvalues.tsv \
  --protax protax_asvs/scores_genus.tsv --out results/soil_crosstab.tsv
```

## Reading the outputs

- `model_*/build_report.txt`: "genus has refs" must be all queries for LSO and none for
  LGO; anything else means a name collision to fix before running. It also reports how
  many held-out genera PROTAX's taxonomy names, and how many named genera without
  sequences sit in their families, which is what PROTAX's absence probability competes with.
- `protax_*/scores_genus.tsv`: `p_absent` is PROTAX's probability that the query's genus
  has no reference sequences (named or not); `truncated` is mass classify.pl did not print.
- `results/*_summary.txt`: the five measurements. `*_mixtures.pdf`: realised FDR and
  detection against the share of held-out-genus queries, the figure for the response.

## Choices to state in the methods

- **Parameters** are PROTAX-fungi's published per-rank ITS2 parameters, not refitted. They
  were fitted on UNITE including the held-out genera; they are per rank, not per taxon, so
  the leakage should be small, and if anything it favours PROTAX.
- **Taxonomy and priors** are PROTAX-fungi's own. This matters: the prior enters every
  "known" branch's score, and named genera without sequences compete with the known ones.
  In the synthetic self-test, moving the priors from a small tree's scale to
  PROTAX-fungi's per-species scale moved held-out-genus queries' mean absence probability
  from 0.04 to 0.50. The builder's fallback without `--protax-taxonomy` is for testing only.
- **Held-out taxa** stay in the tree as named taxa without references, exactly where
  PROTAX puts any named-but-unsequenced taxon.
- **Search** follows UNITE's PlutoF deployment: vsearch `--usearch_global --id 0.75
  --maxaccepts 1000 --strand both` and vsearch `--sintax`. `PROB_TH` is 1e-4 (PlutoF 0.01);
  it only truncates the printed list, not the probabilities. SINTAX gets a fixed seed
  (`SINTAX_SEED`, default 1): vsearch otherwise seeds it randomly, as on PlutoF, so PROTAX's
  probabilities would shift slightly between runs.
- **Queries without an identity hit** at the 50% floor are kept at identity 0, in
  calibration and evaluation. `conformal_novelty.py` dropped them from both.
- **References**: known (LSO) queries are scored against the LSO split's reference and
  held-out-genus (LGO) queries against the LGO split's, for both methods, as in the paper.

classify.pl is single-threaded and visits every reference sequence of each candidate taxon
at each rank, so split the work: `CHUNKS` (default `THREADS`) classifies query chunks in
parallel and concatenates them in input order.
