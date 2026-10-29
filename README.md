# processor-diffcyt

Which immune cell populations change between conditions, and which change what they express?
This step of the **Immune Profiling & Differential Abundance** workflow answers both with
**diffcyt** (Weber et al. 2019), following the Nowicka et al. (2017) CyTOF workflow:

- **Differential abundance (DA)**: is a metacluster (a cell population found by FlowSOM) more
  or less frequent in one condition than in the reference? E.g. Tregs expand after 7 days of
  anti-CD3/CD28 + IL-2 stimulation.
- **Differential state (DS)**: within a metacluster, does a functional or activation marker
  (CD25, HLA-DR, CD69, PD-1, CD39, ...) shift between conditions? The test is on each sample's
  median expression in that cluster, so it is about the population, not single cells.

Both are tests across **samples**, not cells: each sample (donor, tube) is one observation, the
variability between samples within a condition is what a change is judged against, and a
batch (or donor) can be accounted for. That is the point of diffcyt over per-cell tests,
which treat a million cells as a million replicates.

The app runs after `processor-flowsom` and reads its per-sample tables (cell counts and median
expression per metacluster) rather than per-cell data, so it is light: a couple of seconds of
model fitting after R has loaded, in under 1 GB, as a container or a Lambda function.

`workflow: fcs-prep -> peacoqc -> cleanet -> flowsom -> {cytometry-umap, diffcyt} -> cytometry-report`

## Ports

| Port | Kind | Files | Notes |
|---|---|---|---|
| `clusters` | input (`placement: subfolder`, required) | `counts.csv` (sample_id, metacluster, label, n, frac), `medians_by_sample.csv` (sample_id, metacluster, marker, median), `labels.csv` (metacluster, label) | From flowsom. Medians are on the arcsinh scale of fcs-prep. `medians_by_sample.csv` is needed for DS only; labels come from `counts.csv`'s `label` column, else `labels.csv`, else the id. `frac` is not read: proportions are recomputed from `n`. |
| `metadata` | input (`placement: subfolder`, required) | `samples.csv` (sample_id, batch, condition, ...), `panel.csv` (marker, type, clustering, state) | From fcs-prep. `panel.csv` supplies the default state markers (`state == true`). |
| `da` | output | `da/da.csv`, `da/volcano.png`, `da/proportions.png`, `da/summary.json`, `da/provenance.json` | |
| `ds` | output | `ds/ds.csv`, `ds/heatmap.png`, `ds/logfc.png`, `ds/summary.json`, `ds/provenance.json` | |

Linked from a workflow, the inputs arrive in `INPUT_DIR/clusters/` and `INPUT_DIR/metadata/`.
Without links (everything side by side under `INPUT_DIR`, possibly nested, entries symlinked)
the files are found by name anywhere under `INPUT_DIR`.

## Parameters

| Parameter | Type | Default | What it does |
|---|---|---|---|
| `CONDITION_COLUMN` | string | `condition` | Column of `samples.csv` with the groups to compare. |
| `REFERENCE_LEVEL` | string | *(empty)* | The level the others are compared with; logFC > 0 means higher than in the reference. Empty: the first level alphabetically (for the demo that would be `stimulated`, so set `unstimulated`). |
| `BLOCK_COLUMN` | string | `batch` | Column with a blocking factor (batch, donor). **Fixed effect** in the edgeR / voom / limma design (`~ condition + block`), **random intercept** for GLMM / LMM (`(1 \| block)`). Empty (or `none`): no block. Dropped with a warning (in the log and `provenance.json`) when it has one level, is confounded with the condition, or leaves no residual degrees of freedom. |
| `DA_METHOD` | `edgeR` \| `voom` \| `GLMM` | `edgeR` | DA test: `testDA_edgeR` (negative binomial GLM, likelihood ratio test; diffcyt's recommendation), `testDA_voom` (limma-voom) or `testDA_GLMM` (binomial GLMM with an observation-level random effect `(1 \| sample_id)` for overdispersion). |
| `DS_METHOD` | `limma` \| `LMM` | `limma` | DS test: `testDS_limma` (moderated linear model on per-sample medians, weighted by the cluster's cell counts, with a mean-variance trend; diffcyt's recommendation) or `testDS_LMM` (linear mixed model, weighted by each sample's total cells). |
| `MIN_CELLS` | integer | `3` | A cluster is tested only when it has at least `MIN_CELLS` cells in at least `MIN_SAMPLES` samples (diffcyt's filter). Filtered clusters stay in the tables with `tested = false` and empty p-values. |
| `MIN_SAMPLES` | number | `0` | See `MIN_CELLS`. `0` / empty: half the number of samples in the smallest condition level (6 for the demo's 12 vs 12; diffcyt's own default is half of all samples). |
| `FDR` | number | `0.05` | Benjamini-Hochberg adjusted p-value below which a result is `significant`. |
| `STATE_MARKERS` | string | *(empty)* | Comma-separated markers for DS (names as in `panel.csv`, matched ignoring case; unknown names are reported in `ds/summary.json` and skipped). Empty: the panel's `state == true` markers that have medians. |

Parameters are read case-insensitively where they are choices (`edger` works).

**More than two condition levels.** One model is fitted to all included samples, and each
non-reference level is tested against the reference with its own contrast
(`<level> vs <reference>`). The tables then have one row per cluster (and marker) **per
comparison**, told apart by the `comparison` column; Benjamini-Hochberg is applied within each
comparison, not across them. Mean columns are written for every level. The volcano and
log-fold-change figures get one panel per comparison. Comparisons between two non-reference
levels are not run: change `REFERENCE_LEVEL` for those.

**Samples.** A sample is used when it is in both `counts.csv` and `samples.csv`, has cells, has
a condition (empty and `NA` count as none) and, with a block, a block value. Levels with fewer
than 2 such samples are dropped. Each sample's status and the reason it was excluded are in
`provenance.json` and the log (one line per sample). The step fails when fewer than two levels
remain.

## Outputs

`da/da.csv`: one row per metacluster and comparison.

| Column | |
|---|---|
| `metacluster`, `label` | FlowSOM metacluster id and its label |
| `comparison` | `<level> vs <reference>` |
| `logFC` | log2 fold change of the cluster's abundance, level vs reference. edgeR / voom: the model coefficient (batch-adjusted). GLMM: diffcyt's GLMM gives only p-values, so this is the descriptive `log2(mean_frac level / mean_frac reference)` (`summary.json` says which). |
| `logCPM` | average log2 counts per million (cells in the cluster per million cells of the sample). edgeR: its own; voom / GLMM: `edgeR::aveLogCPM`. |
| `p_val`, `p_adj` | p-value and Benjamini-Hochberg adjusted p-value (empty when not tested) |
| `tested`, `significant` | `true`/`false`: passed the `MIN_CELLS` filter; `p_adj < FDR` |
| `mean_frac_<level>` | mean over the level's samples of the cluster's fraction of the sample's cells (0-1), one column per level |

`da/volcano.png` (logFC vs -log10 p, significant clusters coloured and labelled),
`da/proportions.png` (per metacluster, % of cells per condition with a point per sample, point
shape = block; panels ordered by adjusted p-value), `da/summary.json` (method, block and how
it entered the model, `min_samples` used, per comparison: clusters tested / significant / up /
down and the top hits).

`ds/ds.csv`: one row per metacluster, state marker and comparison: `metacluster`, `label`,
`marker`, `comparison`, `logFC`, `p_val`, `p_adj`, `tested`, `significant`,
`mean_median_<level>`. `logFC` here is, as in diffcyt, a **difference** in median arcsinh
expression (level minus reference): the limma coefficient (batch-adjusted), or for LMM the
descriptive difference of the group means of the medians.

`ds/heatmap.png`: the significant (cluster, marker) pairs (up to 60, by adjusted p-value) x
samples, the per-sample medians z-scored per row, samples ordered by condition then block,
with a condition strip on top. When nothing is significant it says so instead.
`ds/logfc.png`: metaclusters x state markers coloured by logFC, significant pairs ringed,
untested clusters gray. `ds/summary.json`: method, state markers and where they came from,
per comparison: pairs tested / significant / up / down, markers with hits, top hits.

Each output folder has a `provenance.json`: app name and version, R and package versions,
the parameters as given and as resolved (reference level, block used and how, `min_samples`),
the input files found (relative to `INPUT_DIR`), timings, and per-sample status. When one test
fails (e.g. no `medians_by_sample.csv`, so no DS) its `summary.json` and `provenance.json` say
`"status": "failed"` with the error; the step only fails when both do, or when the inputs are
missing or unusable (no `counts.csv` / `samples.csv`, unknown `REFERENCE_LEVEL` or
`BLOCK_COLUMN`, a bad parameter value).

## Methods in brief

`d_counts` (a `SummarizedExperiment`, clusters x samples, assay `counts`) and `d_medians` (one
assay per marker, clusters x samples, with `id_type_markers` / `id_state_markers` in its
metadata) are built directly from flowsom's tables in the exact shape diffcyt's `calcCounts()`
and `calcMedians()` would produce from per-cell data, and passed to the `testDA_*` /
`testDS_*` functions unchanged. A cluster with no cells in a sample has count 0 and median NA
there. Normalisation for composition (`normalize = TRUE`, TMM) is not used, as in diffcyt's
defaults: with a strong expansion (Tregs 4x) the other clusters' proportions drop a little,
and a very sensitive test (GLMM on large cell counts) can call those small shifts.

## Running it locally

```sh
make test       # Lambda handler tests (Python) + R end-to-end tests on the synthetic fixture
make run        # builds the image (amd64) and runs it over the fixture into data/output/
make run INPUT=./data/real   # or over real flowsom + prep outputs (see sample-input/README.md)
```

`make test` runs the R tests with a local `Rscript` when diffcyt is installed, otherwise inside
the image (`docker run ... Rscript tests/test_main.R`). The R tests check that the planted DA
(Treg 4x up, B cells 3x down) and DS (CD25 in CD4 memory, HLA-DR in CD8 memory) effects are
found with every method, both input layouts, no block, three condition levels, excluded
samples, a missing medians file and failing inputs.

## Runtime

`computeTypes: [standard, lambda]`, 1 vCPU / 2 GB, 15 min timeout. One image serves both:
`processor/entry.py` runs `Rscript R/main.R` as a container, or, when
`AWS_LAMBDA_RUNTIME_API` is set, starts the Lambda runtime interface client with
`processor/handler.py`, which turns the payload's `inputDir`, `outputDir` and parameters into
the same environment variables and runs the same script as a subprocess (non-zero exit -> the
invocation raises, so the step fails). R writes only to `OUTPUT_DIR` and `/tmp` (`HOME` and
`TMPDIR` point there), so the read-only Lambda filesystem is fine; checked with
`docker run --read-only --tmpfs /tmp` through the handler.

Measured with `docker run --platform linux/amd64 --cpus 1` (amd64 emulation on an M-series
Mac), cgroup `memory.peak`, which includes page cache:

| Data | DA / DS | wall time | peak memory |
|---|---|---|---|
| demo study (flowsom output): 24 samples, 20 metaclusters, 10 state markers, 10.2 M cells upstream | edgeR / limma | 7.7 s (R start-up and package loading ~5 s) | 1.12 GB |
| demo study | voom / limma | 7.0 s (log) | 1.03 GB |
| demo study | GLMM / LMM | 10.3 s (log) | 1.07 GB |
| synthetic fixture: 24 samples, 12 metaclusters, 9 state markers | edgeR / limma | 7.3 s | 0.94 GB; also completes within a 768 MB limit |

The peaks are cgroup `memory.peak` and are mostly page cache from reading ~150 R packages;
anonymous memory at the end of a run was a few MB, and the fixture run completes under a 768 MB
limit. 1 vCPU / 2 GB is therefore generous; it also suits Lambda (2 GB gets a little over one
vCPU there).

The model fits are on clusters x samples tables, so memory and time barely grow with the
number of cells; they grow with samples x clusters x markers (LMM and GLMM fit one model per
cluster, or per cluster and marker). Image: 591 MB compressed, 2.3 GB unpacked (R + ~150
packages, mostly diffcyt's imports).

## Results on the demo study

Run on flowsom's output for the CytoBatchFlagR demo (24 samples: 6 batches x tcu1, tcu2
unstimulated, tcs3, tcs4 stimulated 7 days with anti-CD3/CD28 + IL-2), `REFERENCE_LEVEL=unstimulated`,
batch as a fixed effect, edgeR / limma, FDR 0.05:

- **DA: 13 of 20 metaclusters** change (8 up, 5 down). Up in stimulated: NK CD56dim CD16+
  (mc 15, 1.8 % -> 19 % of cells, logFC 3.5), the CD25+FOXP3+ CD4 T cluster (mc 1, Treg or
  activated, 1.2 % -> 5.1 %, logFC 2.3), CD45RA- CD4 and CD8 T clusters (mc 3, 5, 8, 2), CD8 T
  CD45RA+ CD56+ (mc 11), NK CD56bright (mc 13). Down: the two DN T (gd-like) clusters (mc 6,
  10), NK CD56dim CD16int (mc 14), CD16+ HLA-DR+ monocyte-like (mc 20) and B cells (mc 19).
  voom calls 16 and GLMM 14, the same strongest hits.
- **DS: 142 of 200 cluster-marker pairs** change (limma; LMM: 146). The strongest: CD25 in the
  CD25+FOXP3+ CD4 T cluster (+1.6 arcsinh units), HLA-DR in NK cells, TIGIT in DN and CD8 T
  cells, GZMB up in nearly every T and NK cluster (18 of 19 GZMB hits up), HLA-DR (14/15 up),
  CD69 (13/15), CD25 (12/15) and TIGIT (12/12) up; PFN mostly down.
- Batch 06's unstimulated samples stand out in `proportions.png` (e.g. ~14 % DN T and ~12 %
  CD16+ HLA-DR+ cells where the other batches have 1-3 % and ~0 %); the mc 20 "down" call rests
  on those two samples alone, so read it as a batch artefact rather than biology. That is the
  dataset's purpose (CytoBatchFlagR); the harmonisation workflow addresses it.

## Limitations

- DA is about proportions within each sample: a large expansion of one population shifts every
  other population's proportion a little (compositional effect). The tables give logFC so such
  small shifts can be read as what they are.
- GLMM gives no logFC of its own and LMM no coefficient; their logFC columns are descriptive
  group differences (named in `summary.json`).
- Only "each level vs the reference" contrasts; no interaction terms or continuous covariates.
- A block is either a fixed effect (edgeR/voom/limma) or a random intercept (GLMM/LMM), as the
  method allows; voom's and limma's `duplicateCorrelation` route is not used.
- DS needs flowsom's per-sample medians; markers not in `medians_by_sample.csv` cannot be tested.

## Licences

Our code: Apache-2.0 (`LICENSE`). The image bundles R and GPL-licensed R packages (edgeR,
limma, lme4, multcomp, ...) that our script loads in-process, so the image is distributed under
the GPL; diffcyt itself is MIT. See `THIRD_PARTY.md`.

## Citations

Checked against PubMed (E-utilities `esummary`) and, for lme4, Crossref; the demo data via
DataCite.

- Weber LM, Nowicka M, Soneson C, Robinson MD. diffcyt: Differential discovery in
  high-dimensional cytometry via high-resolution clustering. *Commun Biol* 2:183 (2019).
  doi:10.1038/s42003-019-0415-5. PMID 31098416.
- Nowicka M, Krieg C, Crowell HL, Weber LM, Hartmann FJ, Guglietta S, Becher B, Levesque MP,
  Robinson MD. CyTOF workflow: differential discovery in high-throughput high-dimensional
  cytometry datasets. *F1000Research* 6:748 (2017). doi:10.12688/f1000research.11622.3.
  PMID 28663787.
- Van Gassen S, Callebaut B, Van Helden MJ, Lambrecht BN, Demeester P, Dhaene T, Saeys Y.
  FlowSOM: Using self-organizing maps for visualization and interpretation of cytometry data.
  *Cytometry A* 87(7):636-645 (2015). doi:10.1002/cyto.a.22625. PMID 25573116. (the clusters)
- Robinson MD, McCarthy DJ, Smyth GK. edgeR: a Bioconductor package for differential
  expression analysis of digital gene expression data. *Bioinformatics* 26(1):139-140 (2010).
  doi:10.1093/bioinformatics/btp616. PMID 19910308.
- McCarthy DJ, Chen Y, Smyth GK. Differential expression analysis of multifactor RNA-Seq
  experiments with respect to biological variation. *Nucleic Acids Res* 40(10):4288-4297
  (2012). doi:10.1093/nar/gks042. PMID 22287627. (edgeR GLMs)
- Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK. limma powers differential
  expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Res* 43(7):e47
  (2015). doi:10.1093/nar/gkv007. PMID 25605792.
- Law CW, Chen Y, Shi W, Smyth GK. voom: Precision weights unlock linear model analysis tools
  for RNA-seq read counts. *Genome Biol* 15(2):R29 (2014). doi:10.1186/gb-2014-15-2-r29.
  PMID 24485249.
- Bates D, Mächler M, Bolker B, Walker S. Fitting Linear Mixed-Effects Models Using lme4.
  *J Stat Softw* 67(1) (2015). doi:10.18637/jss.v067.i01. (not in PubMed; checked on Crossref)
- Hothorn T, Bretz F, Westfall P. Simultaneous inference in general parametric models.
  *Biom J* 50(3):346-363 (2008). doi:10.1002/bimj.200810425. PMID 18481363. (multcomp `glht`)

Demo data: Koenig Z, Andorf S, Tilburgs T, Eswar S. CytoBatchFlagR: Spectral flow cytometry
dataset. Zenodo (2025). doi:10.5281/zenodo.15388817. CC BY 4.0.
