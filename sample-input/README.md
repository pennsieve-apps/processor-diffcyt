# sample-input

This app reads the *outputs* of two upstream steps of the Immune Profiling & Differential
Abundance workflow, so there is no public file to download for it directly:

| Input folder | Files | Produced by |
|---|---|---|
| `clusters/` | `counts.csv`, `medians_by_sample.csv`, `labels.csv` | processor-flowsom (`clusters` output) |
| `metadata/` | `samples.csv`, `panel.csv` | processor-fcs-prep (`metadata` output) |

## Synthetic fixture

```sh
make fixture      # writes data/fixture/linked/ and data/fixture/flat/ (gitignored)
make run          # runs the Docker image over data/fixture/linked/ into data/output/
```

`tests/make_fixture.R` writes a study shaped like the demo data (6 batches x {tcu1, tcu2
unstimulated, tcs3, tcs4 stimulated}, 12 metaclusters, 20 markers, 15-60k cells per sample) with
planted effects: Treg 4x and B cells 1/3 in stimulated samples (DA), CD25 +1.2 in CD4 memory and
HLA-DR +1.0 in CD8 memory (DS, arcsinh units), and a rare cluster that the `MIN_CELLS` filter
drops. Everything else carries only batch shifts and per-sample noise. `linked/` has the two
input folders; `flat/` has the same files as symlinks side by side in nested folders (how
`INPUT_DIR` looks when a workflow uses no links).

## Real data

Run the upstream steps on the demo study (CytoBatchFlagR spectral flow dataset, Koenig et al.,
Zenodo 2025, doi:10.5281/zenodo.15388817, CC BY 4.0): fcs-prep -> peacoqc -> cleanet -> flowsom.
Then put flowsom's `clusters` output and fcs-prep's `metadata` output side by side and run:

```sh
mkdir -p data/real
cp -r <flowsom-output>/clusters data/real/clusters
cp -r <prep-output>/metadata data/real/metadata
make run INPUT=./data/real
```

The demo study's reference level is `unstimulated` (`dev.env` sets `REFERENCE_LEVEL`); without it
the reference would be the alphabetically first level, `stimulated`, and every sign flips.
