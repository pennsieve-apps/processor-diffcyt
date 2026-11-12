"""
Synthetic stand-in for the flowsom `clusters` output and the fcs-prep `metadata` output,
for `review.generate`: the platform's generator image has Python, not R, so this is
tests/make_fixture.R in Python (the R tests keep using the R version). Same study shape
(6 batches x {tcu1, tcu2 unstimulated, tcs3, tcs4 stimulated}) and the same planted effects:

  DA  metacluster 3 "Treg"          4x more frequent in stimulated samples
      metacluster 8 "B cells"       3x less frequent in stimulated samples
  DS  metacluster 2 "CD4 memory"    CD25 +1.2 (arcsinh) in stimulated samples
      metacluster 5 "CD8 memory"    HLADR +1.0 in stimulated samples
  metacluster 12 "rare"             0-2 cells per sample: filtered by MIN_CELLS

Everything else is null: per-batch shifts (abundance and expression) and per-sample noise
only. Deterministic (seeded).

    python -m tests.make_fixture OUT [--layout linked|flat|both] [--three]

  linked: OUT/linked/clusters/*.csv, OUT/linked/metadata/*.csv
  flat:   OUT/flat/ with the same files as relative symlinks, side by side in nested folders
  three:  a third condition level ("costim", batches 01-06 x one sample, metacluster 3 2x)
"""
from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path

import numpy as np

LABELS = ["CD4 naive", "CD4 memory", "Treg", "CD8 naive", "CD8 memory", "NK CD56bright", "NK CD56dim",
          "B cells", "gd T", "MAIT", "DN T", "rare"]
BASE_FRAC = np.array([0.20, 0.18, 0.04, 0.12, 0.10, 0.02, 0.10, 0.09, 0.05, 0.04, 0.06, 0.00008])
MARKERS = ["CD3", "CD4", "CD8", "CD19", "CD56", "CD16", "CD45RA", "CCR7", "CD27", "CD95",
           "CD25", "CD69", "HLADR", "PD1", "CD39", "GITR", "TIGIT", "GZMB", "PFN", "FOXP3"]
STATE = {"CD25", "CD69", "HLADR", "PD1", "CD39", "GITR", "TIGIT", "GZMB", "PFN"}
CLUSTERING = {"CD3", "CD4", "CD8", "CD19", "CD56", "CD16", "CD45RA", "CCR7", "CD27", "CD95", "FOXP3"}


def write_csv(path: Path, header: list[str], rows) -> None:
    with open(path, "w", newline="") as f:
        w = csv.writer(f, quoting=csv.QUOTE_NONNUMERIC)
        w.writerow(header)
        w.writerows(rows)


def make_fixture(out_dir: Path, layout: str = "both", three_levels: bool = False, seed: int = 1) -> dict:
    rng = np.random.default_rng(seed)
    frac = BASE_FRAC / BASE_FRAC.sum()
    n_mc = len(LABELS)
    samples = [(ctl, b) for b in range(1, 7) for ctl in ("tcu1", "tcu2", "tcs3", "tcs4")]
    if three_levels:
        samples += [("tcc5", b) for b in range(1, 7)]
    base_med = rng.uniform(0.2, 3.5, (n_mc, len(MARKERS)))
    batch_abund = rng.normal(0, 0.15, (6, n_mc))
    batch_expr = rng.normal(0, 0.10, (6, len(MARKERS)))

    counts, medians, sample_rows = [], [], []
    for ctl, b in samples:
        batch = f"batch{b:02d}"
        sample_id = f"{batch}_{ctl}"
        condition = "unstimulated" if ctl.startswith("tcu") else "stimulated" if ctl.startswith("tcs") else "costim"
        total = int(round(rng.uniform(15000, 60000)))
        w = frac * np.exp(batch_abund[b - 1] + rng.normal(0, 0.12, n_mc))
        if condition == "stimulated":
            w[2] *= 4
            w[7] /= 3
        if condition == "costim":
            w[2] *= 2
        n = rng.multinomial(total, w / w.sum())
        n[11] = rng.integers(0, 3)
        tot = int(n.sum())
        for k in range(n_mc):
            counts.append([sample_id, k + 1, LABELS[k], int(n[k]), n[k] / tot])
            if n[k] == 0:
                continue
            noise_sd = 0.06 + 0.5 / np.sqrt(n[k])
            m = base_med[k] + batch_expr[b - 1] + rng.normal(0, noise_sd, len(MARKERS))
            if condition == "stimulated" and k == 1:
                m[MARKERS.index("CD25")] += 1.2
            if condition == "stimulated" and k == 4:
                m[MARKERS.index("HLADR")] += 1.0
            medians.extend([sample_id, k + 1, mk, round(float(v), 5)] for mk, v in zip(MARKERS, m))
        sample_rows.append([sample_id, f"B{b:02d}_{ctl}_CD45+CD14-.fcs", batch, condition, tot, ctl])
    panel = [["FSC-A", "FSC-A", "scatter", "", "false", "false"], ["SSC-A", "SSC-A", "scatter", "", "false", "false"]]
    panel += [[f"CH{i + 1}-A", mk, "marker", "6000", str(mk in CLUSTERING).lower(), str(mk in STATE).lower()]
              for i, mk in enumerate(MARKERS)]
    panel.append(["Time", "Time", "time", "", "false", "false"])

    paths = {}
    linked = out_dir / "linked"
    (linked / "clusters").mkdir(parents=True, exist_ok=True)
    (linked / "metadata").mkdir(parents=True, exist_ok=True)
    write_csv(linked / "clusters" / "counts.csv", ["sample_id", "metacluster", "label", "n", "frac"], counts)
    write_csv(linked / "clusters" / "medians_by_sample.csv", ["sample_id", "metacluster", "marker", "median"], medians)
    write_csv(linked / "clusters" / "labels.csv", ["metacluster", "label"], [[k + 1, l] for k, l in enumerate(LABELS)])
    write_csv(linked / "metadata" / "samples.csv", ["sample_id", "file", "batch", "condition", "n_events", "control"],
              sample_rows)
    write_csv(linked / "metadata" / "panel.csv", ["channel", "marker", "type", "cofactor", "clustering", "state"], panel)
    paths["linked"] = str(linked.resolve())

    if layout in ("flat", "both"):
        flat = out_dir / "flat"
        (flat / "run-1" / "tables").mkdir(parents=True, exist_ok=True)
        # Relative links, so the fixture still resolves when mounted elsewhere.
        for f in ("counts.csv", "medians_by_sample.csv", "labels.csv"):
            os.symlink(Path("..", "..", "..", "linked", "clusters", f), flat / "run-1" / "tables" / f)
        for f in ("samples.csv", "panel.csv"):
            os.symlink(Path("..", "linked", "metadata", f), flat / f)
        paths["flat"] = str(flat.resolve())
    return paths


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out", type=Path)
    ap.add_argument("--layout", choices=("linked", "flat", "both"), default="both")
    ap.add_argument("--three", action="store_true", help="add a third condition level")
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    for name, p in make_fixture(a.out, a.layout, a.three, a.seed).items():
        print(f"{name}: {p}")
