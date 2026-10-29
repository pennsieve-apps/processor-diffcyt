# Third-party software

This app's own code (`R/`, `processor/`, `tests/`) is Apache-2.0 (see `LICENSE`). The image
it builds bundles the following; their licences apply to the image. Our R script loads
diffcyt and its dependencies in-process, so the running program combines Apache-2.0 code
with GPL-licensed R packages (and R itself); the image as a whole is therefore distributed
under the terms of the GPL (version 2 or later, as the components allow). Every component
is distributed unmodified, with its licence text in its installed package folder
(`/usr/local/lib/R/site-library/<pkg>/`, `/usr/local/lib/R/library/<pkg>/`).

| Component | Version | What | Licence | Source |
|---|---|---|---|---|
| diffcyt | 1.32.1 | DA / DS tests (`testDA_edgeR`, `testDA_voom`, `testDA_GLMM`, `testDS_limma`, `testDS_LMM`, `createFormula`) | MIT | https://bioconductor.org/packages/diffcyt |
| edgeR | 4.10.5 | negative binomial GLMs (DA, edgeR) and `aveLogCPM` | GPL (>= 2) | https://bioconductor.org/packages/edgeR |
| limma | 3.68.5 | voom, moderated linear models (DA voom, DS limma) | GPL (>= 2) | https://bioconductor.org/packages/limma |
| lme4 | 2.0-6 | GLMM / LMM fits | GPL (>= 2) | https://cran.r-project.org/package=lme4 |
| multcomp | 1.4-32 | `glht` contrast tests for GLMM / LMM | GPL-2 | https://cran.r-project.org/package=multcomp |
| SummarizedExperiment, S4Vectors | 1.42.0, 0.50.3 | diffcyt's input/output containers | Artistic-2.0 | Bioconductor |
| jsonlite | 2.0.0 | JSON output | MIT | https://cran.r-project.org/package=jsonlite |
| Other diffcyt imports (installed, loaded as namespaces) | Bioconductor 3.23 / CRAN snapshot 2026-10-01 | FlowSOM (GPL >= 2), flowCore (Artistic-2.0), cytolib (**AGPL-3.0-only**, a flowCore dependency), ComplexHeatmap (MIT), circlize (MIT), dplyr / tidyr / ggplot2 (MIT), igraph (GPL >= 2), Matrix (GPL >= 2), ... 149 packages in all; list them with `installed.packages(fields = "License")` in the image | as listed in each package's DESCRIPTION | Posit Public Package Manager binaries |
| R | 4.6.1 | interpreter, base graphics (PNG via Cairo) | GPL-2 \| GPL-3 | https://www.r-project.org |
| rocker/r-ver | 4.6.1 (digest pinned in `Dockerfile`) | base image (Ubuntu 24.04 + R) | GPL-2 (Rocker scripts); Ubuntu packages under their own free licences | https://github.com/rocker-org/rocker-versioned2 |
| awslambdaric | 3.1.1 | AWS Lambda runtime interface client (Lambda runs only) | Apache-2.0 | https://github.com/aws/aws-lambda-python-runtime-interface-client |
| Python 3.12 (Ubuntu) | 3.12 | runs the entry point / Lambda handler | PSF | Ubuntu 24.04 |

cytolib is AGPL-3.0: the image ships it unmodified, as installed from Bioconductor, and this
app does not call it directly (it is loaded because diffcyt imports flowCore).

## Demo data

The tests use only synthetic data. The demo study the workflow runs on is the CytoBatchFlagR
spectral flow cytometry dataset: Koenig Z, Andorf S, Tilburgs T, Eswar S. Zenodo (2025),
doi:10.5281/zenodo.15388817, CC BY 4.0. Not redistributed here.
