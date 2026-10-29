#!/usr/bin/env Rscript
# Synthetic stand-in for the flowsom `clusters` output and the prep `metadata`
# output, shaped like the demo study (6 batches x {tcu1, tcu2 unstimulated,
# tcs3, tcs4 stimulated}), with effects planted where the tests look for them:
#
#   DA  metacluster 3 "Treg"          4x more frequent in stimulated samples
#       metacluster 8 "B cells"       3x less frequent in stimulated samples
#   DS  metacluster 2 "CD4 memory"    CD25 +1.2 (arcsinh) in stimulated samples
#       metacluster 5 "CD8 memory"    HLADR +1.0 in stimulated samples
#   metacluster 12 "rare"             0-2 cells per sample: filtered by MIN_CELLS
#
# Everything else is null: per-batch shifts (abundance and expression) and
# per-sample noise only. Deterministic (set.seed).
#
# Usage: Rscript tests/make_fixture.R <out_dir> [linked|flat|both] [three]
#   linked: <out>/linked/clusters/*.csv, <out>/linked/metadata/*.csv
#   flat:   <out>/flat/ with the same files as symlinks, side by side in nested
#           folders (how INPUT_DIR looks when no links are used)
#   three:  add a third condition level ("costim", batches 01-06 x one sample,
#           with metacluster 3 2x) to exercise several comparisons

make_fixture <- function(out_dir, layout = "both", three_levels = FALSE, seed = 1L) {
  set.seed(seed)
  labels <- c("CD4 naive", "CD4 memory", "Treg", "CD8 naive", "CD8 memory",
              "NK CD56bright", "NK CD56dim", "B cells", "gd T", "MAIT", "DN T", "rare")
  base_frac <- c(0.20, 0.18, 0.04, 0.12, 0.10, 0.02, 0.10, 0.09, 0.05, 0.04, 0.06, 0.00008)
  base_frac <- base_frac / sum(base_frac)
  n_mc <- length(labels)

  samples <- expand.grid(ctl = c("tcu1", "tcu2", "tcs3", "tcs4"), b = 1:6, stringsAsFactors = FALSE)
  if (three_levels) {
    samples <- rbind(samples, data.frame(ctl = "tcc5", b = 1:6))
  }
  samples$batch <- sprintf("batch%02d", samples$b)
  samples$sample_id <- paste0(samples$batch, "_", samples$ctl)
  samples$condition <- ifelse(grepl("^tcu", samples$ctl), "unstimulated",
                              ifelse(grepl("^tcs", samples$ctl), "stimulated", "costim"))
  samples$file <- paste0(sprintf("B%02d_", samples$b), samples$ctl, "_CD45+CD14-.fcs")
  samples$control <- samples$ctl

  markers <- c("CD3", "CD4", "CD8", "CD19", "CD56", "CD16", "CD45RA", "CCR7", "CD27", "CD95",
               "CD25", "CD69", "HLADR", "PD1", "CD39", "GITR", "TIGIT", "GZMB", "PFN", "FOXP3")
  state <- c("CD25", "CD69", "HLADR", "PD1", "CD39", "GITR", "TIGIT", "GZMB", "PFN")
  clustering <- c("CD3", "CD4", "CD8", "CD19", "CD56", "CD16", "CD45RA", "CCR7", "CD27", "CD95", "FOXP3")
  base_med <- matrix(runif(n_mc * length(markers), 0.2, 3.5), n_mc, dimnames = list(NULL, markers))
  batch_abund <- matrix(rnorm(6 * n_mc, 0, 0.15), 6)
  batch_expr <- matrix(rnorm(6 * length(markers), 0, 0.10), 6)

  counts <- list()
  medians <- list()
  for (i in seq_len(nrow(samples))) {
    s <- samples[i, ]
    total <- round(runif(1, 15000, 60000))
    w <- base_frac * exp(batch_abund[s$b, ] + rnorm(n_mc, 0, 0.12))
    if (s$condition == "stimulated") w[3] <- w[3] * 4
    if (s$condition == "stimulated") w[8] <- w[8] / 3
    if (s$condition == "costim") w[3] <- w[3] * 2
    n <- as.vector(rmultinom(1, total, w / sum(w)))
    n[12] <- sample(0:2, 1)
    counts[[i]] <- data.frame(sample_id = s$sample_id, metacluster = seq_len(n_mc), label = labels,
                              n = n, frac = n / sum(n))
    for (k in seq_len(n_mc)) {
      if (n[k] == 0) next
      noise_sd <- 0.06 + 0.5 / sqrt(n[k])
      m <- base_med[k, ] + batch_expr[s$b, ] + rnorm(length(markers), 0, noise_sd)
      if (s$condition == "stimulated" && k == 2) m["CD25"] <- m["CD25"] + 1.2
      if (s$condition == "stimulated" && k == 5) m["HLADR"] <- m["HLADR"] + 1.0
      medians[[length(medians) + 1]] <- data.frame(sample_id = s$sample_id, metacluster = k,
                                                   marker = markers, median = round(m, 5))
    }
  }
  counts <- do.call(rbind, counts)
  medians <- do.call(rbind, medians)
  n_events <- tapply(counts$n, counts$sample_id, sum)
  samples_csv <- data.frame(sample_id = samples$sample_id, file = samples$file, batch = samples$batch,
                            condition = samples$condition, n_events = as.integer(n_events[samples$sample_id]),
                            control = samples$control)
  panel <- rbind(
    data.frame(channel = c("FSC-A", "SSC-A"), marker = c("FSC-A", "SSC-A"), type = "scatter",
               cofactor = "", clustering = "false", state = "false"),
    data.frame(channel = paste0("CH", seq_along(markers), "-A"), marker = markers, type = "marker",
               cofactor = "6000", clustering = ifelse(markers %in% clustering, "true", "false"),
               state = ifelse(markers %in% state, "true", "false")),
    data.frame(channel = "Time", marker = "Time", type = "time", cofactor = "", clustering = "false",
               state = "false"))

  paths <- list()
  linked <- file.path(out_dir, "linked")
  dir.create(file.path(linked, "clusters"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(linked, "metadata"), recursive = TRUE, showWarnings = FALSE)
  wcsv <- function(x, p) write.csv(x, p, row.names = FALSE, quote = TRUE)
  wcsv(counts, file.path(linked, "clusters", "counts.csv"))
  wcsv(medians, file.path(linked, "clusters", "medians_by_sample.csv"))
  wcsv(data.frame(metacluster = seq_len(n_mc), label = labels), file.path(linked, "clusters", "labels.csv"))
  wcsv(samples_csv, file.path(linked, "metadata", "samples.csv"))
  wcsv(panel, file.path(linked, "metadata", "panel.csv"))
  paths$linked <- normalizePath(linked)

  if (layout %in% c("flat", "both")) {
    flat <- file.path(out_dir, "flat")
    dir.create(file.path(flat, "run-1", "tables"), recursive = TRUE, showWarnings = FALSE)
    # Relative links, so the fixture still resolves when mounted elsewhere.
    for (f in c("counts.csv", "medians_by_sample.csv", "labels.csv")) {
      file.symlink(file.path("..", "..", "..", "linked", "clusters", f), file.path(flat, "run-1", "tables", f))
    }
    for (f in c("samples.csv", "panel.csv")) {
      file.symlink(file.path("..", "linked", "metadata", f), file.path(flat, f))
    }
    paths$flat <- normalizePath(flat)
  }
  invisible(paths)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 1) stop("usage: make_fixture.R <out_dir> [linked|flat|both] [three]")
  p <- make_fixture(args[1], if (length(args) >= 2) args[2] else "both", three_levels = "three" %in% args)
  cat(paste(names(p), unlist(p), sep = ": "), sep = "\n")
}
