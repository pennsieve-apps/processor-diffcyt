#!/usr/bin/env Rscript
# processor-diffcyt: differential abundance (DA) and differential state (DS)
# tests between conditions on FlowSOM metaclusters with diffcyt
# (Weber et al. 2019), following the Nowicka et al. (2017) CyTOF workflow.
#
#   INPUT_DIR/clusters/  counts.csv, medians_by_sample.csv, labels.csv  (flowsom)
#   INPUT_DIR/metadata/  samples.csv, panel.csv                          (fcs-prep)
#   OUTPUT_DIR/da/       da.csv, volcano.png, proportions.png, summary.json, provenance.json
#   OUTPUT_DIR/ds/       ds.csv, heatmap.png, logfc.png, summary.json, provenance.json
#
# The step fails (non-zero exit) when inputs are missing or unusable, or when
# both tests failed; a failed DS (e.g. no medians) with a good DA is logged and
# recorded in ds/summary.json.

local({
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args[grepl("^--file=", args)])
  here <- if (length(f)) dirname(normalizePath(f[1])) else file.path(Sys.getenv("APP_ROOT", "."), "R")
  assign("R_DIR", here, envir = globalenv())
})
source(file.path(R_DIR, "io.R"))

APP_NAME <- "processor-diffcyt"
TOOLS <- c("diffcyt", "edgeR", "limma", "lme4", "multcomp", "SummarizedExperiment", "jsonlite")

app_version <- function() {
  p <- file.path(dirname(R_DIR), "VERSION")
  if (file.exists(p)) trimws(readLines(p, n = 1L, warn = FALSE)) else "unknown"
}

top_hits <- function(tab, cols, n = 10L) {
  s <- tab[tab$significant, , drop = FALSE]
  s <- head(s[order(s$p_adj), cols, drop = FALSE], n)
  if (!nrow(s)) return(list())
  lapply(seq_len(nrow(s)), function(i) as.list(s[i, , drop = FALSE]))
}

run <- function() {
  t_start <- Sys.time()
  timings <- list()
  tick <- function(name, t0) timings[[name]] <<- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 2)

  params <- read_params()
  input_dir <- Sys.getenv("INPUT_DIR")
  output_dir <- Sys.getenv("OUTPUT_DIR")
  if (input_dir == "" || output_dir == "") stop("INPUT_DIR and OUTPUT_DIR must be set")
  if (!dir.exists(input_dir)) stop(sprintf("INPUT_DIR %s does not exist", input_dir))

  files <- list(
    counts = find_input(input_dir, "clusters", "counts.csv"),
    medians = find_input(input_dir, "clusters", "medians_by_sample.csv"),
    labels = find_input(input_dir, "clusters", "labels.csv"),
    samples = find_input(input_dir, "metadata", "samples.csv"),
    panel = find_input(input_dir, "metadata", "panel.csv")
  )
  for (k in names(files)) log_line(sprintf("input %s: %s", k, if (is.null(files[[k]])) "(not found)" else rel_path(files[[k]], input_dir)))
  if (is.null(files$counts)) stop("counts.csv (flowsom clusters output) not found under INPUT_DIR")
  if (is.null(files$samples)) stop("samples.csv (prep metadata output) not found under INPUT_DIR")

  t0 <- Sys.time()
  log_line("loading R packages")
  source(file.path(R_DIR, "analysis.R"))
  source(file.path(R_DIR, "plots.R"))
  tick("load_packages_s", t0)

  t0 <- Sys.time()
  study <- load_study(files, params)
  panel <- if (!is.null(files$panel)) read_table(files$panel) else NULL
  tick("load_inputs_s", t0)
  for (i in seq_len(nrow(study$status))) {
    s <- study$status[i, ]
    log_line(sprintf("sample %s: %s=%s%s cells=%d clusters=%d %s%s", s$sample_id, params$CONDITION_COLUMN,
                     s$condition, if (params$BLOCK_COLUMN != "") sprintf(" %s=%s", params$BLOCK_COLUMN, s$block) else "",
                     as.integer(s$n_cells), sum(study$raw_counts$n[study$raw_counts$sample_id == s$sample_id] > 0),
                     s$status, if (s$reason != "") paste0(" (", s$reason, ")") else ""))
  }
  block_used <- "block" %in% colnames(study$info)
  if (!is.null(study$block_note)) log_line(sprintf("warning: BLOCK_COLUMN '%s' not used: %s", params$BLOCK_COLUMN, study$block_note))
  n_per <- table(study$info$condition)
  log_line(sprintf("%d samples, %d metaclusters; %s levels %s (reference %s); block: %s; min_samples=%s",
                   nrow(study$info), length(study$clusters), params$CONDITION_COLUMN,
                   paste(sprintf("%s (%d)", names(n_per), as.integer(n_per)), collapse = ", "), study$reference,
                   if (block_used) params$BLOCK_COLUMN else "none", format(study$min_samples)))

  resolved <- c(params, list(
    REFERENCE_LEVEL_USED = study$reference,
    BLOCK_USED = if (block_used) params$BLOCK_COLUMN else "",
    BLOCK_EFFECT = if (!block_used) "none" else c(edgeR = "fixed", voom = "fixed", GLMM = "random")[[params$DA_METHOD]],
    BLOCK_EFFECT_DS = if (!block_used) "none" else c(limma = "fixed", LMM = "random")[[params$DS_METHOD]],
    MIN_SAMPLES_USED = study$min_samples
  ))
  resolved$STATE_MARKERS <- as.list(params$STATE_MARKERS)
  tools <- c(list(R = paste(R.version$major, R.version$minor, sep = ".")),
             lapply(setNames(TOOLS, TOOLS), function(p) as.character(utils::packageVersion(p))))
  inputs_seen <- lapply(Filter(Negate(is.null), files), rel_path, root = input_dir)
  samples_status <- lapply(seq_len(nrow(study$status)), function(i) as.list(study$status[i, ]))
  group_sizes <- as.list(setNames(as.integer(n_per), names(n_per)))

  provenance <- function(out, extra) {
    c(list(app = APP_NAME, version = app_version(), output = out,
           generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
           tools = tools, parameters = resolved, inputs = inputs_seen,
           conditions = list(levels = study$levels, reference = study$reference, samples_per_level = group_sizes),
           block_note = study$block_note, timings = timings, samples = samples_status), extra)
  }

  ok <- c(da = FALSE, ds = FALSE)

  # ---- differential abundance ----
  da_dir <- file.path(output_dir, "da")
  dir.create(da_dir, recursive = TRUE, showWarnings = FALSE)
  t0 <- Sys.time()
  da_err <- NULL
  da <- tryCatch(run_da(study, params), error = function(e) { da_err <<- conditionMessage(e); NULL })
  tick("da_test_s", t0)
  if (!is.null(da)) {
    t0 <- Sys.time()
    write_csv(da$table, file.path(da_dir, "da.csv"))
    plot_volcano(da$table, params$FDR, file.path(da_dir, "volcano.png"))
    plot_proportions(study, da$table, da$frac, params$FDR, file.path(da_dir, "proportions.png"))
    tick("da_outputs_s", t0)
    cmp_summary <- lapply(unique(da$table$comparison), function(cm) {
      d <- da$table[da$table$comparison == cm, ]
      log_line(sprintf("DA %s (%s): %d/%d clusters tested, %d significant (%d up, %d down) at FDR %g",
                       cm, params$DA_METHOD, sum(d$tested), nrow(d), sum(d$significant),
                       sum(d$significant & d$logFC > 0, na.rm = TRUE), sum(d$significant & d$logFC < 0, na.rm = TRUE), params$FDR))
      list(comparison = cm, clusters = nrow(d), tested = sum(d$tested), significant = sum(d$significant),
           up = sum(d$significant & d$logFC > 0, na.rm = TRUE), down = sum(d$significant & d$logFC < 0, na.rm = TRUE),
           top = top_hits(d, c("metacluster", "label", "logFC", "p_adj")))
    })
    write_json(list(status = "ok", method = params$DA_METHOD, fdr = params$FDR, reference = study$reference,
                    block = resolved$BLOCK_USED, block_effect = resolved$BLOCK_EFFECT,
                    min_cells = params$MIN_CELLS, min_samples = study$min_samples, samples_per_level = group_sizes,
                    logFC = da$logfc_source, logCPM = da$logcpm_source, comparisons = cmp_summary),
               file.path(da_dir, "summary.json"))
    ok[["da"]] <- TRUE
  } else {
    log_line(sprintf("DA failed: %s", da_err))
    write_json(list(status = "failed", method = params$DA_METHOD, error = da_err), file.path(da_dir, "summary.json"))
  }
  write_json(provenance("da", list(status = if (ok[["da"]]) "ok" else "failed", error = da_err)),
             file.path(da_dir, "provenance.json"))

  # ---- differential state ----
  ds_dir <- file.path(output_dir, "ds")
  dir.create(ds_dir, recursive = TRUE, showWarnings = FALSE)
  t0 <- Sys.time()
  ds_err <- NULL
  ds <- tryCatch({
    if (is.null(files$medians)) stop("medians_by_sample.csv (flowsom clusters output) not found under INPUT_DIR")
    run_ds(study, params, read_table(files$medians), panel)
  }, error = function(e) { ds_err <<- conditionMessage(e); NULL })
  tick("ds_test_s", t0)
  if (!is.null(ds)) {
    for (n in ds$state$notes) log_line(paste("warning:", n))
    t0 <- Sys.time()
    write_csv(ds$table, file.path(ds_dir, "ds.csv"))
    plot_ds_heatmap(study, ds$table, ds$d_medians, params$FDR, file.path(ds_dir, "heatmap.png"))
    plot_ds_logfc(study, ds$table, file.path(ds_dir, "logfc.png"))
    tick("ds_outputs_s", t0)
    cmp_summary <- lapply(unique(ds$table$comparison), function(cm) {
      d <- ds$table[ds$table$comparison == cm, ]
      log_line(sprintf("DS %s (%s): %d/%d cluster-marker pairs tested, %d significant at FDR %g",
                       cm, params$DS_METHOD, sum(d$tested), nrow(d), sum(d$significant), params$FDR))
      list(comparison = cm, pairs = nrow(d), tested = sum(d$tested), significant = sum(d$significant),
           up = sum(d$significant & d$logFC > 0, na.rm = TRUE), down = sum(d$significant & d$logFC < 0, na.rm = TRUE),
           markers_with_hits = unique(d$marker[d$significant]),
           top = top_hits(d, c("metacluster", "label", "marker", "logFC", "p_adj")))
    })
    write_json(list(status = "ok", method = params$DS_METHOD, fdr = params$FDR, reference = study$reference,
                    block = resolved$BLOCK_USED, block_effect = resolved$BLOCK_EFFECT_DS,
                    state_markers = ds$state$markers, state_markers_source = ds$state$source,
                    notes = ds$state$notes, logFC = ds$logfc_source, comparisons = cmp_summary),
               file.path(ds_dir, "summary.json"))
    ok[["ds"]] <- TRUE
  } else {
    log_line(sprintf("DS failed: %s", ds_err))
    write_json(list(status = "failed", method = params$DS_METHOD, error = ds_err), file.path(ds_dir, "summary.json"))
  }
  timings$total_s <- round(as.numeric(difftime(Sys.time(), t_start, units = "secs")), 2)
  write_json(provenance("ds", list(status = if (ok[["ds"]]) "ok" else "failed", error = ds_err,
                                   state_markers = if (!is.null(ds)) ds$state$markers)),
             file.path(ds_dir, "provenance.json"))
  log_line(sprintf("done in %.1f s: da %s, ds %s", timings$total_s, if (ok[["da"]]) "ok" else "FAILED",
                   if (ok[["ds"]]) "ok" else "FAILED"))
  if (any(ok)) 0L else 1L
}

if (sys.nframe() == 0L) {
  code <- tryCatch(run(), error = function(e) {
    log_line(sprintf("error: %s", conditionMessage(e)))
    1L
  })
  quit(save = "no", status = code)
}
