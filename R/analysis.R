# The study (samples, clusters, counts, medians) and the diffcyt tests.
#
# diffcyt's testDA_* / testDS_* functions take two SummarizedExperiments that
# its own calcCounts() / calcMedians() would compute from per-cell data. Here
# they are built directly from flowsom's per-sample tables, in exactly the
# shape those functions read:
#
#   d_counts   assay "counts": clusters x samples (n cells); rownames = cluster
#              ids; rowData$cluster_id a factor whose levels are the rows in
#              order (the tests index result rows by as.numeric(cluster_id));
#              colData$sample_id.
#   d_medians  one assay per marker, each clusters x samples (median arcsinh
#              expression, NA where the sample has no cells in the cluster),
#              same rows and columns as d_counts; metadata(d_medians) has
#              id_type_markers / id_state_markers, logical over
#              names(assays(d_medians)) (testDS_* default markers_to_test).

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(S4Vectors)
  library(diffcyt)
})

load_study <- function(files, params) {
  counts <- read_table(files$counts)
  require_columns(counts, c("sample_id", "metacluster", "n"), "counts.csv")
  counts$n <- as.numeric(counts$n)
  if (anyNA(counts$n) || any(counts$n < 0)) stop("counts.csv: n must be a non-negative number in every row")
  counts$metacluster <- as.character(counts$metacluster)

  samples <- read_table(files$samples)
  require_columns(samples, c("sample_id", params$CONDITION_COLUMN), "samples.csv")
  if (anyDuplicated(samples$sample_id)) stop("samples.csv: sample_id is not unique")
  use_block <- params$BLOCK_COLUMN != ""
  if (use_block && !params$BLOCK_COLUMN %in% colnames(samples)) {
    stop(sprintf("BLOCK_COLUMN '%s' is not a column of samples.csv (columns: %s); set it empty for no block",
                 params$BLOCK_COLUMN, paste(colnames(samples), collapse = ", ")))
  }

  # Cluster labels: counts.csv's label column, else labels.csv, else the id.
  clusters <- order_ids(counts$metacluster)
  label <- setNames(rep(NA_character_, length(clusters)), clusters)
  if ("label" %in% colnames(counts)) {
    l <- counts[!is.na(counts$label), c("metacluster", "label")]
    label[l$metacluster] <- l$label
  }
  if (!is.null(files$labels) && anyNA(label)) {
    lab <- read_table(files$labels)
    if (all(c("metacluster", "label") %in% colnames(lab))) {
      lab <- lab[lab$metacluster %in% names(label) & is.na(label[lab$metacluster]), ]
      label[lab$metacluster] <- lab$label
    }
  }
  label[is.na(label) | label == ""] <- names(label)[is.na(label) | label == ""]

  # Per-sample status.
  ids <- union(samples$sample_id, unique(counts$sample_id))
  status <- data.frame(sample_id = ids, condition = NA_character_, block = NA_character_,
                       n_cells = 0, status = "included", reason = "", stringsAsFactors = FALSE)
  rownames(status) <- ids
  n_by <- tapply(counts$n, counts$sample_id, sum)
  status[names(n_by), "n_cells"] <- n_by
  m <- match(ids, samples$sample_id)
  status$condition <- samples[[params$CONDITION_COLUMN]][m]
  if (use_block) status$block <- samples[[params$BLOCK_COLUMN]][m]
  exclude <- function(which, reason) {
    which <- which & status$status == "included"
    status$status[which] <<- "excluded"
    status$reason[which] <<- reason
  }
  exclude(is.na(m), "not in samples.csv")
  exclude(!ids %in% counts$sample_id, "not in counts.csv")
  exclude(status$n_cells <= 0, "no cells")
  exclude(is.na(status$condition) | toupper(status$condition) == "NA",
          sprintf("no %s", params$CONDITION_COLUMN))
  if (use_block) exclude(is.na(status$block), sprintf("no %s", params$BLOCK_COLUMN))
  inc <- status$status == "included"
  sizes <- table(status$condition[inc])
  small <- names(sizes)[sizes < 2]
  exclude(status$condition %in% small, sprintf("%s level has fewer than 2 samples", params$CONDITION_COLUMN))

  inc <- status$status == "included"
  levels_all <- sort(unique(status$condition[inc]))
  if (length(levels_all) < 2) {
    stop(sprintf("need at least two %s levels with 2+ samples each; found: %s", params$CONDITION_COLUMN,
                 paste(sprintf("%s (%d)", names(sizes), as.integer(sizes)), collapse = ", ")))
  }
  ref <- params$REFERENCE_LEVEL
  if (ref == "") {
    ref <- levels_all[1]
  } else if (!ref %in% levels_all) {
    stop(sprintf("REFERENCE_LEVEL '%s' is not a %s level with 2+ samples (levels: %s)", ref,
                 params$CONDITION_COLUMN, paste(levels_all, collapse = ", ")))
  }
  lv <- c(ref, setdiff(levels_all, ref))

  keep <- status[inc, ]
  keep <- keep[order(factor(keep$condition, levels = lv), keep$block, keep$sample_id), ]
  info <- data.frame(sample_id = factor(keep$sample_id, levels = keep$sample_id),
                     condition = factor(keep$condition, levels = lv))
  block_note <- NULL
  if (use_block) {
    info$block <- factor(keep$block, levels = sort(unique(keep$block)))
    block_note <- check_block(info)
    if (!is.null(block_note)) info$block <- NULL
  }

  mat <- matrix(0, length(clusters), nrow(info), dimnames = list(clusters, as.character(info$sample_id)))
  cc <- counts[counts$sample_id %in% colnames(mat), ]
  cc <- aggregate(n ~ metacluster + sample_id, data = cc, FUN = sum)
  mat[cbind(cc$metacluster, cc$sample_id)] <- cc$n

  min_samples <- params$MIN_SAMPLES
  min_samples_auto <- min_samples <= 0
  if (min_samples_auto) min_samples <- min(table(info$condition)) / 2

  list(info = info, status = status, clusters = clusters, label = label, levels = lv, reference = ref,
       counts = mat, block_note = block_note, min_samples = min_samples, min_samples_auto = min_samples_auto,
       raw_counts = counts)
}

# A block (batch) is usable when it has 2+ levels and adding it leaves a
# full-rank design with residual degrees of freedom; otherwise it is dropped
# (e.g. a batch per sample, or batches confounded with the condition).
check_block <- function(info) {
  if (nlevels(droplevels(info$block)) < 2) return("block has fewer than 2 levels")
  X <- model.matrix(~ condition + block, data = info)
  if (qr(X)$rank < ncol(X)) return("block is confounded with the condition (design not full rank)")
  if (nrow(X) - ncol(X) < 1) return("block leaves no residual degrees of freedom")
  NULL
}

build_d_counts <- function(study) {
  SummarizedExperiment(
    assays = list(counts = study$counts),
    rowData = DataFrame(cluster_id = factor(study$clusters, levels = study$clusters)),
    colData = DataFrame(study$info)
  )
}

build_d_medians <- function(study, medians, panel, state_markers) {
  require_columns(medians, c("sample_id", "metacluster", "marker", "median"), "medians_by_sample.csv")
  medians$median <- as.numeric(medians$median)
  medians$metacluster <- as.character(medians$metacluster)
  markers <- unique(medians$marker)
  if (!is.null(panel)) {
    ord <- match(markers, panel$marker)
    markers <- markers[order(is.na(ord), ord)]
  }
  smp <- as.character(study$info$sample_id)
  med <- medians[medians$sample_id %in% smp & medians$metacluster %in% study$clusters, ]
  assays <- lapply(setNames(markers, markers), function(mk) {
    a <- matrix(NA_real_, length(study$clusters), length(smp), dimnames = list(study$clusters, smp))
    x <- med[med$marker == mk, ]
    a[cbind(x$metacluster, x$sample_id)] <- x$median
    a
  })
  is_type <- if (!is.null(panel) && "clustering" %in% colnames(panel)) {
    markers %in% panel$marker[tolower(panel$clustering) == "true"]
  } else rep(FALSE, length(markers))
  d <- SummarizedExperiment(
    assays = assays,
    rowData = DataFrame(cluster_id = factor(study$clusters, levels = study$clusters)),
    colData = DataFrame(study$info),
    metadata = list(id_type_markers = is_type, id_state_markers = markers %in% state_markers)
  )
  d
}

# STATE_MARKERS, else panel state == true. Names are matched exactly, then
# ignoring case; unknown names are reported and dropped.
resolve_state_markers <- function(requested, panel, available) {
  notes <- character()
  if (length(requested)) {
    exact <- requested[requested %in% available]
    rest <- setdiff(requested, exact)
    ci <- available[match(toupper(rest), toupper(available))]
    if (any(is.na(ci))) notes <- c(notes, sprintf("STATE_MARKERS not in medians_by_sample.csv: %s",
                                                  paste(rest[is.na(ci)], collapse = ", ")))
    chosen <- unique(c(exact, ci[!is.na(ci)]))
    source <- "STATE_MARKERS"
  } else {
    if (is.null(panel) || !"state" %in% colnames(panel)) {
      stop("no STATE_MARKERS given and no panel.csv with a state column to take them from")
    }
    wanted <- panel$marker[tolower(panel$state) == "true"]
    missing <- setdiff(wanted, available)
    if (length(missing)) notes <- c(notes, sprintf("state markers without medians (skipped): %s",
                                                   paste(missing, collapse = ", ")))
    chosen <- intersect(wanted, available)
    source <- "panel.csv state == true"
  }
  if (!length(chosen)) stop("no state markers to test (", source, ")")
  list(markers = available[available %in% chosen], source = source, notes = notes)
}

# Comparisons: each non-reference level against the reference, one contrast
# each, all from the same model.
comparisons <- function(study) {
  lapply(study$levels[-1], function(l) list(level = l, reference = study$reference,
                                             name = sprintf("%s vs %s", l, study$reference)))
}

fixed_design <- function(info) {
  if ("block" %in% colnames(info)) model.matrix(~ condition + block, data = info)
  else model.matrix(~ condition, data = info)
}

contrast_for <- function(X, level) {
  ct <- matrix(0, ncol(X), 1, dimnames = list(colnames(X), NULL))
  ct[paste0("condition", level), 1] <- 1
  ct
}

# Rows kept by diffcyt's filter: >= MIN_CELLS cells in >= min_samples samples.
kept_clusters <- function(study, min_cells) {
  rowSums(study$counts >= min_cells) >= study$min_samples
}

group_means <- function(mat, info, levels, prefix) {
  out <- sapply(levels, function(l) rowMeans(mat[, info$condition == l, drop = FALSE], na.rm = TRUE))
  out <- matrix(out, nrow = nrow(mat), dimnames = list(rownames(mat), paste0(prefix, safe_name(levels))))
  out[is.nan(out)] <- NA
  as.data.frame(out, check.names = FALSE)
}

run_da <- function(study, params) {
  d_counts <- build_d_counts(study)
  info <- study$info
  frac <- sweep(study$counts, 2, colSums(study$counts), "/")
  means <- group_means(frac, info, study$levels, "mean_frac_")
  kept <- kept_clusters(study, params$MIN_CELLS)
  if (sum(kept) < 2) stop(sprintf("only %d cluster(s) pass MIN_CELLS=%d in MIN_SAMPLES=%s samples",
                                  sum(kept), params$MIN_CELLS, format(study$min_samples)))
  avelogcpm <- edgeR::aveLogCPM(edgeR::DGEList(study$counts))
  out <- list()
  for (cmp in comparisons(study)) {
    res <- switch(params$DA_METHOD,
      edgeR = {
        X <- fixed_design(info)
        testDA_edgeR(d_counts, X, contrast_for(X, cmp$level),
                     min_cells = params$MIN_CELLS, min_samples = study$min_samples)
      },
      voom = {
        X <- fixed_design(info)
        testDA_voom(d_counts, X, contrast_for(X, cmp$level),
                    min_cells = params$MIN_CELLS, min_samples = study$min_samples)
      },
      GLMM = {
        # Observation-level random effect (sample_id) for overdispersion, as in
        # the diffcyt vignette, plus the block as a random intercept.
        f <- createFormula(as.data.frame(info), cols_fixed = "condition",
                           cols_random = c("sample_id", if ("block" %in% colnames(info)) "block"))
        X <- model.matrix(~ condition, data = info)
        testDA_GLMM(d_counts, f, t(contrast_for(X, cmp$level)),
                    min_cells = params$MIN_CELLS, min_samples = study$min_samples)
      })
    rd <- as.data.frame(rowData(res))
    stopifnot(identical(as.character(rd$cluster_id), study$clusters))
    ref_col <- paste0("mean_frac_", safe_name(cmp$reference))
    lvl_col <- paste0("mean_frac_", safe_name(cmp$level))
    descriptive <- log2((means[[lvl_col]] + 1e-6) / (means[[ref_col]] + 1e-6))
    logfc <- if (!is.null(rd$logFC)) rd$logFC else ifelse(kept, descriptive, NA)
    logcpm <- if (!is.null(rd$logCPM)) rd$logCPM else avelogcpm
    df <- data.frame(metacluster = study$clusters, label = unname(study$label[study$clusters]),
                     comparison = cmp$name, logFC = logfc, logCPM = logcpm,
                     p_val = rd$p_val, p_adj = rd$p_adj, tested = kept,
                     stringsAsFactors = FALSE)
    df$significant <- !is.na(df$p_adj) & df$p_adj < params$FDR
    out[[cmp$name]] <- cbind(df, means)
  }
  da <- do.call(rbind, out)
  rownames(da) <- NULL
  list(table = da, frac = frac,
       logfc_source = if (params$DA_METHOD == "GLMM") "descriptive: log2(mean_frac level / mean_frac reference)"
                      else sprintf("%s model coefficient (log2 fold change)", params$DA_METHOD),
       logcpm_source = if (params$DA_METHOD == "edgeR") "edgeR" else "edgeR::aveLogCPM over all samples")
}

run_ds <- function(study, params, medians, panel) {
  available <- unique(medians$marker)
  sm <- resolve_state_markers(params$STATE_MARKERS, panel, available)
  d_counts <- build_d_counts(study)
  d_medians <- build_d_medians(study, medians, panel, sm$markers)
  info <- study$info
  kept <- kept_clusters(study, params$MIN_CELLS)
  if (sum(kept) < 1) stop(sprintf("no cluster passes MIN_CELLS=%d in MIN_SAMPLES=%s samples",
                                  params$MIN_CELLS, format(study$min_samples)))
  state_names <- names(assays(d_medians))[metadata(d_medians)$id_state_markers]
  out <- list()
  for (cmp in comparisons(study)) {
    res <- switch(params$DS_METHOD,
      limma = {
        X <- fixed_design(info)
        testDS_limma(d_counts, d_medians, X, contrast_for(X, cmp$level),
                     min_cells = params$MIN_CELLS, min_samples = study$min_samples)
      },
      LMM = {
        f <- createFormula(as.data.frame(info), cols_fixed = "condition",
                           cols_random = if ("block" %in% colnames(info)) "block" else NULL)
        X <- model.matrix(~ condition, data = info)
        testDS_LMM(d_counts, d_medians, f, t(contrast_for(X, cmp$level)),
                   min_cells = params$MIN_CELLS, min_samples = study$min_samples)
      })
    rd <- as.data.frame(rowData(res))
    rd$cluster_id <- as.character(rd$cluster_id)
    rd$marker_id <- as.character(rd$marker_id)
    # Group means of the per-sample medians, for each (cluster, marker) row.
    gm <- do.call(rbind, lapply(state_names, function(mk) {
      g <- group_means(assay(d_medians, mk), info, study$levels, "mean_median_")
      g
    }))
    ref_col <- paste0("mean_median_", safe_name(cmp$reference))
    lvl_col <- paste0("mean_median_", safe_name(cmp$level))
    tested <- rep(kept, length(state_names))
    logfc <- if (!is.null(rd$logFC)) rd$logFC else ifelse(tested, gm[[lvl_col]] - gm[[ref_col]], NA)
    df <- data.frame(metacluster = rd$cluster_id, label = unname(study$label[rd$cluster_id]),
                     marker = rd$marker_id, comparison = cmp$name, logFC = logfc,
                     p_val = rd$p_val, p_adj = rd$p_adj, tested = tested, stringsAsFactors = FALSE)
    df$significant <- !is.na(df$p_adj) & df$p_adj < params$FDR
    out[[cmp$name]] <- cbind(df, gm)
  }
  ds <- do.call(rbind, out)
  rownames(ds) <- NULL
  list(table = ds, d_medians = d_medians, state = sm,
       logfc_source = if (params$DS_METHOD == "limma") "limma model coefficient (difference in median arcsinh expression)"
                      else "descriptive: mean_median level - mean_median reference (arcsinh)")
}
