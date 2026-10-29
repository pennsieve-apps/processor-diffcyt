# Inputs, parameters and small output helpers for processor-diffcyt.

# Files are found by name: first in INPUT_DIR/<port>/ (placement: subfolder),
# else anywhere under INPUT_DIR (no links used). INPUT_DIR entries are symlinks,
# so the walk follows them.
find_input <- function(input_dir, port, name) {
  linked <- file.path(input_dir, port)
  if (dir.exists(linked)) {
    hits <- list_files_follow(linked, name)
    if (length(hits)) return(hits[[1]])
  }
  hits <- list_files_follow(input_dir, name)
  if (length(hits) > 1) {
    log_line(sprintf("warning: %d files named %s under INPUT_DIR, using %s", length(hits), name, hits[[1]]))
  }
  if (length(hits)) hits[[1]] else NULL
}

# list.files(recursive = TRUE) does not descend into symlinked directories, so
# walk the tree ourselves (dir.exists/file.exists follow links). Sorted, so the
# shallowest, alphabetically first match wins.
list_files_follow <- function(root, name, depth = 0L, max_depth = 12L) {
  if (depth > max_depth || !dir.exists(root)) return(character())
  entries <- sort(list.files(root, all.files = FALSE, full.names = TRUE, no.. = TRUE))
  here <- entries[basename(entries) == name & !dir.exists(entries)]
  sub <- entries[dir.exists(entries)]
  deeper <- unlist(lapply(sub, list_files_follow, name = name, depth = depth + 1L, max_depth = max_depth))
  c(here, deeper)
}

read_table <- function(path) {
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("", "NA"),
           colClasses = "character")
}

require_columns <- function(df, cols, what) {
  missing <- setdiff(cols, colnames(df))
  if (length(missing)) stop(sprintf("%s lacks column(s): %s", what, paste(missing, collapse = ", ")))
}

env_str <- function(name, default = "") {
  v <- Sys.getenv(name, unset = NA)
  if (is.na(v)) default else trimws(v)
}

# A parameter that the platform may pass as "" when the user clears it.
env_num <- function(name, default, integer = FALSE) {
  v <- env_str(name, "")
  if (v == "") return(default)
  x <- suppressWarnings(as.numeric(v))
  if (is.na(x)) stop(sprintf("%s must be a number, got '%s'", name, v))
  if (integer) x <- as.integer(round(x))
  x
}

read_params <- function() {
  block <- Sys.getenv("BLOCK_COLUMN", unset = NA)
  block <- if (is.na(block)) "batch" else trimws(block)
  if (tolower(block) %in% c("none", "null", "false")) block <- ""
  da <- env_str("DA_METHOD", "edgeR")
  ds <- env_str("DS_METHOD", "limma")
  da <- c(edger = "edgeR", voom = "voom", glmm = "GLMM")[tolower(da)]
  ds <- c(limma = "limma", lmm = "LMM")[tolower(ds)]
  if (is.na(da)) stop(sprintf("DA_METHOD must be edgeR, voom or GLMM, got '%s'", env_str("DA_METHOD")))
  if (is.na(ds)) stop(sprintf("DS_METHOD must be limma or LMM, got '%s'", env_str("DS_METHOD")))
  fdr <- env_num("FDR", 0.05)
  if (fdr <= 0 || fdr >= 1) stop("FDR must be between 0 and 1")
  state <- env_str("STATE_MARKERS", "")
  list(
    CONDITION_COLUMN = env_str("CONDITION_COLUMN", "condition"),
    REFERENCE_LEVEL = env_str("REFERENCE_LEVEL", ""),
    BLOCK_COLUMN = block,
    DA_METHOD = unname(da),
    DS_METHOD = unname(ds),
    MIN_CELLS = env_num("MIN_CELLS", 3L, integer = TRUE),
    MIN_SAMPLES = env_num("MIN_SAMPLES", 0, integer = FALSE),
    FDR = fdr,
    STATE_MARKERS = if (state == "") character() else trimws(strsplit(state, "[,;[:space:]]+")[[1]])
  )
}

log_line <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", ..., "\n", sep = " ")
  flush(stdout())
}

# Logical columns as true/false (like panel.csv), numbers to 6 significant digits.
write_csv <- function(df, path) {
  for (k in names(df)) {
    if (is.logical(df[[k]])) df[[k]] <- ifelse(df[[k]], "true", "false")
    else if (is.double(df[[k]])) df[[k]] <- signif(df[[k]], 6)
  }
  write.csv(df, path, row.names = FALSE, na = "")
}

write_json <- function(x, path) {
  writeLines(jsonlite::toJSON(x, auto_unbox = TRUE, pretty = TRUE, na = "null", null = "null", digits = NA), path)
}

rel_path <- function(path, root) {
  path <- normalizePath(path, mustWork = FALSE)
  root <- normalizePath(root, mustWork = FALSE)
  if (startsWith(path, paste0(root, "/"))) substring(path, nchar(root) + 2L) else path
}

# Metacluster ids sort numerically when they are all numbers.
order_ids <- function(ids) {
  ids <- unique(as.character(ids))
  num <- suppressWarnings(as.numeric(ids))
  if (!anyNA(num)) ids[order(num)] else sort(ids)
}

safe_name <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)
