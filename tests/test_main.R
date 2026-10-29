#!/usr/bin/env Rscript
# End-to-end tests: the synthetic fixture (tests/make_fixture.R) through
# R/main.R, run as the platform runs it (Rscript with INPUT_DIR/OUTPUT_DIR and
# parameters in the environment). Plain base R, no test framework.
#
#   Rscript tests/test_main.R          (from the repo root; needs the R packages
#                                       of the image: make test runs it in Docker
#                                       when Rscript is not installed)

args <- commandArgs(trailingOnly = FALSE)
f <- sub("^--file=", "", args[grepl("^--file=", args)])
REPO <- normalizePath(file.path(dirname(normalizePath(f[1])), ".."))
source(file.path(REPO, "tests", "make_fixture.R"))
MAIN <- file.path(REPO, "R", "main.R")
PARAMS <- c("CONDITION_COLUMN", "REFERENCE_LEVEL", "BLOCK_COLUMN", "DA_METHOD", "DS_METHOD",
            "MIN_CELLS", "MIN_SAMPLES", "FDR", "STATE_MARKERS")

tmp <- tempfile("diffcyt-test-")
dir.create(tmp)
fx <- make_fixture(file.path(tmp, "fx"), "both")
fx3 <- make_fixture(file.path(tmp, "fx3"), "linked", three_levels = TRUE)

run_main <- function(input, name, ...) {
  out <- file.path(tmp, name)
  env <- c(INPUT_DIR = input, OUTPUT_DIR = out, unlist(list(...)))
  # Clear parameters this run doesn't set, so the defaults are exercised.
  unset <- setdiff(PARAMS, names(env))
  old <- Sys.getenv(unset, unset = NA, names = TRUE)
  Sys.unsetenv(unset)
  on.exit(for (k in names(old)) if (!is.na(old[[k]])) do.call(Sys.setenv, setNames(list(old[[k]]), k)))
  envs <- paste0(names(env), "=", shQuote(env))
  log <- suppressWarnings(system2("env", c(envs, file.path(R.home("bin"), "Rscript"), MAIN), stdout = TRUE, stderr = TRUE))
  code <- attr(log, "status")
  if (is.null(code)) code <- 0L
  list(code = code, out = out, log = log)
}

results <- list()
check <- function(name, expr) {
  ok <- tryCatch({ force(expr); TRUE }, error = function(e) { message(sprintf("FAIL %s: %s", name, conditionMessage(e))); FALSE })
  if (ok) message(sprintf("ok   %s", name))
  results[[name]] <<- ok
}
expect <- function(cond, msg) if (!isTRUE(cond)) stop(msg, call. = FALSE)
read_out <- function(r, port, file) {
  d <- read.csv(file.path(r$out, port, file), stringsAsFactors = FALSE, check.names = FALSE)
  for (k in intersect(c("tested", "significant"), colnames(d))) {
    expect(all(d[[k]] %in% c("true", "false")), paste(k, "is written as true/false"))
    d[[k]] <- d[[k]] == "true"
  }
  d
}
read_json <- function(r, port, file) jsonlite::fromJSON(file.path(r$out, port, file), simplifyVector = FALSE)
show_log <- function(r) message(paste(r$log, collapse = "\n"))

DA_PLANTED <- c("3", "8")
DS_PLANTED <- list(c("2", "CD25"), c("5", "HLADR"))

check_planted <- function(r, da_method = "edgeR", ds_method = "limma", ds = TRUE) {
  if (r$code != 0) show_log(r)
  expect(r$code == 0, sprintf("exit code %s", r$code))
  da <- read_out(r, "da", "da.csv")
  for (col in c("metacluster", "label", "comparison", "logFC", "logCPM", "p_val", "p_adj",
                "mean_frac_stimulated", "mean_frac_unstimulated"))
    expect(col %in% colnames(da), paste("da.csv lacks", col))
  sig <- da$metacluster[da$significant]
  expect(all(DA_PLANTED %in% sig), sprintf("%s: planted DA clusters not found (significant: %s)", da_method, paste(sig, collapse = ",")))
  # Planting 4x / 0.33x shifts every other cluster's proportion by a few per
  # cent (compositional effect); a test with high power (GLMM) may call those,
  # so extra hits are allowed only with a small effect.
  extra <- da[da$significant & !da$metacluster %in% DA_PLANTED, ]
  expect(all(abs(extra$logFC) < 0.35), sprintf("%s: false DA hits with large effect: %s", da_method,
                                               paste(sprintf("%s (%.2f)", extra$metacluster, extra$logFC), collapse = ", ")))
  expect(da$logFC[da$metacluster == "3"] > 1, "Treg logFC should be clearly positive (4x planted)")
  expect(da$logFC[da$metacluster == "8"] < -1, "B cell logFC should be clearly negative (1/3 planted)")
  expect(!da$tested[da$metacluster == "12"] && is.na(da$p_val[da$metacluster == "12"]), "rare cluster should be filtered")
  expect(da$label[da$metacluster == "3"] == "Treg", "labels carried")
  for (f in c("volcano.png", "proportions.png", "summary.json", "provenance.json"))
    expect(file.size(file.path(r$out, "da", f)) > 0, paste("missing da/", f))
  if (ds) {
    d <- read_out(r, "ds", "ds.csv")
    for (col in c("metacluster", "label", "marker", "logFC", "p_val", "p_adj")) expect(col %in% colnames(d), paste("ds.csv lacks", col))
    hits <- paste(d$metacluster, d$marker)[d$significant]
    for (p in DS_PLANTED) expect(paste(p, collapse = " ") %in% hits,
                                 sprintf("%s: planted DS %s not found (hits: %s)", ds_method, paste(p, collapse = "/"), paste(hits, collapse = ", ")))
    expect(length(setdiff(hits, sapply(DS_PLANTED, paste, collapse = " "))) <= 2,
           sprintf("%s: too many false DS hits: %s", ds_method, paste(hits, collapse = ", ")))
    expect(setequal(unique(d$marker), c("CD25", "CD69", "HLADR", "PD1", "CD39", "GITR", "TIGIT", "GZMB", "PFN")),
           "state markers default to panel state == true")
    for (f in c("heatmap.png", "logfc.png", "summary.json", "provenance.json"))
      expect(file.size(file.path(r$out, "ds", f)) > 0, paste("missing ds/", f))
  }
  invisible(TRUE)
}

check("linked layout, defaults (edgeR + limma, block=batch fixed)", {
  r <- run_main(fx$linked, "linked", REFERENCE_LEVEL = "unstimulated")
  check_planted(r)
  prov <- read_json(r, "da", "provenance.json")
  expect(prov$app == "processor-diffcyt" && prov$version == "0.1.0", "provenance app/version")
  expect(prov$parameters$BLOCK_USED == "batch" && prov$parameters$MIN_SAMPLES_USED == 6, "block and min_samples resolved")
  expect(length(prov$samples) == 24 && all(sapply(prov$samples, `[[`, "status") == "included"), "24 samples included")
  expect(all(c("clusters/counts.csv", "metadata/samples.csv") %in% unlist(prov$inputs)), "inputs recorded relative")
  s <- read_json(r, "da", "summary.json")
  expect(s$comparisons[[1]]$comparison == "stimulated vs unstimulated", "comparison name")
  expect(any(grepl("^.* sample batch01_tcu1: condition=unstimulated batch=batch01 cells=[0-9]+ clusters=1[12] included", r$log)), "per-sample log line")
})

check("flat layout (no links, symlinked files, nested)", {
  r <- run_main(fx$flat, "flat", REFERENCE_LEVEL = "unstimulated")
  check_planted(r)
})

check("voom + LMM (block random for LMM)", {
  r <- run_main(fx$linked, "voom-lmm", REFERENCE_LEVEL = "unstimulated", DA_METHOD = "voom", DS_METHOD = "LMM")
  check_planted(r, "voom", "LMM")
  expect(read_json(r, "ds", "summary.json")$block_effect == "random", "LMM block is random")
})

check("GLMM (block random), no DS markers given -> panel", {
  r <- run_main(fx$linked, "glmm", REFERENCE_LEVEL = "unstimulated", DA_METHOD = "GLMM", DS_METHOD = "limma")
  check_planted(r, "GLMM")
  expect(grepl("descriptive", read_json(r, "da", "summary.json")$logFC), "GLMM logFC is descriptive")
})

check("no block (BLOCK_COLUMN empty) and explicit STATE_MARKERS", {
  r <- run_main(fx$linked, "noblock", REFERENCE_LEVEL = "unstimulated", BLOCK_COLUMN = "",
                STATE_MARKERS = "cd25, HLADR,CD69,NOPE")
  if (r$code != 0) show_log(r)
  expect(r$code == 0, "exit code")
  d <- read_out(r, "ds", "ds.csv")
  expect(setequal(unique(d$marker), c("CD25", "HLADR", "CD69")), "STATE_MARKERS matched ignoring case, unknown dropped")
  s <- read_json(r, "ds", "summary.json")
  expect(any(grepl("NOPE", unlist(s$notes))), "unknown marker reported")
  expect(read_json(r, "da", "provenance.json")$parameters$BLOCK_USED == "", "no block")
  da <- read_out(r, "da", "da.csv")
  expect(all(DA_PLANTED %in% da$metacluster[da$significant]), "planted DA without block")
})

check("default reference is the first level alphabetically", {
  r <- run_main(fx$linked, "defref")
  expect(r$code == 0, "exit code")
  da <- read_out(r, "da", "da.csv")
  expect(unique(da$comparison) == "unstimulated vs stimulated", "reference = stimulated")
  expect(da$logFC[da$metacluster == "3"] < -1, "sign follows the reference")
})

check("three condition levels: each vs reference", {
  r <- run_main(fx3$linked, "three", REFERENCE_LEVEL = "unstimulated")
  if (r$code != 0) show_log(r)
  expect(r$code == 0, "exit code")
  da <- read_out(r, "da", "da.csv")
  expect(setequal(unique(da$comparison), c("costim vs unstimulated", "stimulated vs unstimulated")), "two comparisons")
  expect("mean_frac_costim" %in% colnames(da), "mean frac per level")
  st <- da[da$comparison == "stimulated vs unstimulated", ]
  expect(all(DA_PLANTED %in% st$metacluster[st$significant]), "planted DA in stimulated vs unstimulated")
  co <- da[da$comparison == "costim vs unstimulated", ]
  expect(co$significant[co$metacluster == "3"] && co$logFC[co$metacluster == "3"] > 0.5, "costim Treg 2x found")
  d <- read_out(r, "ds", "ds.csv")
  expect(nrow(d) == 2 * 12 * 9, "ds rows per comparison")
})

check("samples without a condition are excluded, not fatal", {
  src <- file.path(tmp, "fx-na")
  file.copy(fx$linked, tmp, recursive = TRUE)
  file.rename(file.path(tmp, "linked"), src)
  s <- read.csv(file.path(src, "metadata", "samples.csv"), stringsAsFactors = FALSE)
  s$condition[s$sample_id == "batch02_tcs3"] <- NA
  s <- s[s$sample_id != "batch03_tcu1", ]
  write.csv(s, file.path(src, "metadata", "samples.csv"), row.names = FALSE, na = "")
  r <- run_main(src, "na", REFERENCE_LEVEL = "unstimulated")
  expect(r$code == 0, "exit code")
  prov <- read_json(r, "da", "provenance.json")
  st <- setNames(sapply(prov$samples, `[[`, "status"), sapply(prov$samples, `[[`, "sample_id"))
  expect(st[["batch02_tcs3"]] == "excluded" && st[["batch03_tcu1"]] == "excluded" && sum(st == "included") == 22,
         "two excluded")
})

check("DS without medians fails alone; DA still succeeds", {
  src <- file.path(tmp, "fx-nomed")
  dir.create(file.path(src, "clusters"), recursive = TRUE)
  file.copy(file.path(fx$linked, "clusters", "counts.csv"), file.path(src, "clusters"))
  file.copy(file.path(fx$linked, "metadata"), src, recursive = TRUE)
  r <- run_main(src, "nomed", REFERENCE_LEVEL = "unstimulated")
  expect(r$code == 0, "exit code 0 when DA succeeded")
  expect(read_json(r, "ds", "summary.json")$status == "failed", "ds failed")
  expect(file.exists(file.path(r$out, "da", "da.csv")), "da written")
})

check("missing inputs and bad parameters fail the step", {
  empty <- file.path(tmp, "empty")
  dir.create(empty)
  expect(run_main(empty, "empty")$code != 0, "no input -> non-zero")
  expect(run_main(fx$linked, "badref", REFERENCE_LEVEL = "nope")$code != 0, "unknown reference -> non-zero")
  expect(run_main(fx$linked, "badmethod", DA_METHOD = "DESeq2")$code != 0, "unknown method -> non-zero")
  expect(run_main(fx$linked, "badblock", BLOCK_COLUMN = "donor")$code != 0, "unknown block column -> non-zero")
})

n_fail <- sum(!unlist(results))
message(sprintf("\n%d passed, %d failed", sum(unlist(results)), n_fail))
unlink(tmp, recursive = TRUE)
quit(save = "no", status = if (n_fail) 1L else 0L)
