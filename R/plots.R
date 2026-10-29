# Figures (base R graphics, PNG). Colours: categorical slots in fixed order for
# condition levels; a blue <-> gray <-> red diverging scale for signed effects.

PAL <- list(
  groups = c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"),
  up = "#d03b3b", down = "#2a78d6", ns = "#b5b4ae", ink = "#0b0b0b", ink2 = "#52514e",
  grid = "#e6e5e0", surface = "#fcfcfb", mid = "#f0efec"
)

group_colors <- function(levels) {
  cols <- PAL$groups[(seq_along(levels) - 1L) %% length(PAL$groups) + 1L]
  setNames(cols, levels)
}

diverging <- function(n = 101) {
  colorRampPalette(c("#104281", PAL$down, "#9ec5f4", PAL$mid, "#f4b0ae", PAL$up, "#7a1a1a"))(n)
}

open_png <- function(path, width, height) {
  png(path, width = width, height = height, units = "in", res = 120, type = "cairo", bg = PAL$surface)
  par(fg = PAL$ink, col.axis = PAL$ink2, col.lab = PAL$ink2, col.main = PAL$ink, family = "sans",
      las = 1, mgp = c(2.2, 0.6, 0), tcl = -0.25)
}

message_png <- function(path, text) {
  open_png(path, 7, 2.5)
  par(mar = c(0, 0, 0, 0))
  plot.new()
  text(0.5, 0.5, text, col = PAL$ink2, cex = 1.1)
  dev.off()
}

# Deterministic spread of points around a centre (no random jitter).
spread <- function(n, width = 0.28) {
  if (n <= 1) return(0)
  seq(-width, width, length.out = n)[order(rep_len(c(1, -1), n) * seq_len(n))]
}

plot_volcano <- function(da, fdr, path) {
  cmps <- unique(da$comparison)
  open_png(path, 5.2 * length(cmps), 4.6)
  par(mfrow = c(1, length(cmps)), mar = c(4, 4, 3, 1))
  for (cm in cmps) {
    d <- da[da$comparison == cm & da$tested & !is.na(da$p_val), ]
    y <- -log10(pmax(d$p_val, 1e-300))
    lim <- max(1, abs(d$logFC), na.rm = TRUE) * 1.15
    plot(NA, xlim = c(-lim, lim), ylim = c(0, max(1, y) * 1.12), xlab = "log2 fold change",
         ylab = expression(-log[10]~p), main = cm, font.main = 1, cex.main = 1, bty = "n")
    abline(v = 0, col = PAL$grid)
    col <- ifelse(!d$significant, PAL$ns, ifelse(d$logFC > 0, PAL$up, PAL$down))
    points(d$logFC, y, pch = 21, bg = col, col = PAL$surface, cex = 1.6, lwd = 1.2)
    s <- d$significant
    if (any(s)) text(d$logFC[s], y[s], d$label[s], pos = ifelse(d$logFC[s] > 0, 4, 2), cex = 0.75,
                     col = PAL$ink, xpd = NA)
    legend("topleft", legend = c(sprintf("up (FDR < %g)", fdr), "down", "not significant"),
           pt.bg = c(PAL$up, PAL$down, PAL$ns), pch = 21, col = PAL$surface, bty = "n", cex = 0.75,
           text.col = PAL$ink2)
  }
  dev.off()
}

plot_proportions <- function(study, da, frac, fdr, path) {
  info <- study$info
  cols <- group_colors(study$levels)
  best <- tapply(ifelse(is.na(da$p_adj), 2, da$p_adj), da$metacluster, min)
  ids <- names(sort(best[study$clusters]))
  ncol <- min(4L, length(ids))
  nrow <- ceiling(length(ids) / ncol)
  pch_block <- if ("block" %in% colnames(info)) (as.integer(info$block) - 1L) %% 5L + 21L else rep(21L, nrow(info))
  open_png(path, 3.1 * ncol, 2.9 * nrow + 0.5)
  layout(rbind(matrix(seq_len(nrow * ncol), nrow, byrow = TRUE), nrow * ncol + 1L), heights = c(rep(2.9, nrow), 0.5))
  par(mar = c(2.2, 4.4, 3.2, 0.6), mgp = c(3, 0.6, 0))
  for (id in ids) {
    v <- 100 * frac[id, ]
    plot(NA, xlim = c(0.4, length(study$levels) + 0.6), ylim = c(0, max(v, 1e-3) * 1.08), xaxt = "n",
         xlab = "", ylab = "% of cells", bty = "n", main = "")
    axis(1, at = seq_along(study$levels), labels = study$levels, cex.axis = 0.75, lwd = 0, line = -0.6)
    rows <- da[da$metacluster == id, ]
    sub <- paste(ifelse(is.na(rows$p_adj), "not tested",
                        sprintf("%s p_adj %s", ifelse(rows$significant, "*", ""), format(signif(rows$p_adj, 2)))),
                 collapse = "; ")
    title(main = sprintf("%s (%s)", study$label[[id]], id), cex.main = 0.95, font.main = 2, line = 1.6)
    mtext(sub, side = 3, line = 0.5, cex = 0.62, col = PAL$ink2)
    for (k in seq_along(study$levels)) {
      sel <- info$condition == study$levels[k]
      x <- v[sel]
      q <- quantile(x, c(0.25, 0.5, 0.75), na.rm = TRUE)
      rect(k - 0.32, q[1], k + 0.32, q[3], border = PAL$ink2, col = NA, lwd = 1)
      segments(k - 0.32, q[2], k + 0.32, q[2], col = PAL$ink, lwd = 2)
      points(k + spread(length(x)), x, pch = pch_block[sel], bg = cols[k], col = PAL$surface, cex = 1.3, lwd = 1)
    }
  }
  for (i in seq_len(nrow * ncol - length(ids))) plot.new()
  par(mar = c(0, 0, 0, 0))
  plot.new()
  leg <- names(cols)
  legend("center", legend = leg, pt.bg = cols, pch = 21, col = PAL$surface, horiz = TRUE, bty = "n",
         cex = 0.9, text.col = PAL$ink2, pt.cex = 1.5,
         title = if ("block" %in% colnames(info)) "condition (point shape = block)" else "condition")
  dev.off()
}

# Significant (cluster, marker) pairs x samples: per-sample medians, z-scored
# per row; samples ordered by condition then block.
plot_ds_heatmap <- function(study, ds, d_medians, fdr, path, max_rows = 60L) {
  sig <- ds[ds$significant, ]
  if (!nrow(sig)) {
    message_png(path, sprintf("No cluster-marker pair is significant at FDR < %g.", fdr))
    return(invisible(0L))
  }
  sig <- sig[order(sig$p_adj), ]
  sig <- sig[!duplicated(paste(sig$metacluster, sig$marker)), ]
  n_sig <- nrow(sig)
  sig <- head(sig, max_rows)
  info <- study$info
  mat <- t(sapply(seq_len(nrow(sig)), function(i) assay(d_medians, sig$marker[i])[sig$metacluster[i], ]))
  if (nrow(sig) == 1) mat <- matrix(mat, nrow = 1)
  z <- t(apply(mat, 1, function(r) (r - mean(r, na.rm = TRUE)) / max(sd(r, na.rm = TRUE), 1e-9)))
  if (nrow(sig) == 1) z <- matrix(z, nrow = 1)
  z <- pmax(pmin(z, 2.5), -2.5)
  rlab <- sprintf("%s (%s) %s", sig$label, sig$metacluster, sig$marker)
  h <- max(3.6, 1.9 + 0.17 * nrow(z))
  w <- 3.6 + 0.24 * ncol(z)
  open_png(path, w, h)
  layout(matrix(c(1, 2), 1), widths = c(w - 1.0, 1.0))
  par(mar = c(5.5, max(nchar(rlab)) * 0.42 + 1, 2.6, 0.5))
  image(seq_len(ncol(z)), seq_len(nrow(z)), t(z[nrow(z):1, , drop = FALSE]), col = diverging(), zlim = c(-2.5, 2.5),
        axes = FALSE, xlab = "", ylab = "", useRaster = FALSE)
  axis(2, at = seq_len(nrow(z)), labels = rev(rlab), lwd = 0, cex.axis = 0.7)
  axis(1, at = seq_len(ncol(z)), labels = as.character(info$sample_id), lwd = 0, las = 2, cex.axis = 0.6)
  cols <- group_colors(study$levels)
  y0 <- nrow(z) + 0.6
  rect(seq_len(ncol(z)) - 0.5, y0, seq_len(ncol(z)) + 0.5, y0 + max(0.5, nrow(z) * 0.03),
       col = cols[as.character(info$condition)], border = PAL$surface, xpd = NA)
  title(main = sprintf("Significant cluster-marker pairs (FDR < %g)%s", fdr,
                       if (n_sig > nrow(sig)) sprintf(", top %d of %d", nrow(sig), n_sig) else ""),
        cex.main = 0.9, font.main = 1, line = 1.4)
  par(mar = c(5.5, 0.2, 2.6, 2.6))
  plot.new()
  legend("top", legend = names(cols), fill = cols, border = NA, bty = "n", cex = 0.7, text.col = PAL$ink2)
  zs <- seq(-2.5, 2.5, length.out = 101)
  ry <- c(0.08, 0.55)
  rect(0.35, head(seq(ry[1], ry[2], length.out = 102), -1), 0.6, tail(seq(ry[1], ry[2], length.out = 102), -1),
       col = diverging(), border = NA)
  text(0.65, c(ry[1], mean(ry), ry[2]), c("-2.5", "0", "2.5"), adj = 0, cex = 0.65, col = PAL$ink2)
  text(0.47, ry[2] + 0.04, "row z-score", cex = 0.65, col = PAL$ink2)
  dev.off()
  invisible(nrow(sig))
}

# Clusters x state markers: logFC per comparison, significant pairs ringed.
plot_ds_logfc <- function(study, ds, path) {
  cmps <- unique(ds$comparison)
  markers <- unique(ds$marker)
  ids <- study$clusters
  rlab <- sprintf("%s (%s)", study$label[ids], ids)
  lim <- max(abs(ds$logFC[ds$tested]), 0.1, na.rm = TRUE)
  w <- 1.2 + max(nchar(rlab)) * 0.075 + length(cmps) * (0.32 * length(markers) + 0.4)
  open_png(path, w, 1.8 + 0.26 * length(ids))
  par(mfrow = c(1, length(cmps)), mar = c(5, max(nchar(rlab)) * 0.45 + 1, 2.6, 0.5))
  for (cm in cmps) {
    d <- ds[ds$comparison == cm, ]
    m <- matrix(NA_real_, length(ids), length(markers), dimnames = list(ids, markers))
    m[cbind(d$metacluster, d$marker)] <- d$logFC
    s <- matrix(FALSE, length(ids), length(markers), dimnames = list(ids, markers))
    s[cbind(d$metacluster, d$marker)] <- d$significant
    zcol <- diverging()
    image(seq_along(markers), seq_along(ids), t(m[rev(ids), , drop = FALSE]), col = zcol, zlim = c(-lim, lim),
          axes = FALSE, xlab = "", ylab = "", main = cm, font.main = 1, cex.main = 0.9)
    na <- which(is.na(t(m[rev(ids), , drop = FALSE])), arr.ind = TRUE)
    if (nrow(na)) rect(na[, 1] - 0.5, na[, 2] - 0.5, na[, 1] + 0.5, na[, 2] + 0.5, col = PAL$grid, border = NA)
    sg <- which(t(s[rev(ids), , drop = FALSE]), arr.ind = TRUE)
    if (nrow(sg)) points(sg[, 1], sg[, 2], pch = 21, cex = 1.1, col = PAL$ink, lwd = 1.4)
    axis(1, at = seq_along(markers), labels = markers, las = 2, lwd = 0, cex.axis = 0.7)
    axis(2, at = seq_along(ids), labels = rev(rlab), lwd = 0, cex.axis = 0.7)
    mtext(sprintf("logFC, colour range +/-%.2g; ring = significant; gray = not tested", lim),
          side = 1, line = 3.8, cex = 0.55, col = PAL$ink2)
  }
  dev.off()
}
