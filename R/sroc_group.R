# ============================================================================
# sroc_group.R  -  Subgroup layer of the sROC plots
#
# When dta_sroc() / dta_sroc_pair() receive a `group` (one subgroup level
# per study) the figure changes shape:
#   * the in-plot summary box is dropped;
#   * one sROC curve + summary point per subgroup level is drawn in the
#     level's colour (each level refitted on its own studies);
#   * no legend: each level's symbol sits in its table header;
#   * beneath the panel: a centred title "By <group.name> (p = X)", X being
#     the test for subgroup differences, over one summary table per level,
#     spanning exactly the panel width (level name + symbol in a white header row).
#   * pairwise: one such figure per arm, side by side; p is per arm.
#
# Test for subgroup differences
#   single ..... likelihood-ratio test of the bivariate meta-regression with
#                the subgroup as covariate (Cochrane Handbook DTA ch. 10;
#                dta_compare_tests(), overall 2-df test), as in dta_forest().
#   pairwise ... the same per-arm test on each arm's studies (one figure
#                per arm).  .sroc_group_test_pair() (2-df test x subgroup
#                interaction on model B) is kept for a future joint p.
# ============================================================================

# Subgroup value per study from a vector named by study label or given in
# study order; NA marks studies outside every subgroup.
.sroc_group_values <- function(group, studlab) {
  if (is.null(group)) return(rep(NA_character_, length(studlab)))
  g <- if (!is.null(names(group))) as.character(group[studlab]) else as.character(group)
  if (length(g) != length(studlab))
    stop("`group` must have one entry per study in the fit.")
  g
}

# Refit a dta_single per subgroup level (levels with fewer than two studies
# are skipped and come back NULL).
.sroc_subgroup_fits <- function(fit, grp, levs, conf) {
  long  <- fit$long
  studs <- sort(unique(as.character(long$studlab)))
  nAGQ  <- fit$call_args$nAGQ
  if (is.null(nAGQ)) nAGQ <- 1L
  fits <- lapply(levs, function(lv) {
    keep <- studs[!is.na(grp) & grp == lv]
    if (length(keep) < 2) return(NULL)
    tryCatch(suppressWarnings(
      dta_fit_single(long[long$studlab %in% keep, , drop = FALSE],
                     wide = FALSE, nAGQ = nAGQ, conf = conf)),
      error = function(e) NULL)
  })
  stats::setNames(fits, levs)
}

# Harbord sROC curve + summary point of every non-NULL subgroup fit.
.sroc_subgroup_curves <- function(fits, fpr_grid) {
  lsp <- -stats::qlogis(fpr_grid)
  curves <- NULL; summ <- NULL
  for (lv in names(fits)) {
    f0 <- fits[[lv]]
    if (is.null(f0)) next
    f  <- .fixed_se_sp(f0$fit)
    g  <- .arm_auc_geom(f0, length(fpr_grid))
    curves <- rbind(curves, data.frame(
      grp = lv, fpr = fpr_grid,
      tpr = stats::plogis(f$lsens - g$slope * (lsp - f$lspec)),
      stringsAsFactors = FALSE))
    summ <- rbind(summ, data.frame(
      grp = lv, fpr = 1 - stats::plogis(f$lspec), tpr = stats::plogis(f$lsens),
      stringsAsFactors = FALSE))
  }
  list(curves = curves, summary = summ)
}

# -- Tests for subgroup differences ----------------------------------------

# Single test: overall LR test of the subgroup covariate (needs exactly two
# levels with at least two studies each; NA otherwise).
.sroc_group_test_single <- function(fit, grp, conf) {
  long  <- fit$long
  studs <- sort(unique(as.character(long$studlab)))
  cnt   <- .extract_counts(fit, NULL)
  ok    <- !is.na(grp)
  levs  <- unique(grp[ok])
  if (length(levs) != 2 || any(table(grp[ok]) < 2)) return(NA_real_)
  wd <- data.frame(studlab = studs[ok], TP = cnt$TP[ok], FP = cnt$FP[ok],
                   FN = cnt$FN[ok], TN = cnt$TN[ok], subgroup = grp[ok],
                   stringsAsFactors = FALSE)
  tryCatch(suppressWarnings(
    dta_compare_tests(wd, test_var = "subgroup", conf = conf)$compare$
      lr_tests$p_value[1]),
    error = function(e) NA_real_)
}

# Pairwise: 2-df LR test of the test x subgroup interaction on model B.
# `grp` is named by study label.
.sroc_group_test_pair <- function(long, tv, grp, nAGQ = 1L) {
  g <- as.character(grp[as.character(long$studlab)])
  d <- long[!is.na(g), , drop = FALSE]
  g <- g[!is.na(g)]
  levs_g <- unique(g)
  levs_t <- sort(unique(as.character(d[[tv]])))
  if (length(levs_g) != 2 || length(levs_t) != 2) return(NA_real_)
  tB <- as.integer(as.character(d[[tv]]) == levs_t[2])
  g2 <- as.integer(g == levs_g[2])
  d$seA <- d$sens * (1 - tB); d$seB <- d$sens * tB
  d$spA <- d$spec * (1 - tB); d$spB <- d$spec * tB
  d$se_g  <- d$sens * g2;     d$sp_g  <- d$spec * g2
  d$seB_g <- d$seB  * g2;     d$spB_g <- d$spB  * g2
  re <- "(0 + sens + spec | studlab)"
  f0 <- stats::as.formula(paste(
    "cbind(true, n - true) ~ 0 + seA + seB + spA + spB + se_g + sp_g +", re))
  f1 <- stats::as.formula(paste(
    "cbind(true, n - true) ~ 0 + seA + seB + spA + spB + se_g + sp_g +",
    "seB_g + spB_g +", re))
  tryCatch(suppressWarnings({
    m0 <- lme4::glmer(f0, data = d, family = stats::binomial, nAGQ = nAGQ)
    m1 <- lme4::glmer(f1, data = d, family = stats::binomial, nAGQ = nAGQ)
    lmtest::lrtest(m0, m1)[["Pr(>Chisq)"]][2]
  }), error = function(e) NA_real_)
}

# -- Summary tables ----------------------------------------------------------

.sroc_fmt_ci <- function(e, l, u, digits = 2) {
  if (is.na(e)) return("n/a")
  if (is.na(l) || is.na(u)) return(sprintf(sprintf("%%.%df", digits), e))
  sprintf(sprintf("%%.%df (%%.%df, %%.%df)", digits, digits, digits), e, l, u)
}

.sroc_fmt_p <- function(p) {
  if (is.null(p) || is.na(p)) return("n/a")
  if (p < 0.001) "<0.001" else sprintf("%.3f", p)
}

# Single test, one subgroup level: Studies, Sens, Spec, AUC, LR+, LR-, I2.
.sroc_group_table_single <- function(fit, k, conf, auc_ci, B, digits = 2,
                                     lr.show = FALSE) {
  r  <- .loo_single_row(fit, conf, auc_ci, B)$row
  i2 <- if (!is.null(fit) && !is.null(fit$heterogeneity))
    sprintf("%.1f%%", 100 * fit$heterogeneity$I2_biv) else "n/a"
  out <- data.frame(
    Measure = c("Studies", "Sensitivity", "Specificity", "AUC", "LR+", "LR-",
                "I²"),
    Estimate = c(
      as.character(k),
      .sroc_fmt_ci(r$sens, r$sens_lci, r$sens_uci, digits),
      .sroc_fmt_ci(r$spec, r$spec_lci, r$spec_uci, digits),
      .sroc_fmt_ci(r$auc,  r$auc_lci,  r$auc_uci,  digits),
      .sroc_fmt_ci(r$lrp,  r$lrp_lci,  r$lrp_uci,  digits),
      .sroc_fmt_ci(r$lrn,  r$lrn_lci,  r$lrn_uci,  digits),
      i2),
    stringsAsFactors = FALSE, check.names = FALSE)
  # LR+/LR- are always computed; shown only on request
  if (!isTRUE(lr.show)) out <- out[!out$Measure %in% c("LR+", "LR-"), ]
  rownames(out) <- NULL
  out
}

# Pairwise, one subgroup level: rows Sens / Spec / AUC; columns arm.e,
# arm.c, Diff (.e - .c) and p.  The comparison is refitted on the level's
# studies only.
.sroc_group_table_pair <- function(x, studs, arm.e, arm.c,
                                   label.e, label.c, conf, auc_ci, B,
                                   digits = 2) {
  long <- x$long[x$long$studlab %in% studs, , drop = FALSE]
  tv   <- x$test_var
  nAGQ <- x$pair$call_args$nAGQ; if (is.null(nAGQ)) nAGQ <- 1L
  variance <- x$pair$variance;  if (is.null(variance)) variance <- "equal"
  n_arm <- vapply(c(arm.e, arm.c), function(a)
    length(unique(long$studlab[long[[tv]] == a])), integer(1))
  res <- if (all(n_arm >= 2)) tryCatch(suppressWarnings(list(
    pair = dta_fit_pairwise(long, test_var = tv, variance = variance,
                            nAGQ = nAGQ),
    e = dta_fit_single(long[long[[tv]] == arm.e, , drop = FALSE],
                       wide = FALSE, nAGQ = nAGQ, conf = conf),
    c = dta_fit_single(long[long[[tv]] == arm.c, , drop = FALSE],
                       wide = FALSE, nAGQ = nAGQ, conf = conf))),
    error = function(e) NULL) else NULL
  d  <- .loo_pair_row(res$pair, res$e, res$c, arm.e, conf, auc_ci, B)$row
  re <- .loo_single_row(res$e, conf, auc_ci, B)$row
  rc <- .loo_single_row(res$c, conf, auc_ci, B)$row
  out <- data.frame(
    Measure = c("Sensitivity", "Specificity", "AUC"),
    e = c(.sroc_fmt_ci(re$sens, re$sens_lci, re$sens_uci, digits),
          .sroc_fmt_ci(re$spec, re$spec_lci, re$spec_uci, digits),
          .sroc_fmt_ci(re$auc,  re$auc_lci,  re$auc_uci,  digits)),
    c = c(.sroc_fmt_ci(rc$sens, rc$sens_lci, rc$sens_uci, digits),
          .sroc_fmt_ci(rc$spec, rc$spec_lci, rc$spec_uci, digits),
          .sroc_fmt_ci(rc$auc,  rc$auc_lci,  rc$auc_uci,  digits)),
    Diff = c(.sroc_fmt_ci(d$dsens, d$dsens_lci, d$dsens_uci, digits),
             .sroc_fmt_ci(d$dspec, d$dspec_lci, d$dspec_uci, digits),
             .sroc_fmt_ci(d$dauc,  d$dauc_lci,  d$dauc_uci,  digits)),
    p = c(.sroc_fmt_p(d$p_sens), .sroc_fmt_p(d$p_spec), .sroc_fmt_p(d$p_auc)),
    stringsAsFactors = FALSE, check.names = FALSE)
  names(out) <- c("Measure", label.e, label.c, "Diff", "P-value")
  out
}

# -- Grobs ------------------------------------------------------------------

# Bordered tableGrob (same look as the dta_sroc_pair differences table);
# an optional bold `title` row spans the columns inside the border.
.sroc_table_grob <- function(df, title = NULL, cex = 0.75,
                             pch = NULL, col = NULL) {
  # zebra body rows (fill recycles down each column), light column header
  zebra <- rep(c("grey93", "white"), length.out = nrow(df))
  tg <- gridExtra::tableGrob(
    df, rows = NULL,
    theme = gridExtra::ttheme_minimal(
      core    = list(fg_params = list(cex = cex),
                     bg_params = list(fill = zebra, col = NA)),
      colhead = list(fg_params = list(cex = cex, fontface = "bold"),
                     bg_params = list(fill = "grey85", col = NA))))
  if (!is.null(title)) {
    # level name as a dark header band with white text, inside the border
    tg <- gtable::gtable_add_rows(tg, grid::unit(18, "pt"), pos = 0)
    tg <- gtable::gtable_add_grob(
      tg, grobs = grid::rectGrob(gp = grid::gpar(fill = "white", col = NA)),
      t = 1, b = 1, l = 1, r = ncol(tg), name = "level-band")
    ttl <- grid::textGrob(title, gp = grid::gpar(fontface = "bold",
                                                 fontsize = 10, col = "black"))
    band <- if (!is.null(pch)) {
      # the level's plot symbol (as in the panel) just left of the name
      sym <- grid::pointsGrob(
        x = grid::unit(0.5, "npc") - grid::grobWidth(ttl) * 0.5 - grid::unit(9, "pt"),
        y = grid::unit(0.5, "npc"), pch = pch, size = grid::unit(9, "pt"),
        gp = grid::gpar(col = if (is.null(col)) "black" else col, lwd = 2.2))
      grid::gTree(children = grid::gList(sym, ttl))
    } else ttl
    tg <- gtable::gtable_add_grob(tg, grobs = band, t = 1, b = 1, l = 1,
                                  r = ncol(tg), name = "level-title")
  }
  gtable::gtable_add_grob(
    tg, grobs = grid::rectGrob(gp = grid::gpar(fill = NA, col = "black",
                                               lwd = 1.2)),
    t = 1, b = nrow(tg), l = 1, r = ncol(tg), name = "container-border")
}

# Subgroup figure: the sROC panel on top and, beneath it, the title
# "By <name> (p = X)" over the per-level tables (level name inside each
# border).  The tables row is placed in the same layout columns as the
# ggplot panel and its columns are relative units, so the first table
# starts at the panel's left edge, the last ends at its right edge, and all
# width but the gaps is used.  The panel is `panel.scale` times the natural
# width of the tables.  Returns a gtable with "height" / "width" attributes.
.sroc_group_figure <- function(p, tables, group.name, p_value,
                               panel.scale = 1.4, gap_in = 0.3,
                               shapes = NULL, colors = NULL) {
  title_txt <- sprintf("By %s (p %s)", group.name,
                       if (is.na(p_value)) "= n/a"
                       else if (p_value < 0.001) "< 0.001"
                       else sprintf("= %.3f", p_value))
  title <- grid::textGrob(title_txt, gp = grid::gpar(fontface = "bold",
                                                     fontsize = 12))
  grobs <- lapply(names(tables), function(lv)
    .sroc_table_grob(tables[[lv]], lv,
                     pch = if (!is.null(shapes)) shapes[[lv]] else NULL,
                     col = if (!is.null(colors)) colors[[lv]] else NULL))
  nat_w <- vapply(grobs, function(g)
    grid::convertWidth(sum(g$widths), "in", valueOnly = TRUE), numeric(1))
  tab_h <- Reduce(function(a, b) max(a, b), lapply(grobs, function(g) sum(g$heights)))
  n     <- length(grobs)

  # stretchable tables: column widths become proportions of the cell
  grobs <- lapply(grobs, function(g) {
    w <- grid::convertWidth(g$widths, "in", valueOnly = TRUE)
    g$widths <- grid::unit(w / sum(w), "null")
    g
  })
  cells <- list(); widths <- list()
  for (i in seq_len(n)) {
    if (i > 1) {
      cells  <- c(cells, list(grid::nullGrob()))
      widths <- c(widths, list(grid::unit(gap_in, "in")))
    }
    cells  <- c(cells, list(grobs[[i]]))
    widths <- c(widths, list(grid::unit(nat_w[i], "null")))
  }
  row <- gridExtra::arrangeGrob(grobs = cells, ncol = length(cells),
                                widths = do.call(grid::unit.c, widths),
                                heights = tab_h)

  # figure width from the tables; room for the y axis on the left
  axis_in <- 0.75
  fig_w   <- grid::unit((sum(nat_w) + gap_in * (n - 1)) * panel.scale + axis_in, "in")
  pan_h   <- fig_w + grid::unit(1.3, "in")   # + title (3 lines) and x axis

  # tables in the ggplot's own panel columns so both share one x extent
  gp  <- ggplot2::ggplotGrob(p)
  pan <- gp$layout[gp$layout$name == "panel", ][1, ]
  title_h <- grid::unit(22, "pt")
  tb <- gtable::gtable(widths = gp$widths, heights = grid::unit.c(title_h, tab_h))
  tb <- gtable::gtable_add_grob(tb, title, t = 1, l = pan$l, r = pan$r,
                                name = "group-title")
  tb <- gtable::gtable_add_grob(tb, row,   t = 2, l = pan$l, r = pan$r,
                                name = "group-tables")

  g <- gridExtra::arrangeGrob(gp, tb, nrow = 2,
                              heights = grid::unit.c(pan_h, title_h + tab_h),
                              widths  = fig_w)
  attr(g, "height") <- pan_h + title_h + tab_h + grid::unit(12, "pt")
  attr(g, "width")  <- fig_w
  g
}

#' @export
print.dta_sroc_group <- function(x, ...) {
  grid::grid.newpage()
  h <- attr(x, "height"); w <- attr(x, "width")
  if (!is.null(h)) {
    if (is.null(w)) w <- grid::unit(1, "npc")
    grid::pushViewport(grid::viewport(height = h, width = w))
    grid::grid.draw(x)
    grid::popViewport()
  } else {
    grid::grid.draw(x)
  }
  invisible(x)
}
