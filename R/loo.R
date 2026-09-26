# ============================================================================
# loo.R  -  Leave-one-out sensitivity analysis
#
# For every study the model is refitted with that study removed, and the
# pooled measures are recomputed, in the spirit of meta::metainf().
#
#   * dta_single ........... Sens, Spec, AUC, LR+, LR- (each with CI) after
#                            omitting each study, plus the pooled estimate.
#   * dta_pairwise_result .. the differences (.e - .c) in Sens, Spec, AUC,
#                            LR+ and LR- with CI and p-value after omitting
#                            each study (both arms of a paired design), plus
#                            the full-data comparison.
#
# Inference matches the rest of the package:
#   Sens / Spec ....... Wald CI on the logit scale (Cochrane Appendix 5);
#                       pairwise p from the likelihood-ratio ladder
#                       (Cochrane Appendix 12, via dta_compare()).
#   AUC ............... trapezoidal integral of the Harbord sROC; CI and
#                       the dAUC p-value from the parametric MVN bootstrap
#                       used by dta_sroc() / dta_sroc_pair().
#   LR+ / LR- ......... delta method (msm::deltamethod); pairwise diff CI
#                       and Wald p from the joint 4-parameter model B.
#
# dta_loo_forest() draws the meta-style forest (one row per omitted
# study, letter-coded) and, beneath it, an unlabelled sROC panel with one
# curve per omission; each study point is drawn as its letter in the
# colour of the curve that omits it.
# ============================================================================

#' Leave-one-out sensitivity analysis
#'
#' Refits the model once per study with that study removed and collects
#' the pooled measures each time, like `meta::metainf()`.  For a single
#' test the table holds sensitivity, specificity, AUC, LR+ and LR- (each
#' with a confidence interval).  For a pairwise comparison it holds the
#' difference `.e - .c` in each of those measures with its confidence
#' interval and p-value.  The last row is the estimate from all studies.
#'
#' @param x       A `dta_single` fit or a `dta_pairwise_result` (from
#'   [dta_pairwise()] / [dta_compare_tests()]).  In a paired design both
#'   arms of the omitted study are dropped together.
#' @param conf    Confidence level (default 0.95).
#' @param auc_ci  Logical. Bootstrap a CI for the AUC (and, for a
#'   comparison, the dAUC CI and p-value)?  Default `TRUE`; `FALSE` keeps
#'   only the AUC point estimates and is much faster.
#' @param B       MVN-bootstrap replicates per refit (default 2000).
#' @param verbose Logical. Print one line per refit (default `FALSE`).
#'
#' @return An S3 object of class `"dta_loo"` with `$table` (one row per
#'   omitted study plus the pooled row, flagged by `$table$pooled`),
#'   `$type` (`"single"` or `"pair"`), `$letters` (the row codes used by
#'   the plot) and the sROC geometry of every refit.  `print()` shows the
#'   table; `plot()` / [dta_loo_forest()] draws the forest and sROC figure.
#' @examples
#' \donttest{
#' data(anti_ccp2)
#' fit <- dta_fit_single(anti_ccp2, wide = TRUE)
#' loo <- dta_loo(fit, auc_ci = FALSE)
#' print(loo)
#' plot(loo)
#' }
#' @export
dta_loo <- function(x, conf = 0.95, auc_ci = TRUE, B = 2000,
                    verbose = FALSE) {
  if (inherits(x, "dta_single")) {
    .loo_single(x, conf, auc_ci, B, verbose)
  } else if (inherits(x, "dta_pairwise_result")) {
    .loo_pair(x, conf, auc_ci, B, verbose)
  } else {
    stop("`x` must be a dta_single fit or a dta_pairwise_result.")
  }
}

# Row codes: A..Z, then AA, AB, ...
.loo_letters <- function(k) {
  codes <- c(LETTERS, as.vector(t(outer(LETTERS, LETTERS, paste0))))
  codes[seq_len(k)]
}

# -- Single test ------------------------------------------------------------

.loo_single <- function(fit, conf, auc_ci, B, verbose) {
  long  <- fit$long
  studs <- unique(as.character(long$studlab))
  k     <- length(studs)
  if (k < 3) stop("Leave-one-out needs at least three studies.")
  nAGQ  <- fit$call_args$nAGQ
  if (is.null(nAGQ)) nAGQ <- 1L

  rows <- vector("list", k)
  geom <- vector("list", k)
  for (i in seq_len(k)) {
    if (verbose) message("Omitting ", studs[i], " (", i, "/", k, ")")
    sub <- long[long$studlab != studs[i], , drop = FALSE]
    f_i <- tryCatch(
      suppressWarnings(dta_fit_single(sub, wide = FALSE, nAGQ = nAGQ,
                                      conf = conf)),
      error = function(e) NULL)
    r <- .loo_single_row(f_i, conf, auc_ci, B)
    rows[[i]] <- r$row
    geom[[i]] <- r$geom
  }
  full <- .loo_single_row(fit, conf, auc_ci, B)

  tab <- rbind(do.call(rbind, rows), full$row)
  tab <- cbind(
    data.frame(letter  = c(.loo_letters(k), ""),
               studlab = c(studs, "Pooled estimate"),
               pooled  = c(rep(FALSE, k), TRUE),
               stringsAsFactors = FALSE),
    tab)
  rownames(tab) <- NULL

  out <- list(type = "single", table = tab, letters = .loo_letters(k),
              studies = studs, geom = geom, geom_full = full$geom,
              points = .loo_points(fit), conf = conf, auc_ci = auc_ci)
  class(out) <- "dta_loo"
  out
}

# One row of pooled measures (+ sROC geometry) for a dta_single fit; all
# NA when the refit failed.
.loo_single_row <- function(fit, conf, auc_ci, B, n_grid = 200) {
  na  <- NA_real_
  row <- data.frame(
    sens = na, sens_lci = na, sens_uci = na,
    spec = na, spec_lci = na, spec_uci = na,
    auc  = na, auc_lci  = na, auc_uci  = na,
    lrp  = na, lrp_lci  = na, lrp_uci  = na,
    lrn  = na, lrn_lci  = na, lrn_uci  = na)
  if (is.null(fit)) return(list(row = row, geom = NULL))

  f   <- .fixed_se_sp(fit$fit)
  se  <- .logit_ci(f$lsens, f$se_lsens, conf)
  sp  <- .logit_ci(f$lspec, f$se_lspec, conf)
  der <- .derived_from_logits(f$lsens, f$lspec, f$vcov_fixed, conf)
  g   <- .arm_auc_geom(fit, n_grid)
  aci <- if (isTRUE(auc_ci)) .trapz_auc_ci(fit, g$slope, g$fpr_grid, B, conf)
         else NULL
  if (is.null(aci)) aci <- c(na, na)
  lrp <- der[der$measure == "LR+", ]
  lrn <- der[der$measure == "LR-", ]

  row[] <- list(
    se[["estimate"]], se[["lci"]], se[["uci"]],
    sp[["estimate"]], sp[["lci"]], sp[["uci"]],
    g$AUC, aci[1], aci[2],
    lrp$estimate, lrp$lci, lrp$uci,
    lrn$estimate, lrn$lci, lrn$uci)
  list(row  = row,
       geom = list(lsens = f$lsens, lspec = f$lspec, slope = g$slope))
}

# Study points (FPR, TPR) of a dta_single fit, in study order.
.loo_points <- function(fit) {
  long  <- fit$long
  studs <- unique(as.character(long$studlab))
  rs <- long[long$sens == 1, , drop = FALSE]
  rp <- long[long$spec == 1, , drop = FALSE]
  is <- match(studs, as.character(rs$studlab))
  ip <- match(studs, as.character(rp$studlab))
  data.frame(studlab = studs,
             tpr = rs$true[is] / rs$n[is],
             fpr = 1 - rp$true[ip] / rp$n[ip],
             stringsAsFactors = FALSE)
}

# -- Pairwise comparison ----------------------------------------------------

.loo_pair <- function(x, conf, auc_ci, B, verbose) {
  long  <- x$long
  tv    <- x$test_var
  arm.e <- x$labels$intervention
  arm.c <- x$labels$control
  studs <- unique(as.character(long$studlab))
  k     <- length(studs)
  if (k < 3) stop("Leave-one-out needs at least three studies.")
  nAGQ     <- x$pair$call_args$nAGQ
  variance <- x$pair$variance
  if (is.null(nAGQ)) nAGQ <- 1L
  if (is.null(variance)) variance <- "equal"

  arm_fit <- function(sub, arm) {
    dta_fit_single(sub[sub[[tv]] == arm, , drop = FALSE], wide = FALSE,
                   nAGQ = nAGQ, conf = conf)
  }

  rows <- vector("list", k)
  geom <- vector("list", k)
  for (i in seq_len(k)) {
    if (verbose) message("Omitting ", studs[i], " (", i, "/", k, ")")
    sub <- long[long$studlab != studs[i], , drop = FALSE]
    n_lev <- vapply(c(arm.e, arm.c), function(a)
      length(unique(sub$studlab[sub[[tv]] == a])), integer(1))
    res_i <- if (all(n_lev >= 2)) tryCatch(suppressWarnings(list(
      pair = dta_fit_pairwise(sub, test_var = tv, variance = variance,
                              nAGQ = nAGQ),
      e    = arm_fit(sub, arm.e),
      c    = arm_fit(sub, arm.c))), error = function(e) NULL) else NULL
    r <- .loo_pair_row(res_i$pair, res_i$e, res_i$c, arm.e, conf, auc_ci, B)
    rows[[i]] <- r$row
    geom[[i]] <- r$geom
  }
  full <- .loo_pair_row(x$pair, x$arms[[arm.e]], x$arms[[arm.c]], arm.e,
                        conf, auc_ci, B)

  tab <- rbind(do.call(rbind, rows), full$row)
  tab <- cbind(
    data.frame(letter  = c(.loo_letters(k), ""),
               studlab = c(studs, "Pooled estimate"),
               pooled  = c(rep(FALSE, k), TRUE),
               stringsAsFactors = FALSE),
    tab)
  rownames(tab) <- NULL

  out <- list(type = "pair", table = tab, letters = .loo_letters(k),
              studies = studs, geom = geom, geom_full = full$geom,
              points = list(e = .loo_points(x$arms[[arm.e]]),
                            c = .loo_points(x$arms[[arm.c]])),
              arms = list(e = arm.e, c = arm.c),
              conf = conf, auc_ci = auc_ci)
  class(out) <- "dta_loo"
  out
}

# One row of differences (.e - .c) with CI and p-value from a dta_pairwise
# fit plus the two per-arm dta_single fits; all NA when a refit failed.
.loo_pair_row <- function(pair, fit_e, fit_c, arm.e, conf, auc_ci, B,
                          n_grid = 200) {
  na  <- NA_real_
  row <- data.frame(
    dsens = na, dsens_lci = na, dsens_uci = na, p_sens = na,
    dspec = na, dspec_lci = na, dspec_uci = na, p_spec = na,
    dauc  = na, dauc_lci  = na, dauc_uci  = na, p_auc  = na,
    dlrp  = na, dlrp_lci  = na, dlrp_uci  = na, p_lrp  = na,
    dlrn  = na, dlrn_lci  = na, dlrn_uci  = na, p_lrn  = na)
  if (is.null(pair) || is.null(fit_e) || is.null(fit_c))
    return(list(row = row, geom = NULL))

  z   <- stats::qnorm(1 - (1 - conf) / 2)
  M   <- pair$models$B
  nm  <- c("seA", "seB", "spA", "spB")
  th  <- stats::setNames(summary(M)$coefficients[nm, 1], nm)
  V   <- as.matrix(stats::vcov(M))[nm, nm]
  # .e - .c: the pairwise levels are sorted alphabetically, so flip when the
  # intervention arm is level B.
  s   <- if (identical(pair$levels[1], arm.e)) 1 else -1
  pl  <- stats::plogis

  # Sens / Spec: delta-method CI and LR-test p (Cochrane Appendix 12).
  cmp    <- dta_compare(pair, conf = conf)
  d      <- cmp$differences
  lr     <- cmp$lr_tests
  flip   <- function(r) c(s * r$estimate,
                          if (s == 1) r$lci else -r$uci,
                          if (s == 1) r$uci else -r$lci)
  dse    <- flip(d[d$measure == "Absolute diff Sens", ][1, ])
  dsp    <- flip(d[d$measure == "Absolute diff Spec", ][1, ])
  p_sens <- lr$p_value[grepl("Se", lr$comparison)][1]
  p_spec <- lr$p_value[grepl("Sp", lr$comparison)][1]

  # LR+ / LR-: difference of the two arms' ratios, delta method on the joint
  # 4-parameter model, Wald p.
  pA <- "(exp(x1)/(1+exp(x1)))"; pB <- "(exp(x2)/(1+exp(x2)))"
  qA <- "(exp(x3)/(1+exp(x3)))"; qB <- "(exp(x4)/(1+exp(x4)))"
  f_lrp <- stats::as.formula(sprintf("~ %s/(1-%s) - %s/(1-%s)", pA, qA, pB, qB))
  f_lrn <- stats::as.formula(sprintf("~ (1-%s)/%s - (1-%s)/%s", pA, qA, pB, qB))
  lrp_A <- pl(th["seA"]) / (1 - pl(th["spA"]))
  lrp_B <- pl(th["seB"]) / (1 - pl(th["spB"]))
  lrn_A <- (1 - pl(th["seA"])) / pl(th["spA"])
  lrn_B <- (1 - pl(th["seB"])) / pl(th["spB"])
  dlrp  <- s * unname(lrp_A - lrp_B)
  dlrn  <- s * unname(lrn_A - lrn_B)
  se_dlrp <- msm::deltamethod(f_lrp, mean = unname(th), cov = V)
  se_dlrn <- msm::deltamethod(f_lrn, mean = unname(th), cov = V)
  wald_p  <- function(est, se) 2 * stats::pnorm(-abs(est / se))

  # AUC: per-arm trapezoidal AUC; dAUC CI and p from the MVN bootstrap.
  if (isTRUE(auc_ci)) {
    ap   <- .compute_auc_pair(fit_e, fit_c, B = B, conf = conf, n_grid = n_grid)
    dauc <- c(ap$diff$est, ap$diff$ci, ap$diff$p)
  } else {
    ge <- .arm_auc_geom(fit_e, n_grid); gc <- .arm_auc_geom(fit_c, n_grid)
    dauc <- c(ge$AUC - gc$AUC, na, na, na)
  }

  row[] <- list(
    dse[1], dse[2], dse[3], p_sens,
    dsp[1], dsp[2], dsp[3], p_spec,
    dauc[1], dauc[2], dauc[3], dauc[4],
    dlrp, dlrp - z * se_dlrp, dlrp + z * se_dlrp, wald_p(dlrp, se_dlrp),
    dlrn, dlrn - z * se_dlrn, dlrn + z * se_dlrn, wald_p(dlrn, se_dlrn))

  arm_geom <- function(fit) {
    f <- .fixed_se_sp(fit$fit)
    list(lsens = f$lsens, lspec = f$lspec, slope = .arm_auc_geom(fit)$slope)
  }
  list(row = row, geom = list(e = arm_geom(fit_e), c = arm_geom(fit_c)))
}

# -- print ------------------------------------------------------------------

.loo_fmt_ci <- function(e, l, u, digits) {
  fstr <- sprintf("%%.%df (%%.%df, %%.%df)", digits, digits, digits)
  ifelse(is.na(e), "n/a",
         ifelse(is.na(l) | is.na(u), sprintf(sprintf("%%.%df", digits), e),
                sprintf(fstr, e, l, u)))
}

.loo_fmt_p <- function(p) {
  ifelse(is.na(p), "n/a", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

# Display table: one text column per measure (and per p-value).
.loo_display <- function(x, digits = 2) {
  t <- x$table
  lab <- ifelse(t$pooled, t$studlab, paste("Omitting", t$studlab))
  if (x$type == "single") {
    data.frame(
      Study  = lab,
      Sens   = .loo_fmt_ci(t$sens, t$sens_lci, t$sens_uci, digits),
      Spec   = .loo_fmt_ci(t$spec, t$spec_lci, t$spec_uci, digits),
      AUC    = .loo_fmt_ci(t$auc,  t$auc_lci,  t$auc_uci,  digits),
      `LR+`  = .loo_fmt_ci(t$lrp,  t$lrp_lci,  t$lrp_uci,  digits),
      `LR-`  = .loo_fmt_ci(t$lrn,  t$lrn_lci,  t$lrn_uci,  digits),
      check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    data.frame(
      Study    = lab,
      `dSens`  = .loo_fmt_ci(t$dsens, t$dsens_lci, t$dsens_uci, digits),
      `p`      = .loo_fmt_p(t$p_sens),
      `dSpec`  = .loo_fmt_ci(t$dspec, t$dspec_lci, t$dspec_uci, digits),
      `p `     = .loo_fmt_p(t$p_spec),
      `dAUC`   = .loo_fmt_ci(t$dauc,  t$dauc_lci,  t$dauc_uci,  digits),
      `p  `    = .loo_fmt_p(t$p_auc),
      `dLR+`   = .loo_fmt_ci(t$dlrp,  t$dlrp_lci,  t$dlrp_uci,  digits),
      `p   `   = .loo_fmt_p(t$p_lrp),
      `dLR-`   = .loo_fmt_ci(t$dlrn,  t$dlrn_lci,  t$dlrn_uci,  digits),
      `p    `  = .loo_fmt_p(t$p_lrn),
      check.names = FALSE, stringsAsFactors = FALSE)
  }
}

#' @export
print.dta_loo <- function(x, digits = 2, ...) {
  cat("<dta_loo>  Leave-one-out sensitivity analysis\n")
  if (x$type == "single") {
    cat("  Studies: ", length(x$studies), "   ", 100 * x$conf,
        "% CI in brackets\n\n", sep = "")
  } else {
    cat("  Arms: ", x$arms$e, " vs ", x$arms$c,
        "   differences are ", x$arms$e, " - ", x$arms$c, "\n",
        "  Studies: ", length(x$studies), "   ", 100 * x$conf,
        "% CI in brackets\n\n", sep = "")
  }
  print(.loo_display(x, digits), row.names = FALSE, right = FALSE)
  invisible(x)
}

#' @export
plot.dta_loo <- function(x, ...) dta_loo_forest(x, ...)

# -- Forest + sROC figure ---------------------------------------------------

#' Leave-one-out forest plot with sROC curves
#'
#' Draws the leave-one-out table of [dta_loo()] as a forest plot in the
#' style of `meta::forest(metainf())`: one row per omitted study
#' ("Omitting ..."), the pooled estimate last, a letter code on every row.
#' For a single test the sensitivity and specificity carry the CI panels
#' and AUC, LR+ and LR- are text columns; for a comparison the
#' differences in sensitivity and specificity carry the CI panels and each
#' difference is followed by its p-value.
#'
#' Beneath the forest, an unlabelled sROC panel (two for a comparison, one
#' per arm) draws one thin coloured curve per omission and the full-data
#' curve in black.  Each study point is drawn as its row letter, in the
#' colour of the curve fitted without that study, so an influential study
#' is read off directly: its letter sits away from the cloud and its curve
#' departs from the black one.
#'
#' @param x       A `dta_loo` object, or a `dta_single` /
#'   `dta_pairwise_result` (then [dta_loo()] is run first; `...` is
#'   forwarded to it).
#' @param sroc    Logical. Draw the sROC panel below the forest?
#'   Default `TRUE`.
#' @param digits  Display digits for the value columns (default 2).
#' @param title   Optional title, left-aligned above the forest.
#' @param row_in  Height of one forest row in inches (default 0.2).
#' @param sroc_in Height of the sROC band in inches (default 4 for one
#'   panel or two).
#' @param ...     Forwarded to [dta_loo()] when `x` is not a `dta_loo`.
#'
#' @return A `gtable` of class `"dta_loo_plot"`; its `print()` method draws
#'   it, so it renders when returned at top level.
#' @examples
#' \donttest{
#' data(anti_ccp2)
#' fit <- dta_fit_single(anti_ccp2, wide = TRUE)
#' dta_loo_forest(fit, auc_ci = FALSE)
#' }
#' @export
dta_loo_forest <- function(x, sroc = TRUE, digits = 2, title = NULL,
                           row_in = 0.2, sroc_in = NULL, ...) {
  if (!inherits(x, "dta_loo")) x <- dta_loo(x, ...)
  loo <- x
  t   <- loo$table
  n   <- nrow(t)
  k   <- n - 1
  conf_pct <- paste0(100 * loo$conf, "% CI")

  # Rows top to bottom: header, k omitted rows, pooled row.
  ypos     <- n + 1 - seq_len(n)
  header_y <- n + 1
  ylim     <- c(0.5, header_y + 0.5)
  face     <- ifelse(t$pooled, "bold", "plain")
  omit_y   <- ypos[!t$pooled]
  zebra_df <- data.frame(ypos = omit_y[seq_along(omit_y) %% 2L == 1L])

  hdr_size <- 3.2; cell_size <- 3.0; axis_text_size <- 8

  base_theme <- ggplot2::theme_bw() +
    ggplot2::theme(
      axis.title       = ggplot2::element_blank(),
      axis.text.y      = ggplot2::element_blank(),
      axis.ticks.y     = ggplot2::element_blank(),
      panel.grid       = ggplot2::element_blank(),
      panel.border     = ggplot2::element_blank(),
      plot.margin      = ggplot2::margin(t = 5.5, r = 0, b = 2, l = 0),
      legend.position  = "none")
  invisible_axis <- ggplot2::theme(
    axis.text.x  = ggplot2::element_text(colour = NA, size = axis_text_size),
    axis.ticks.x = ggplot2::element_line(colour = NA),
    axis.line.x  = ggplot2::element_line(colour = NA))
  visible_axis <- ggplot2::theme(
    axis.text.x  = ggplot2::element_text(colour = "black", size = axis_text_size),
    axis.ticks.x = ggplot2::element_line(colour = "black", linewidth = 0.4),
    axis.line.x  = ggplot2::element_line(colour = "black", linewidth = 0.4),
    axis.ticks.length = grid::unit(3, "pt"))
  zebra <- function() ggplot2::geom_rect(
    data = zebra_df, inherit.aes = FALSE,
    ggplot2::aes(xmin = -Inf, xmax = Inf, ymin = ypos - 0.5, ymax = ypos + 0.5),
    fill = "grey92")

  # ----- label panel: letter + "Omitting <study>" ---------------------------
  lab <- ifelse(t$pooled, t$studlab, paste("Omitting", t$studlab))
  letter_w <- max(nchar(t$letter)) + 2
  lab_w    <- max(nchar(c("Study", lab))) + 2
  tot_w    <- letter_w + lab_w
  x_letter <- 0
  x_lab    <- letter_w / tot_w
  cells <- rbind(
    data.frame(x = x_letter, ypos = ypos, label = t$letter, face = "bold",
               stringsAsFactors = FALSE),
    data.frame(x = x_lab, ypos = ypos, label = lab, face = face,
               stringsAsFactors = FALSE))
  cells <- cells[cells$label != "", , drop = FALSE]
  p_labels <- ggplot2::ggplot(cells, ggplot2::aes(y = ypos)) +
    zebra() +
    ggplot2::geom_text(ggplot2::aes(x = x, label = label, fontface = face),
                       hjust = 0, size = cell_size) +
    ggplot2::annotate("text", x = x_lab, y = header_y, label = "Study",
                      hjust = 0, fontface = "bold", size = hdr_size) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = ylim) +
    base_theme + invisible_axis
  label_w_in <- max(2.2, 3.0 * tot_w / 43)

  # ----- text and CI panels --------------------------------------------------
  text_panel <- function(txt, header, x = 0.5, hjust = 0.5) {
    df <- data.frame(ypos = ypos, txt = txt, face = face, stringsAsFactors = FALSE)
    ggplot2::ggplot(df, ggplot2::aes(y = ypos)) +
      zebra() +
      ggplot2::geom_text(ggplot2::aes(label = txt, fontface = face),
                         x = x, hjust = hjust, size = cell_size) +
      ggplot2::annotate("text", x = x, y = header_y, label = header,
                        hjust = hjust, fontface = "bold", size = hdr_size) +
      ggplot2::coord_cartesian(xlim = c(0, 1), ylim = ylim) +
      base_theme + invisible_axis
  }
  ci_panel <- function(est, lci, uci, xlim, breaks, ref = NULL) {
    df <- data.frame(ypos = ypos, est = est, lci = lci, uci = uci,
                     pooled = t$pooled)
    pooled_est <- est[t$pooled]
    df <- df[!is.na(df$est), , drop = FALSE]
    df$shape <- ifelse(df$pooled, 18L, 15L)
    df$size  <- ifelse(df$pooled, 4, 2)
    p <- ggplot2::ggplot(df, ggplot2::aes(x = est, y = ypos)) + zebra()
    if (!is.null(ref))
      p <- p + ggplot2::geom_vline(xintercept = ref, colour = "grey55",
                                   linewidth = 0.4)
    if (length(pooled_est) && !is.na(pooled_est))
      p <- p + ggplot2::geom_vline(xintercept = pooled_est, linetype = "dashed",
                                   colour = "grey40", linewidth = 0.4)
    p +
      ggplot2::geom_errorbar(ggplot2::aes(xmin = lci, xmax = uci),
                             width = 0.25, orientation = "y") +
      ggplot2::geom_point(ggplot2::aes(shape = shape, size = size)) +
      ggplot2::scale_shape_identity() + ggplot2::scale_size_identity() +
      ggplot2::scale_x_continuous(breaks = breaks) +
      ggplot2::coord_cartesian(xlim = xlim, ylim = ylim) +
      base_theme + visible_axis
  }
  # differences can be negative, so the comparison uses ", " between bounds
  sep <- if (loo$type == "single") "-" else ", "
  fmt <- function(e, l, u) {
    fstr <- sprintf("%%.%df (%%.%df%s%%.%df)", digits, digits, sep, digits)
    ifelse(is.na(e), "n/a",
           ifelse(is.na(l) | is.na(u), sprintf(sprintf("%%.%df", digits), e),
                  sprintf(fstr, e, l, u)))
  }
  diff_lim <- function(l, u) {
    r <- range(c(l, u, 0), na.rm = TRUE)
    r + c(-1, 1) * 0.06 * diff(r)
  }

  if (loo$type == "single") {
    panels <- list(
      p_labels,
      ci_panel(t$sens, t$sens_lci, t$sens_uci, c(0, 1), seq(0, 1, 0.2)),
      text_panel(fmt(t$sens, t$sens_lci, t$sens_uci), paste0("Sens (", conf_pct, ")")),
      ci_panel(t$spec, t$spec_lci, t$spec_uci, c(0, 1), seq(0, 1, 0.2)),
      text_panel(fmt(t$spec, t$spec_lci, t$spec_uci), paste0("Spec (", conf_pct, ")")),
      text_panel(fmt(t$auc, t$auc_lci, t$auc_uci), paste0("AUC (", conf_pct, ")")),
      text_panel(fmt(t$lrp, t$lrp_lci, t$lrp_uci), paste0("LR+ (", conf_pct, ")")),
      text_panel(fmt(t$lrn, t$lrn_lci, t$lrn_uci), paste0("LR- (", conf_pct, ")")))
    widths <- c(label_w_in, 1.6, 1.5, 1.6, 1.5, 1.5, 1.5, 1.5)
  } else {
    pv <- function(p) .loo_fmt_p(p)
    lim_se <- diff_lim(t$dsens_lci, t$dsens_uci)
    lim_sp <- diff_lim(t$dspec_lci, t$dspec_uci)
    panels <- list(
      p_labels,
      ci_panel(t$dsens, t$dsens_lci, t$dsens_uci, lim_se, ggplot2::waiver(), ref = 0),
      text_panel(fmt(t$dsens, t$dsens_lci, t$dsens_uci), paste0("ΔSens (", conf_pct, ")")),
      text_panel(pv(t$p_sens), "p"),
      ci_panel(t$dspec, t$dspec_lci, t$dspec_uci, lim_sp, ggplot2::waiver(), ref = 0),
      text_panel(fmt(t$dspec, t$dspec_lci, t$dspec_uci), paste0("ΔSpec (", conf_pct, ")")),
      text_panel(pv(t$p_spec), "p"),
      text_panel(fmt(t$dauc, t$dauc_lci, t$dauc_uci), paste0("ΔAUC (", conf_pct, ")")),
      text_panel(pv(t$p_auc), "p"),
      text_panel(fmt(t$dlrp, t$dlrp_lci, t$dlrp_uci), paste0("ΔLR+ (", conf_pct, ")")),
      text_panel(pv(t$p_lrp), "p"),
      text_panel(fmt(t$dlrn, t$dlrn_lci, t$dlrn_uci), paste0("ΔLR- (", conf_pct, ")")),
      text_panel(pv(t$p_lrn), "p"))
    widths <- c(label_w_in, 1.5, 1.6, 0.55, 1.5, 1.6, 0.55, 1.6, 0.55, 1.6, 0.55, 1.6, 0.55)
  }

  forest <- gridExtra::arrangeGrob(grobs = panels, ncol = length(panels),
                                   widths = widths,
                                   padding = grid::unit(0, "line"))
  forest_h <- grid::unit((ylim[2] - ylim[1]) * row_in + 0.35, "inches")

  grobs   <- list(forest)
  heights <- list(forest_h)
  if (!is.null(title)) {
    title_grob <- grid::textGrob(title, x = grid::unit(2, "pt"), hjust = 0,
                                 gp = grid::gpar(fontface = "bold", fontsize = 12))
    grobs   <- c(list(title_grob), grobs)
    heights <- c(list(grid::unit(1.8 * 12, "points")), heights)
  }

  # ----- sROC band ------------------------------------------------------------
  if (isTRUE(sroc)) {
    cols <- grDevices::hcl.colors(k, "Dark 3")
    if (loo$type == "single") {
      band <- .loo_sroc_panel(loo$geom, loo$geom_full, loo$points,
                              loo$letters, cols,
                              "sROC omitting each study")
      if (is.null(sroc_in)) sroc_in <- 4
    } else {
      ge <- lapply(loo$geom, function(g) g$e)
      gc <- lapply(loo$geom, function(g) g$c)
      p_e <- .loo_sroc_panel(ge, loo$geom_full$e, loo$points$e, loo$letters,
                             cols, paste0(loo$arms$e, ": sROC omitting each study"))
      p_c <- .loo_sroc_panel(gc, loo$geom_full$c, loo$points$c, loo$letters,
                             cols, paste0(loo$arms$c, ": sROC omitting each study"))
      band <- gridExtra::arrangeGrob(p_e, p_c, ncol = 2)
      if (is.null(sroc_in)) sroc_in <- 4
    }
    grobs   <- c(grobs, list(band))
    heights <- c(heights, list(grid::unit(sroc_in, "inches")))
  }

  total_h <- Reduce(`+`, heights)
  g <- gridExtra::arrangeGrob(grobs = grobs, ncol = 1,
                              heights = do.call(grid::unit.c, heights))
  attr(g, "height") <- total_h
  attr(g, "loo")    <- loo
  class(g) <- c("dta_loo_plot", class(g))
  g
}

#' @export
print.dta_loo_plot <- function(x, ...) {
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(
    x = grid::unit(0.5, "npc"), y = grid::unit(0.5, "npc"), just = "centre",
    height = attr(x, "height"), width = grid::unit(1, "npc")))
  grid::grid.draw(x)
  grid::popViewport()
  invisible(x)
}

# Unlabelled sROC: one thin coloured curve per omission (geom[[i]] = NULL
# when that refit failed), the full-data curve in black, study points drawn
# as their row letter in the colour of the curve that omits them.
.loo_sroc_panel <- function(geom, geom_full, pts, letters, cols, title,
                            n_grid = 200) {
  fpr <- seq(0.001, 0.999, length.out = n_grid)
  lsp <- -stats::qlogis(fpr)
  curve_of <- function(g) {
    if (is.null(g)) return(NULL)
    data.frame(fpr = fpr, tpr = stats::plogis(g$lsens - g$slope * (lsp - g$lspec)))
  }
  loo_curves <- do.call(rbind, lapply(seq_along(geom), function(i) {
    cv <- curve_of(geom[[i]])
    if (is.null(cv)) return(NULL)
    cv$letter <- letters[i]
    cv
  }))
  full <- curve_of(geom_full)
  pts$letter <- letters[seq_len(nrow(pts))]
  names(cols) <- letters

  p <- ggplot2::ggplot()
  if (!is.null(loo_curves))
    p <- p + ggplot2::geom_line(data = loo_curves,
                                ggplot2::aes(x = fpr, y = tpr, group = letter,
                                             colour = letter),
                                linewidth = 0.45, alpha = 0.85)
  if (!is.null(full))
    p <- p + ggplot2::geom_line(data = full, ggplot2::aes(x = fpr, y = tpr),
                                colour = "black", linewidth = 1) +
      ggplot2::annotate("point",
                        x = 1 - stats::plogis(geom_full$lspec),
                        y = stats::plogis(geom_full$lsens),
                        shape = 16, size = 3, colour = "black")
  p +
    ggplot2::geom_text(data = pts,
                       ggplot2::aes(x = fpr, y = tpr, label = letter,
                                    colour = letter),
                       size = 3, fontface = "bold") +
    ggplot2::scale_colour_manual(values = cols, guide = "none") +
    ggplot2::scale_x_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1),
                                expand = ggplot2::expansion(mult = 0.02)) +
    ggplot2::scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1),
                                expand = ggplot2::expansion(mult = 0.02)) +
    ggplot2::coord_fixed(ratio = 1, clip = "off") +
    ggplot2::labs(x = "False positive rate (1 - Specificity)",
                  y = "Sensitivity", title = title) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid      = ggplot2::element_blank(),
                   legend.position = "none",
                   plot.title      = ggplot2::element_text(face = "bold",
                                                           hjust = 0.5,
                                                           size = 10))
}
