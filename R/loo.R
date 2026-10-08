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
# dta_loo_plot() draws, like dta_sroc_pair(), the sROC panel(s) on top
# (one curve per omission; each study point drawn as its row letter in the
# colour of the curve that omits it) and a bordered summary table beneath,
# one letter-coded row per omitted study plus the pooled estimate.
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
#'   the plot) and the sROC geometry of every refit.  The single-test table
#'   also carries `i2_sens`, `i2_spec` and `i2_biv` (Zhou-Dendukuri bivariate
#'   I^2 of each refit); the pairwise table carries `i2_e` and `i2_c`.
#'   `print()` shows the table; `plot()` / [dta_loo_plot()] draws the
#'   sROC panels and table figure.
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
    lrn  = na, lrn_lci  = na, lrn_uci  = na,
    i2_sens = na, i2_spec = na, i2_biv = na)
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
  h   <- fit$heterogeneity

  row[] <- list(
    se[["estimate"]], se[["lci"]], se[["uci"]],
    sp[["estimate"]], sp[["lci"]], sp[["uci"]],
    g$AUC, aci[1], aci[2],
    lrp$estimate, lrp$lci, lrp$uci,
    lrn$estimate, lrn$lci, lrn$uci,
    .loo_i2(h, "I2_sens"), .loo_i2(h, "I2_spec"), .loo_i2(h, "I2_biv"))
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
    dlrn  = na, dlrn_lci  = na, dlrn_uci  = na, p_lrn  = na,
    i2_e = na, i2_c = na)
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
    dlrn, dlrn - z * se_dlrn, dlrn + z * se_dlrn, wald_p(dlrn, se_dlrn),
    .loo_i2(fit_e$heterogeneity, "I2_biv"),
    .loo_i2(fit_c$heterogeneity, "I2_biv"))

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

# Zhou-Dendukuri bivariate I^2 of a fit's heterogeneity block, as a
# proportion; NA when the refit failed or the statistic is undefined.
.loo_i2 <- function(h, what) {
  if (is.null(h) || is.null(h[[what]])) return(NA_real_)
  as.numeric(h[[what]])
}

.loo_fmt_i2 <- function(p) {
  ifelse(is.na(p), "n/a", sprintf("%.0f%%", 100 * p))
}

# Display table: one text column per measure (and per p-value).
.loo_display <- function(x, digits = 2, lr.show = FALSE, i2.show = TRUE) {
  t <- x$table
  lab <- ifelse(t$pooled, t$studlab, paste("Omitting", t$studlab))
  out <- if (x$type == "single") {
    data.frame(
      Study  = lab,
      Sens   = .loo_fmt_ci(t$sens, t$sens_lci, t$sens_uci, digits),
      Spec   = .loo_fmt_ci(t$spec, t$spec_lci, t$spec_uci, digits),
      AUC    = .loo_fmt_ci(t$auc,  t$auc_lci,  t$auc_uci,  digits),
      `LR+`  = .loo_fmt_ci(t$lrp,  t$lrp_lci,  t$lrp_uci,  digits),
      `LR-`  = .loo_fmt_ci(t$lrn,  t$lrn_lci,  t$lrn_uci,  digits),
      `I2`   = .loo_fmt_i2(t$i2_biv),
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
      `I2.e`   = .loo_fmt_i2(t$i2_e),
      `I2.c`   = .loo_fmt_i2(t$i2_c),
      check.names = FALSE, stringsAsFactors = FALSE)
  }
  # LR+/LR- are always computed (x$table); shown only when lr.show = TRUE
  if (!isTRUE(lr.show)) {
    keep <- if (x$type == "single") !names(out) %in% c("LR+", "LR-")
            else !names(out) %in% c("dLR+", "p   ", "dLR-", "p    ")
    out <- out[, keep, drop = FALSE]
  }
  # I^2 is on by default; i2.show = FALSE drops the column(s)
  if (!isTRUE(i2.show)) {
    keep <- if (x$type == "single") names(out) != "I2"
            else !names(out) %in% c("I2.e", "I2.c")
    out <- out[, keep, drop = FALSE]
  } else if (x$type == "pair") {
    names(out)[names(out) == "I2.e"] <- paste0("I2 ", x$arms$e)
    names(out)[names(out) == "I2.c"] <- paste0("I2 ", x$arms$c)
  }
  out
}

#' @export
print.dta_loo <- function(x, digits = 2, lr.show = FALSE, i2.show = TRUE,
                          ...) {
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
  print(.loo_display(x, digits, lr.show, i2.show), row.names = FALSE,
        right = FALSE)
  invisible(x)
}

#' @export
plot.dta_loo <- function(x, ...) dta_loo_plot(x, ...)

# -- sROC panels + summary table figure ------------------------------------

#' Leave-one-out sROC panels with summary table
#'
#' Draws the leave-one-out analysis of [dta_loo()] in the layout of
#' [dta_sroc_pair()]: the sROC panels on top (one per arm for a comparison,
#' a single panel for one test) and a bordered summary table beneath with
#' one row per omitted study ("Omitting ...") and the pooled estimate last.
#' For a single test the table columns are Sens, Spec, AUC, LR+ and LR-
#' (each with CI) and the bivariate I^2; for a comparison each column is the
#' difference `.e - .c` followed by its p-value, then the I^2 of each arm.
#'
#' Each sROC panel is unlabelled: one thin coloured curve per omission and
#' the full-data curve in black.  Every study point is drawn as its row
#' letter, in the colour of the curve fitted without that study, so an
#' influential study is read off directly: its letter sits away from the
#' cloud and its curve departs from the black one.
#'
#' @param x       A `dta_loo` object, or a `dta_single` /
#'   `dta_pairwise_result` (then [dta_loo()] is run first; `...` is
#'   forwarded to it).
#' @param table   Logical. Show the summary table? Default `TRUE`.
#' @param table.position Where the table sits relative to the sROC
#'   panel(s): `"below"` (default), `"right"`, `"left"` or `"above"`.
#' @param digits  Display digits for the table (default 2).
#' @param lr.show Logical. Show the LR+ / LR- columns? They are always
#'   computed (see `x$table`); default `FALSE` keeps the table compact.
#' @param i2.show Logical. Show the Zhou-Dendukuri bivariate I^2 column
#'   (one per arm for a comparison)? Default `TRUE`.
#' @param ...     Forwarded to [dta_loo()] when `x` is not a `dta_loo`.
#'
#' @return A `gtable` of class `"dta_loo_plot"`; its `print()` method draws
#'   it, so it renders when returned at top level.  `attr(., "panels")`
#'   holds the ggplot panel(s) and `attr(., "table")` the rendered
#'   data.frame.
#' @examples
#' \donttest{
#' data(anti_ccp2)
#' fit <- dta_fit_single(anti_ccp2, wide = TRUE)
#' dta_loo_plot(fit, auc_ci = FALSE)
#' }
#' @export
dta_loo_plot <- function(x, table = TRUE,
                         table.position = c("below", "right", "left", "above"),
                         digits = 2, lr.show = FALSE, i2.show = TRUE, ...) {
  # tolerate the informal spellings "bellow" / "up"
  table.position <- switch(as.character(table.position[1]),
                           bellow = "below", up = "above", table.position[1])
  table.position <- match.arg(table.position,
                              c("below", "right", "left", "above"))
  if (!inherits(x, "dta_loo")) x <- dta_loo(x, ...)
  loo  <- x
  k    <- length(loo$studies)
  cols <- grDevices::hcl.colors(k, "Dark 3")

  if (loo$type == "single") {
    p <- .loo_sroc_panel(loo$geom, loo$geom_full, loo$points, loo$letters,
                         cols, "sROC omitting each study")
    panels_list <- list(p)
    panels <- gridExtra::arrangeGrob(p, ncol = 1)
  } else {
    ge  <- lapply(loo$geom, function(g) g$e)
    gc  <- lapply(loo$geom, function(g) g$c)
    p_e <- .loo_sroc_panel(ge, loo$geom_full$e, loo$points$e, loo$letters,
                           cols, paste0("sROC of ", loo$arms$e,
                                        "\nomitting each study"))
    p_c <- .loo_sroc_panel(gc, loo$geom_full$c, loo$points$c, loo$letters,
                           cols, paste0("sROC of ", loo$arms$c,
                                        "\nomitting each study"))
    panels_list <- list(.e = p_e, .c = p_c)
    panels <- gridExtra::arrangeGrob(p_e, p_c, ncol = 2)
  }

  tbl_df <- cbind(data.frame(" " = loo$table$letter, check.names = FALSE),
                  .loo_display(loo, digits, lr.show, i2.show))
  # console print keeps ASCII "d"; the drawn table uses the delta sign
  names(tbl_df) <- sub("^d(Sens|Spec|AUC|LR)", "Δ\\1", names(tbl_df))
  # and the I2 column(s) are drawn with the superscript
  names(tbl_df) <- sub("^I2", "I²", names(tbl_df))
  if (isTRUE(table)) {
    face <- ifelse(loo$table$pooled, "bold", "plain")
    tbl_grob <- gridExtra::tableGrob(
      tbl_df, rows = NULL,
      theme = gridExtra::ttheme_minimal(
        core    = list(fg_params = list(cex = 0.8,
                                        fontface = rep(face, ncol(tbl_df)))),
        colhead = list(fg_params = list(cex = 0.8, fontface = "bold"))))
    tbl_grob <- gtable::gtable_add_grob(
      tbl_grob,
      grobs = grid::rectGrob(gp = grid::gpar(fill = NA, col = "black",
                                             lwd = 1.2)),
      t = 1, b = nrow(tbl_grob), l = 1, r = ncol(tbl_grob),
      name = "container-border")
    # grobHeight() undercounts a gtable; sum the row heights instead so the
    # band always fits the whole table (plus a small margin).
    tbl_h <- sum(tbl_grob$heights) + grid::unit(40, "pt")
    tbl_w <- sum(tbl_grob$widths)  + grid::unit(40, "pt")
    rest_h <- grid::unit(1, "npc") - tbl_h
    rest_w <- grid::unit(1, "npc") - tbl_w
    g <- switch(table.position,
      below = gridExtra::arrangeGrob(panels, tbl_grob, nrow = 2,
                                     heights = grid::unit.c(rest_h, tbl_h)),
      above = gridExtra::arrangeGrob(tbl_grob, panels, nrow = 2,
                                     heights = grid::unit.c(tbl_h, rest_h)),
      right = gridExtra::arrangeGrob(panels, tbl_grob, ncol = 2,
                                     widths = grid::unit.c(rest_w, tbl_w)),
      left  = gridExtra::arrangeGrob(tbl_grob, panels, ncol = 2,
                                     widths = grid::unit.c(tbl_w, rest_w)))
  } else {
    g <- panels
  }

  attr(g, "panels") <- panels_list
  attr(g, "table")  <- tbl_df
  attr(g, "loo")    <- loo
  class(g) <- c("dta_loo_plot", class(g))
  g
}

#' @export
print.dta_loo_plot <- function(x, ...) {
  grid::grid.newpage()
  grid::grid.draw(x)
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
