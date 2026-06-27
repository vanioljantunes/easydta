# ============================================================================
# summary.R  -  Tidy summary + derived measures (DOR, LR+, LR-)
#
# Cochrane ref: Appendix 5, p.13-14 (Takwoingi et al., 2023). Implements the
# canonical delta-method recipe using msm::deltamethod on vcov(fit).
# ============================================================================

#' Tidy summary for a fitted easydta object
#'
#' @param object A `dta_single` or `dta_pairwise` object.
#' @param conf  Confidence level (default 0.95).
#' @param ...   Ignored.
#' @return A data frame of Se/Sp (and derived measures) with CIs.
#' @examples
#' data(anti_ccp2)
#' fit <- dta_fit_single(anti_ccp2, wide = TRUE)
#' dta_summary(fit)
#' @export
dta_summary <- function(object, conf = 0.95, ...) {
  UseMethod("dta_summary")
}

#' @export
dta_summary.dta_single <- function(object, conf = 0.95, ...) {
  f <- .fixed_se_sp(object$fit)
  se_row <- .logit_ci(f$lsens, f$se_lsens, conf)
  sp_row <- .logit_ci(f$lspec, f$se_lspec, conf)
  derived <- dta_derived(object, conf = conf)
  rbind(
    data.frame(measure = "Sensitivity",
               estimate = se_row["estimate"],
               lci = se_row["lci"], uci = se_row["uci"],
               row.names = NULL),
    data.frame(measure = "Specificity",
               estimate = sp_row["estimate"],
               lci = sp_row["lci"], uci = sp_row["uci"],
               row.names = NULL),
    derived
  )
}

#' @export
dta_summary.dta_pairwise <- function(object, conf = 0.95, ...) {
  fit <- object$models$B
  co  <- summary(fit)$coefficients
  V   <- stats::vcov(fit)
  lev <- object$levels
  z   <- stats::qnorm(1 - (1 - conf) / 2)

  build_row <- function(measure, coef_name, test_label) {
    est_logit <- co[coef_name, 1]
    se        <- co[coef_name, 2]
    data.frame(
      test     = test_label,
      measure  = measure,
      estimate = .inv_logit(est_logit),
      lci      = .inv_logit(est_logit - z * se),
      uci      = .inv_logit(est_logit + z * se),
      row.names = NULL
    )
  }

  rows <- rbind(
    build_row("Sensitivity", "seA", lev[1]),
    build_row("Specificity", "spA", lev[1]),
    build_row("Sensitivity", "seB", lev[2]),
    build_row("Specificity", "spB", lev[2])
  )

  derive_for <- function(se_name, sp_name, test_label) {
    lsens <- co[se_name, 1]
    lspec <- co[sp_name, 1]
    Vss <- as.matrix(V[c(se_name, sp_name), c(se_name, sp_name)])
    .derived_from_logits(lsens, lspec, Vss, conf, test_label)
  }
  d_rows <- rbind(
    derive_for("seA", "spA", lev[1]),
    derive_for("seB", "spB", lev[2])
  )

  rbind(rows, d_rows[, c("test", "measure", "estimate", "lci", "uci")])
}

#' Derived measures (DOR, LR+, LR-) with delta-method CIs
#'
#' Cochrane ref: Appendix 5, p.14.
#'
#' @param fit A `dta_single` object.
#' @param conf Confidence level (default 0.95).
#'
#' @return Data frame with one row per measure.
#' @examples
#' data(anti_ccp2)
#' fit <- dta_fit_single(anti_ccp2, wide = TRUE)
#' dta_derived(fit)
#' @export
dta_derived <- function(fit, conf = 0.95) {
  stopifnot(inherits(fit, "dta_single"))
  f <- .fixed_se_sp(fit$fit)
  .derived_from_logits(f$lsens, f$lspec, f$vcov_fixed, conf, test = NA)
}

.derived_from_logits <- function(lsens, lspec, Vss, conf, test = NA) {
  z <- stats::qnorm(1 - (1 - conf) / 2)

  DOR <- exp(lsens + lspec)
  LRp <- .inv_logit(lsens) / (1 - .inv_logit(lspec))
  LRn <- (1 - .inv_logit(lsens)) / .inv_logit(lspec)

  se.logDOR <- msm::deltamethod(
    ~ (x1 + x2), mean = c(lsens, lspec), cov = Vss)
  se.logLRp <- msm::deltamethod(
    ~ log((exp(x1) / (1 + exp(x1))) /
          (1 - exp(x2) / (1 + exp(x2)))),
    mean = c(lsens, lspec), cov = Vss)
  se.logLRn <- msm::deltamethod(
    ~ log((1 - exp(x1) / (1 + exp(x1))) /
          (exp(x2) / (1 + exp(x2)))),
    mean = c(lsens, lspec), cov = Vss)

  out <- data.frame(
    measure  = c("DOR", "LR+", "LR-"),
    estimate = c(DOR, LRp, LRn),
    lci = c(exp(log(DOR) - z * se.logDOR),
            exp(log(LRp) - z * se.logLRp),
            exp(log(LRn) - z * se.logLRn)),
    uci = c(exp(log(DOR) + z * se.logDOR),
            exp(log(LRp) + z * se.logLRp),
            exp(log(LRn) + z * se.logLRn)),
    row.names = NULL
  )
  if (!is.na(test)) out <- cbind(test = test, out)
  out
}

# -- S3 print / summary methods ---------------------------------------------

#' @export
print.dta_single <- function(x, digits = 3, ...) {
  cat("<dta_single>  Bivariate binomial GLMM (Cochrane Appendix 5)\n")
  cat("  Studies: ", length(unique(x$long$studlab)),
      "   Observations: ", nrow(x$long), "\n", sep = "")
  cat("  nAGQ = ", x$call_args$nAGQ, "\n\n", sep = "")
  s <- dta_summary(x)
  print(format(s, digits = digits), row.names = FALSE)
  if (!is.null(x$heterogeneity)) {
    cat("\n")
    .print_heterogeneity_block(x$heterogeneity, digits = digits)
  }
  invisible(x)
}

#' @export
summary.dta_single <- function(object, ...) dta_summary(object, ...)

#' @export
print.dta_pairwise <- function(x, digits = 3, ...) {
  cat("<dta_pairwise>  Bivariate meta-regression (Cochrane Appendix 12)\n")
  cat("  Test variable: ", x$test_var, "\n", sep = "")
  cat("  Levels:        ", paste(x$levels, collapse = " vs "), "\n", sep = "")
  var_lab <- if (identical(x$variance, "unequal")) {
    "unequal (model E)"
  } else {
    "equal (model B)"
  }
  cat("  Variance:      ", var_lab, "\n", sep = "")
  cat("  Studies: ", length(unique(x$long$studlab)),
      "   Observations: ", nrow(x$long), "\n\n", sep = "")
  s <- dta_summary(x)
  print(format(s, digits = digits), row.names = FALSE)
  invisible(x)
}

#' @export
summary.dta_pairwise <- function(object, ...) dta_summary(object, ...)
