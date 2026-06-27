# ============================================================================
# compare.R  -  Pairwise comparison: LR tests + delta-method differences
#
# Cochrane ref: Appendix 12, p.56-58 (Takwoingi et al., 2023).
#   Likelihood-ratio tests: A vs B (overall test effect), B vs C (Se differs?),
#   B vs D (Sp differs?).  Absolute and relative differences in Se and Sp
#   between the two tests computed from model B via msm::deltamethod.
# ============================================================================

#' Likelihood-ratio tests and Se/Sp differences between two tests
#'
#' @param object  A `dta_pairwise` object from [dta_fit_pairwise()].
#' @param conf    Confidence level for difference CIs (default 0.95).
#'
#' @return An S3 `dta_compare` object.
#' @examples
#' data(schuetz)
#' long <- dta_reshape_pairwise(schuetz, studlab = "studlab",
#'                              intervention = "CT", control = "MRI")
#' pair <- dta_fit_pairwise(long, test_var = "test")
#' dta_compare(pair)
#' @export
dta_compare <- function(object, conf = 0.95) {
  stopifnot(inherits(object, "dta_pairwise"))
  M <- object$models
  lev <- object$levels
  z <- stats::qnorm(1 - (1 - conf) / 2)

  lr_tbl <- rbind(
    .lrtest_row(M$A, M$B, "Null (A) vs Full (B) -- overall test effect"),
    .lrtest_row(M$C, M$B, "Se-common (C) vs Full (B) -- does Se differ?"),
    .lrtest_row(M$D, M$B, "Sp-common (D) vs Full (B) -- does Sp differ?")
  )

  co <- summary(M$B)$coefficients
  parts <- .split_pairwise_vcov(M$B)
  seA <- co["seA", 1]; seB <- co["seB", 1]
  spA <- co["spA", 1]; spB <- co["spB", 1]

  d_se <- .inv_logit(seA) - .inv_logit(seB)
  se.d_se <- msm::deltamethod(
    ~ (exp(x1)/(1+exp(x1))) - (exp(x2)/(1+exp(x2))),
    mean = c(seA, seB), cov = parts$se_cov)
  d_sp <- .inv_logit(spA) - .inv_logit(spB)
  se.d_sp <- msm::deltamethod(
    ~ (exp(x1)/(1+exp(x1))) - (exp(x2)/(1+exp(x2))),
    mean = c(spA, spB), cov = parts$sp_cov)
  r_se <- .inv_logit(seA) / .inv_logit(seB)
  se.log_rse <- msm::deltamethod(
    ~ log((exp(x1)/(1+exp(x1))) / (exp(x2)/(1+exp(x2)))),
    mean = c(seA, seB), cov = parts$se_cov)
  r_sp <- .inv_logit(spA) / .inv_logit(spB)
  se.log_rsp <- msm::deltamethod(
    ~ log((exp(x1)/(1+exp(x1))) / (exp(x2)/(1+exp(x2)))),
    mean = c(spA, spB), cov = parts$sp_cov)

  diff_tbl <- data.frame(
    comparison = paste0(lev[1], " vs ", lev[2]),
    measure = c("Absolute diff Sens", "Absolute diff Spec",
                "Relative Sens",      "Relative Spec"),
    estimate = c(d_se, d_sp, r_se, r_sp),
    lci = c(d_se - z * se.d_se,
            d_sp - z * se.d_sp,
            exp(log(r_se) - z * se.log_rse),
            exp(log(r_sp) - z * se.log_rsp)),
    uci = c(d_se + z * se.d_se,
            d_sp + z * se.d_sp,
            exp(log(r_se) + z * se.log_rse),
            exp(log(r_sp) + z * se.log_rsp)),
    row.names = NULL
  )

  out <- list(lr_tests = lr_tbl, differences = diff_tbl, levels = lev)
  class(out) <- "dta_compare"
  out
}

.lrtest_row <- function(m_small, m_big, label) {
  lr <- lmtest::lrtest(m_small, m_big)
  data.frame(
    comparison = label,
    chisq      = lr$Chisq[2],
    df         = lr$Df[2],
    p_value    = lr[["Pr(>Chisq)"]][2],
    row.names = NULL
  )
}

#' @export
print.dta_compare <- function(x, digits = 3, ...) {
  cat("<dta_compare>  Pairwise test comparison (Cochrane Appendix 12)\n")
  cat("  Levels: ", paste(x$levels, collapse = " vs "), "\n\n", sep = "")
  cat("Likelihood-ratio tests:\n")
  print(format(x$lr_tests, digits = digits), row.names = FALSE)
  cat("\nDifferences (Test A vs Test B):\n")
  print(format(x$differences, digits = digits), row.names = FALSE)
  invisible(x)
}
