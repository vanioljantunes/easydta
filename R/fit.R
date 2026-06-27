# ============================================================================
# fit.R  -  Bivariate binomial GLMM fitters
#
# Cochrane ref:
#   Appendix 5, p.12  -- single test
#   Appendix 12, p.55-56 -- pairwise with covariate
#   Appendix 14, p.62-64 -- direct comparison variant
#
# Both fitters wrap lme4::glmer and return an S3 object carrying:
#   $fit         -- the glmerMod (or list of glmerMods for pairwise)
#   $long        -- the long data actually fit
#   $call_args   -- args passed by the user (for reproducibility)
#   $vcov_fixed  -- fixed-effect VCV (cached to avoid recomputation)
#   $Psi         -- between-study VCV (2x2 random-effect covariance)
# ============================================================================

#' Fit the bivariate binomial GLMM for a single diagnostic test
#'
#' Implements Cochrane Handbook Appendix 5 (Takwoingi et al., 2023):
#'
#'   glmer(cbind(true, n - true) ~ 0 + sens + spec +
#'             (0 + sens + spec | study_id),
#'         data = long, family = binomial, nAGQ = 1)
#'
#' @param long   Data frame: either long format (from `dta_reshape()`, with
#'   `sens`/`spec`/`true`/`n` columns) OR a wide one-row-per-study 2x2 table
#'   (with `TP`/`FP`/`FN`/`TN` columns). The format is auto-detected; pass
#'   `wide` explicitly to override.
#' @param nAGQ   Integer. 1 = Laplace (default, fastest, matches the Handbook
#'   baseline). Increase (e.g. 5, 7) for adaptive Gauss-Hermite quadrature.
#' @param wide   Logical, or `NA` (default) to auto-detect: a frame lacking the
#'   long-format columns `sens`/`spec`/`true`/`n` is treated as wide and
#'   reshaped via `dta_reshape()` first.
#' @param conf   Confidence level used for the prediction-region ellipse
#'   computed as part of the heterogeneity summary (default 0.95).
#' @param ...    Forwarded to `dta_reshape()` when `wide = TRUE`.
#'
#' @return An S3 object of class `"dta_single"` (a plain list). The object
#'   carries `$heterogeneity`, a list with `tau_sens`, `tau_spec`, `rho`,
#'   `I2_sens`, `I2_spec`, `I2_biv` (Zhou-Dendukuri bivariate I^2), and
#'   `pred_region` (logit-scale ellipse). `print(fit)` shows it.
#'
#' @examples
#' data(anti_ccp2)
#' fit <- dta_fit_single(anti_ccp2, wide = TRUE)
#' print(fit)
#' @export
dta_fit_single <- function(long, nAGQ = 1L, wide = NA, conf = 0.95, ...) {
  if (is.na(wide)) {
    wide <- !all(c("sens", "spec", "true", "n") %in% names(long))
  }
  if (wide) long <- dta_reshape(long, ...)
  .check_long(long)

  fit <- lme4::glmer(
    formula = cbind(true, n - true) ~ 0 + sens + spec +
              (0 + sens + spec | studlab),
    data    = long,
    family  = stats::binomial,
    nAGQ    = nAGQ
  )

  out <- list(
    fit        = fit,
    long       = long,
    call_args  = list(nAGQ = nAGQ, conf = conf),
    vcov_fixed = stats::vcov(fit),
    Psi        = .random_vcv(fit)
  )
  out$heterogeneity <- .compute_heterogeneity(out, conf = conf)
  class(out) <- "dta_single"
  out
}

#' Fit the bivariate meta-regression for pairwise test comparison
#'
#' Fits four nested models in the order recommended by Cochrane Handbook
#' Appendix 12 (p.56-57):
#'
#'   A. null:                0 + sens + spec
#'   B. both-vary (full):    0 + seA + seB + spA + spB
#'   C. spec-varies only:    0 + sens + spA + spB
#'   D. sens-varies only:    0 + seA + seB + spec
#'
#' The random-effects structure is controlled by `variance`:
#'   * `"equal"`   (default, Cochrane model B): a single shared block
#'     `(0 + sens + spec | study_id)` in every model.
#'   * `"unequal"` (Cochrane model E): two per-test blocks
#'     `(0 + seA + spA | study_id) + (0 + seB + spB | study_id)`, allowing
#'     each test its own between-study variances and covariance.
#'
#' The chosen structure is applied to all four nested models so the
#' likelihood-ratio ladder stays nested.
#'
#' The user-supplied `test_var` is coerced to a two-level factor; the first
#' level (alphabetical) becomes "A", the second "B".
#'
#' @param long      Long-format data frame from `dta_reshape()` (with
#'   `extra` including the test-type column) OR a wide data frame when
#'   `wide = TRUE`.
#' @param test_var  Name (string) of the test-type column -- must have
#'   exactly two distinct non-NA values.
#' @param variance  Random-effects structure: `"equal"` (default, model B,
#'   shared between-study variances) or `"unequal"` (model E, separate
#'   per-test variances).
#' @param nAGQ      Integer; Laplace by default.
#' @param wide      Logical. If TRUE, reshape first.
#' @param ...       Forwarded to `dta_reshape()` when `wide = TRUE`.
#'
#' @return An S3 object of class `"dta_pairwise"`.
#'
#' @examples
#' data(schuetz)
#' long <- dta_reshape_pairwise(schuetz, studlab = "studlab",
#'                              intervention = "CT", control = "MRI")
#' pair <- dta_fit_pairwise(long, test_var = "test")
#' @export
dta_fit_pairwise <- function(long,
                             test_var,
                             variance = c("equal", "unequal"),
                             nAGQ = 1L,
                             wide = FALSE,
                             ...) {
  variance <- match.arg(variance)
  if (wide) long <- dta_reshape(long, extra = test_var, ...)
  .check_long(long)
  if (!test_var %in% names(long)) {
    stop("test_var '", test_var, "' not found in data.")
  }

  tv <- long[[test_var]]
  levs <- sort(unique(stats::na.omit(as.character(tv))))
  if (length(levs) != 2L) {
    stop("test_var '", test_var, "' must have exactly two levels; found ",
         length(levs), ".")
  }

  long$seA <- long$sens * as.integer(tv == levs[1])
  long$seB <- long$sens * as.integer(tv == levs[2])
  long$spA <- long$spec * as.integer(tv == levs[1])
  long$spB <- long$spec * as.integer(tv == levs[2])

  # Random-effects term: shared (model B) or split per test (model E).
  re_term <- if (variance == "equal") {
    "(0 + sens + spec | studlab)"
  } else {
    "(0 + seA + spA | studlab) + (0 + seB + spB | studlab)"
  }

  mk_formula <- function(fixed) {
    stats::as.formula(
      paste("cbind(true, n - true) ~", fixed, "+", re_term)
    )
  }

  f_null   <- mk_formula("0 + sens + spec")
  f_full   <- mk_formula("0 + seA + seB + spA + spB")
  f_spOnly <- mk_formula("0 + sens + spA + spB")
  f_seOnly <- mk_formula("0 + seA + seB + spec")

  mk <- function(formula) {
    lme4::glmer(formula = formula, data = long,
                family = stats::binomial, nAGQ = nAGQ)
  }

  A <- mk(f_null)
  B <- mk(f_full)
  C <- mk(f_spOnly)
  D <- mk(f_seOnly)

  out <- list(
    models     = list(A = A, B = B, C = C, D = D),
    long       = long,
    test_var   = test_var,
    variance   = variance,
    levels     = levs,
    call_args  = list(nAGQ = nAGQ, variance = variance),
    vcov_full  = stats::vcov(B),
    Psi_full   = .random_vcv(B, levels = levs)
  )
  class(out) <- "dta_pairwise"
  out
}

#' One-call pairwise DTA meta-analysis (paired wide design)
#'
#' Takes a single wide data frame -- one row per study, with `.e`
#' (intervention / index) and `.c` (control / comparator) suffixed count
#' columns -- as used for paired head-to-head designs where each study
#' evaluated BOTH tests (e.g. the `schuetz` CT vs MRI data).  Runs the full
#' pipeline:
#'   1. `dta_reshape_pairwise()` -> long format (two arms stacked, shared
#'      `studlab` so the within-study pairing is captured by the random effect)
#'   2. `dta_fit_pairwise()`     -> nested LR-test models
#'   3. `dta_compare()`          -> LR tests + Se/Sp differences
#'   4. `dta_fit_single()` per arm, keyed by the `intervention` / `control`
#'      labels, ready for `dta_forest()` and `dta_sroc()`.
#'
#' @inheritParams dta_reshape_pairwise
#' @param variance  Random-effects structure for the comparison: `"equal"`
#'   (default, model B) or `"unequal"` (model E). See [dta_fit_pairwise()].
#' @param nAGQ  Integer; Laplace by default (passed to both fitters).
#' @param conf  Confidence level used by per-arm heterogeneity summaries
#'   and by `dta_compare()` (default 0.95).
#'
#' @return An S3 object of class `"dta_pairwise_result"` with components:
#'   `long`, `pair` (the `dta_pairwise` fit), `compare` (the `dta_compare`
#'   result), `arms` (named list of two `dta_single` fits), and `labels`
#'   (the intervention / control labels).
#'
#' @examples
#' data(schuetz)
#' res <- dta_pairwise(schuetz, studlab = "studlab",
#'                     intervention = "CT", control = "MRI")
#' print(res)
#' @export
dta_pairwise <- function(data,
                         author       = "author",
                         year         = "year",
                         intervention = "Intervention",
                         control      = "Control",
                         tp.e         = "TP.e",
                         fp.e         = "FP.e",
                         fn.e         = "FN.e",
                         tn.e         = "TN.e",
                         tp.c         = "TP.c",
                         fp.c         = "FP.c",
                         fn.c         = "FN.c",
                         tn.c         = "TN.c",
                         studlab      = NULL,
                         test_var     = "test",
                         variance     = c("equal", "unequal"),
                         nAGQ         = 1L,
                         conf         = 0.95) {
  variance <- match.arg(variance)
  long <- dta_reshape_pairwise(
    data,
    author       = author,
    year         = year,
    intervention = intervention,
    control      = control,
    tp.e = tp.e, fp.e = fp.e, fn.e = fn.e, tn.e = tn.e,
    tp.c = tp.c, fp.c = fp.c, fn.c = fn.c, tn.c = tn.c,
    studlab      = studlab,
    test_var     = test_var
  )

  pair <- dta_fit_pairwise(long, test_var = test_var,
                           variance = variance, nAGQ = nAGQ)
  cmp  <- dta_compare(pair, conf = conf)

  arm_rows <- function(label) long[long[[test_var]] == label, , drop = FALSE]
  fit_e <- dta_fit_single(arm_rows(intervention), nAGQ = nAGQ, conf = conf)
  fit_c <- dta_fit_single(arm_rows(control),      nAGQ = nAGQ, conf = conf)

  arms <- stats::setNames(list(fit_e, fit_c), c(intervention, control))

  out <- list(
    long     = long,
    pair     = pair,
    compare  = cmp,
    arms     = arms,
    labels   = list(intervention = intervention, control = control),
    test_var = test_var
  )
  class(out) <- "dta_pairwise_result"
  out
}

#' One-call test-comparison meta-analysis (covariate / unpaired design)
#'
#' Sister of [dta_pairwise()] for the case where two diagnostic tests are
#' compared via a between-study covariate (Cochrane Handbook chapter 10 /
#' Appendix 14) rather than a paired within-study design.  Each row of
#' `data` is one study with a single 2x2 table; the `test_var` column
#' carries the test type and must have exactly two levels.  The pipeline
#' is the same as `dta_pairwise()`:
#'   1. `dta_reshape()` (carrying `test_var` as `extra`)
#'   2. `dta_fit_pairwise()` -> nested LR-test models
#'   3. `dta_compare()`      -> LR tests + Se/Sp differences
#'   4. `dta_fit_single()` per arm, keyed by the two `test_var` levels.
#'
#' @param data      Data frame, one row per study, with a 2x2 table and a
#'   two-level test-type column.
#' @param test_var  Name (string) of the test-type column.
#' @param studlab   Column name of the study label (default `"studlab"`).
#' @param tp,fp,fn,tn  Column names for the 2x2 counts.
#' @param variance  Random-effects structure for the comparison: `"equal"`
#'   (default, model B) or `"unequal"` (model E). See [dta_fit_pairwise()].
#' @param nAGQ      Integer; Laplace by default.
#' @param conf      Confidence level (default 0.95).
#'
#' @return An S3 object of class `"dta_pairwise_result"` (same shape as
#'   the return of `dta_pairwise()`).
#'
#' @examples
#' data(anti_ccp1); data(anti_ccp2)
#' d <- rbind(anti_ccp1, anti_ccp2)          # between-study covariate design
#' res <- dta_compare_tests(d, test_var = "test")
#' print(res)
#' @export
dta_compare_tests <- function(data,
                              test_var,
                              studlab = "studlab",
                              tp = "TP", fp = "FP", fn = "FN", tn = "TN",
                              variance = c("equal", "unequal"),
                              nAGQ = 1L,
                              conf = 0.95) {
  variance <- match.arg(variance)
  stopifnot(is.data.frame(data))
  if (missing(test_var) || is.null(test_var)) {
    stop("`test_var` is required (the column carrying the test type).")
  }
  if (!test_var %in% names(data)) {
    stop("test_var '", test_var, "' not found in data.")
  }

  long <- dta_reshape(data,
                      tp = tp, fp = fp, fn = fn, tn = tn,
                      studlab = studlab,
                      extra   = test_var)

  pair <- dta_fit_pairwise(long, test_var = test_var,
                           variance = variance, nAGQ = nAGQ)
  cmp  <- dta_compare(pair, conf = conf)

  levs <- pair$levels
  arm_rows <- function(label) long[long[[test_var]] == label, , drop = FALSE]
  fits <- lapply(levs, function(lv) {
    dta_fit_single(arm_rows(lv), nAGQ = nAGQ, conf = conf)
  })
  arms <- stats::setNames(fits, levs)

  out <- list(
    long     = long,
    pair     = pair,
    compare  = cmp,
    arms     = arms,
    labels   = list(intervention = levs[1], control = levs[2]),
    test_var = test_var
  )
  class(out) <- "dta_pairwise_result"
  out
}

#' @export
print.dta_pairwise_result <- function(x, ...) {
  cat("<dta_pairwise_result>  One-call pairwise DTA analysis\n")
  cat("  Arms: ", x$labels$intervention, " vs ", x$labels$control, "\n\n",
      sep = "")
  print(x$pair)
  cat("\n")
  print(x$compare)
  invisible(x)
}

# -- Internals --------------------------------------------------------------

.check_long <- function(long) {
  need <- c("studlab", "sens", "spec", "true", "n")
  miss <- setdiff(need, names(long))
  if (length(miss)) {
    stop("Long data is missing columns: ", paste(miss, collapse = ", "),
         ". Did you run dta_reshape() first?")
  }
}
