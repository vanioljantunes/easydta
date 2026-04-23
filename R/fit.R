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
#' @param long   Long-format data frame produced by `dta_reshape()`, OR a
#'   wide data frame -- in which case `wide = TRUE` will call `dta_reshape()`.
#' @param nAGQ   Integer. 1 = Laplace (default, fastest, matches the Handbook
#'   baseline). Increase (e.g. 5, 7) for adaptive Gauss-Hermite quadrature.
#' @param wide   Logical. If TRUE, `long` is actually wide and will be
#'   reshaped first.
#' @param ...    Forwarded to `dta_reshape()` when `wide = TRUE`.
#'
#' @return An S3 object of class `"dta_single"` (a plain list).
#'
#' @export
dta_fit_single <- function(long, nAGQ = 1L, wide = FALSE, ...) {
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
    call_args  = list(nAGQ = nAGQ),
    vcov_fixed = stats::vcov(fit),
    Psi        = .random_vcv(fit)
  )
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
#' with equal variances: (0 + sens + spec | study_id) in every model.
#' The user-supplied `test_var` is coerced to a two-level factor; the first
#' level (alphabetical) becomes "A", the second "B".
#'
#' @param long      Long-format data frame from `dta_reshape()` (with
#'   `extra` including the test-type column) OR a wide data frame when
#'   `wide = TRUE`.
#' @param test_var  Name (string) of the test-type column -- must have
#'   exactly two distinct non-NA values.
#' @param nAGQ      Integer; Laplace by default.
#' @param wide      Logical. If TRUE, reshape first.
#' @param ...       Forwarded to `dta_reshape()` when `wide = TRUE`.
#'
#' @return An S3 object of class `"dta_pairwise"`.
#'
#' @export
dta_fit_pairwise <- function(long,
                             test_var,
                             nAGQ = 1L,
                             wide = FALSE,
                             ...) {
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

  f_null <- cbind(true, n - true) ~ 0 + sens + spec +
    (0 + sens + spec | studlab)
  f_full <- cbind(true, n - true) ~ 0 + seA + seB + spA + spB +
    (0 + sens + spec | studlab)
  f_spOnly <- cbind(true, n - true) ~ 0 + sens + spA + spB +
    (0 + sens + spec | studlab)
  f_seOnly <- cbind(true, n - true) ~ 0 + seA + seB + spec +
    (0 + sens + spec | studlab)

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
    levels     = levs,
    call_args  = list(nAGQ = nAGQ),
    vcov_full  = stats::vcov(B),
    Psi_full   = .random_vcv(B)
  )
  class(out) <- "dta_pairwise"
  out
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
