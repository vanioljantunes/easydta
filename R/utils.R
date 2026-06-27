# ============================================================================
# utils.R  -  Shared helpers for the easydta package
#
# Cochrane ref: Chapter 10 Supplementary Material 1 (Takwoingi et al., 2023)
#               Appendix 5, Appendix 12. Helpers are adapted from those
#               appendices to avoid duplication across the other files.
# ============================================================================

# -- Logit helpers -----------------------------------------------------------

.logit <- function(p) stats::qlogis(p)
.inv_logit <- function(x) stats::plogis(x)

# Build a Wald CI on the logit scale and back-transform to a probability.
# Used for Se, Sp summary rows.
.logit_ci <- function(logit_est, se, conf = 0.95) {
  z <- stats::qnorm(1 - (1 - conf) / 2)
  est <- .inv_logit(logit_est)
  lo  <- .inv_logit(logit_est - z * se)
  hi  <- .inv_logit(logit_est + z * se)
  c(estimate = est, lci = lo, uci = hi)
}

# Exact (Clopper-Pearson) binomial CI for a single study. Used by forest
# plots for the per-study rows -- matches the RevMan/Cochrane display.
.exact_binom_ci <- function(x, n, conf = 0.95) {
  if (is.na(x) || is.na(n) || n == 0) {
    return(c(estimate = NA_real_, lci = NA_real_, uci = NA_real_))
  }
  bt <- stats::binom.test(x, n, conf.level = conf)
  c(estimate = unname(bt$estimate),
    lci      = bt$conf.int[1],
    uci      = bt$conf.int[2])
}

# -- Ellipse on the logit ROC plane -----------------------------------------
#
# Given a 2x2 covariance matrix V (in the (logit Se, logit Sp) coordinate
# system) and a centre (logitSe, logitSp), return a data frame with columns
# fpr, tpr describing the ellipse boundary on probability scale.
#
# Used twice by dta_sroc():
#   (a) 95% confidence region -> V = vcov(fit)               [fixed only]
#   (b) 95% prediction region -> V = vcov(fit) + Psi_random  [+ between-study]
.logit_ellipse <- function(centre, V, conf = 0.95, n_pts = 200) {
  # chi-sq 2 df critical value for the requested confidence level
  k <- stats::qchisq(conf, df = 2)
  ev <- eigen(V, symmetric = TRUE)
  theta <- seq(0, 2 * pi, length.out = n_pts)
  circle <- rbind(cos(theta), sin(theta))
  ell <- ev$vectors %*% diag(sqrt(ev$values * k)) %*% circle
  logitSe <- centre[1] + ell[1, ]
  logitSp <- centre[2] + ell[2, ]
  data.frame(fpr = 1 - .inv_logit(logitSp),
             tpr =      .inv_logit(logitSe))
}

# -- VCV partitioning for pairwise comparisons ------------------------------
#
# Given a glmer fit with fixed effects ordered (seA, seB, spA, spB), return
# the 2x2 sub-blocks needed by msm::deltamethod for Se and Sp differences.
.split_pairwise_vcov <- function(fit,
                                 se_cols = c("seA", "seB"),
                                 sp_cols = c("spA", "spB")) {
  V <- stats::vcov(fit)
  coef_names <- rownames(V)
  find <- function(name) {
    idx <- which(coef_names == name)
    if (!length(idx)) stop("Coefficient '", name, "' not found in fit.")
    idx
  }
  se_idx <- vapply(se_cols, find, integer(1))
  sp_idx <- vapply(sp_cols, find, integer(1))
  list(se_cov = as.matrix(V[se_idx, se_idx]),
       sp_cov = as.matrix(V[sp_idx, sp_idx]),
       se_idx = se_idx,
       sp_idx = sp_idx)
}

# -- Random-effect VCV extractor --------------------------------------------
#
# Returns the 2x2 between-study covariance matrix Psi for (logit Se, logit Sp).
# For an equal-variance fit (single shared RE block, model B) this is one 2x2
# matrix.  For an unequal-variance fit (two per-test RE blocks, model E) it is a
# named list of two 2x2 matrices, keyed by the test `levels` when supplied.
.random_vcv <- function(fit, levels = NULL) {
  clean <- function(Psi) {
    dim_ok <- attr(Psi, "dim")
    dn_ok  <- attr(Psi, "dimnames")
    Psi <- unclass(Psi)
    attributes(Psi) <- list(dim = dim_ok, dimnames = dn_ok)
    Psi
  }
  vc <- lme4::VarCorr(fit)
  if (length(vc) == 1L) {
    return(clean(vc[[1]]))
  }
  # Unequal-variance (model E): one block per test.
  blocks <- lapply(vc, clean)
  if (!is.null(levels) && length(levels) == length(blocks)) {
    names(blocks) <- levels
  }
  blocks
}

# -- Fixed-effect extractor (lsens, lspec, SEs) for single-test fits --------
.fixed_se_sp <- function(fit) {
  co <- summary(fit)$coefficients
  nm <- rownames(co)
  if (all(c("sens", "spec") %in% nm)) {
    V <- stats::vcov(fit)
    V <- as.matrix(V[c("sens", "spec"), c("sens", "spec")])
    return(list(
      lsens      = co["sens", 1],
      lspec      = co["spec", 1],
      se_lsens   = co["sens", 2],
      se_lspec   = co["spec", 2],
      vcov_fixed = V
    ))
  }
  stop("Fixed-effect labels 'sens' and 'spec' not found; this helper ",
       "is for single-test fits.  Use .split_pairwise_vcov() for pairwise.")
}
