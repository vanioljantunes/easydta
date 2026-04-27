# ============================================================================
# heterogeneity.R  -  Internal between-study variance helpers
#
# Reference:
#   Zhou Y, Dendukuri N (2014). "Statistics for quantifying heterogeneity
#   in univariate and bivariate meta-analyses of binary data: the case of
#   meta-analyses of diagnostic accuracy."  Stat Med. 33(16):2701-2717.
#
# Cochrane ref: Handbook Supp. Mat. 1, Section 10.2.5 advises against the
#   naive Higgins I^2 on Se/Sp separately because of the threshold effect.
#   The Zhou-Dendukuri formulation below is defined on the bivariate model
#   and is methodologically defensible in DTA.
#
# These helpers are called from inside dta_fit_single() so that the result
# of the fit already carries the heterogeneity summary; there is no
# separately exported user-facing dta_heterogeneity() function.
# ============================================================================

.compute_heterogeneity <- function(fit_obj, conf = 0.95) {
  Psi <- fit_obj$Psi
  tau_sens <- sqrt(Psi[1, 1])
  tau_spec <- sqrt(Psi[2, 2])
  rho <- if (tau_sens > 0 && tau_spec > 0) {
    Psi[1, 2] / (tau_sens * tau_spec)
  } else NA_real_

  Sigma_bar <- .avg_within_study_vcov(fit_obj$long)

  I2_sens <- .zd_i2_scalar(Psi[1, 1], Sigma_bar[1, 1])
  I2_spec <- .zd_i2_scalar(Psi[2, 2], Sigma_bar[2, 2])
  I2_biv  <- .zd_i2_scalar(sum(diag(Psi)), sum(diag(Sigma_bar)))

  f <- .fixed_se_sp(fit_obj$fit)
  V_pred <- f$vcov_fixed + Psi
  pred_region <- .logit_ellipse(
    centre = c(f$lsens, f$lspec),
    V = V_pred,
    conf = conf
  )

  list(
    tau_sens    = tau_sens,
    tau_spec    = tau_spec,
    rho         = rho,
    Psi         = Psi,
    Sigma_bar   = Sigma_bar,
    I2_sens     = I2_sens,
    I2_spec     = I2_spec,
    I2_biv      = I2_biv,
    conf        = conf,
    pred_region = pred_region
  )
}

.avg_within_study_vcov <- function(long) {
  studies <- split(long, long$studlab)
  var_se <- numeric(length(studies))
  var_sp <- numeric(length(studies))
  for (i in seq_along(studies)) {
    s <- studies[[i]]
    r_sens <- s[s$sens == 1, ]
    r_spec <- s[s$spec == 1, ]
    tp <- r_sens$true; n1 <- r_sens$n
    tn <- r_spec$true; n0 <- r_spec$n

    if (tp == 0 || tp == n1) { tp <- tp + 0.5; n1 <- n1 + 1 }
    if (tn == 0 || tn == n0) { tn <- tn + 0.5; n0 <- n0 + 1 }

    p <- tp / n1
    q <- tn / n0
    var_se[i] <- 1 / (n1 * p * (1 - p))
    var_sp[i] <- 1 / (n0 * q * (1 - q))
  }
  matrix(c(mean(var_se), 0, 0, mean(var_sp)),
         nrow = 2, ncol = 2,
         dimnames = list(c("sens", "spec"), c("sens", "spec")))
}

.zd_i2_scalar <- function(tau2, sigma2) {
  denom <- tau2 + sigma2
  if (denom <= 0) return(NA_real_)
  tau2 / denom
}

.print_heterogeneity_block <- function(h, digits = 3) {
  cat("Heterogeneity (Zhou-Dendukuri 2014, bivariate I^2):\n")
  cat("  tau_sens = ", format(h$tau_sens, digits = digits),
      "    tau_spec = ", format(h$tau_spec, digits = digits),
      "    rho = ",      format(h$rho,      digits = digits), "\n", sep = "")
  cat("  I^2(sens) = ", format(100 * h$I2_sens, digits = digits), "%",
      "    I^2(spec) = ", format(100 * h$I2_spec, digits = digits), "%",
      "    I^2(biv) = ", format(100 * h$I2_biv,  digits = digits), "%\n",
      sep = "")
  cat("  (joint statistic, unaffected by threshold effect)\n")
  invisible(h)
}
