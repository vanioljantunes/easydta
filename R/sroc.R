# ============================================================================
# sroc.R  -  Summary ROC curve, confidence region, prediction region, AUC
#
# SROC line is derived from the bivariate glmer fit via the Harbord
# equivalence (no HSROC refit required):
#
#   logit(TPR) = lsens - (tau_sens / tau_spec) * (logit(FPR) - (-lspec))
#              [ because logit(FPR) = -logit(Sp) ]
#
# AUC: primary method is the bootstrap from `dmetatools::AUC_boot` (Noma
# et al.), which resamples the bivariate model and yields both a point
# estimate and a bootstrap CI. When `dmetatools` is not installed -- or
# when `auc_method = "trapz"` is requested explicitly -- we still return
# a CI: the trapz path runs a parametric MVN bootstrap on the cached
# fixed-effect VCV (lsens, lspec) and integrates each resampled SROC.
#
# Layout: in-plot summary box and legend box are flush against the
# bottom-right and bottom-left panel corners (axis expansion is disabled).
# Both boxes share the same height; the summary box auto-sizes its width
# to the widest "label value" pair so it never bloats beyond what the
# numbers require.  The sROC curve is solid black and the 95% prediction
# region is dotted red so it is visually distinct from the dashed CI loop.
# ============================================================================

#' SROC plot with 95% confidence and prediction regions + AUC
#'
#' @param fit         A `dta_single` object.
#' @param test        Test name shown in the main title (default "test").
#' @param outcome     Outcome name shown in the main title (default "outcome").
#' @param population  Population name shown in the main title (default
#'   "population"). The title is rendered as
#'   `sROC of <test> to predict <outcome> in <population>`.
#' @param ci          Draw the 95% confidence region?  (default TRUE)
#' @param pred        Draw the 95% prediction region?  (default TRUE; drawn
#'   as a red dotted loop to distinguish it from the dashed CI region).
#' @param labels      Logical. Print study labels next to each triangle?
#'   (default TRUE).  Labels are drawn on top of the panel and are allowed
#'   to overlap each other and the sROC curve -- they never reposition the
#'   underlying study points.
#' @param auc         Compute AUC and attach as attribute?  (default TRUE)
#' @param auc_method  "boot" (default; uses `dmetatools::AUC_boot`) or
#'   "trapz" (numerical integration of the SROC curve with `pracma::trapz`;
#'   automatically used as a fallback if `dmetatools` is not installed).
#'   In both cases an AUC CI is returned.
#' @param B           Number of bootstrap replicates (default 2000).
#' @param conf        Confidence level (default 0.95).
#' @param n_grid      Grid size for the SROC curve (default 200).
#'
#' @return A ggplot object; if `auc = TRUE`, `attr(x, "AUC")` holds the AUC
#'   point estimate and `attr(x, "AUC_CI")` holds the CI.
#' @export
dta_sroc <- function(fit,
                     test       = "test",
                     outcome    = "outcome",
                     population = "population",
                     ci    = TRUE,
                     pred  = TRUE,
                     labels = TRUE,
                     auc   = TRUE,
                     auc_method = c("boot", "trapz"),
                     B     = 2000,
                     conf  = 0.95,
                     n_grid = 200) {
  auc_method <- match.arg(auc_method)
  stopifnot(inherits(fit, "dta_single"))
  f <- .fixed_se_sp(fit$fit)
  Psi <- fit$Psi
  tau_sens <- sqrt(Psi[1, 1])
  tau_spec <- sqrt(Psi[2, 2])

  slope <- if (tau_spec > 0) tau_sens / tau_spec else 0
  fpr_grid <- seq(0.001, 0.999, length.out = n_grid)
  logit_fpr <- stats::qlogis(fpr_grid)
  logit_sp  <- -logit_fpr
  logit_se  <- f$lsens - slope * (logit_sp - f$lspec)
  tpr_grid  <- stats::plogis(logit_se)
  curve <- data.frame(fpr = fpr_grid, tpr = tpr_grid)

  studies <- split(fit$long, fit$long$studlab)
  pts <- do.call(rbind, lapply(names(studies), function(nm) {
    s <- studies[[nm]]
    r_sens <- s[s$sens == 1, ]
    r_spec <- s[s$spec == 1, ]
    tp <- r_sens$true; n1 <- r_sens$n
    tn <- r_spec$true; n0 <- r_spec$n
    data.frame(studlab = nm,
               fpr = (n0 - tn) / n0,
               tpr = tp / n1,
               n_total = n1 + n0,
               stringsAsFactors = FALSE)
  }))
  rownames(pts) <- NULL

  centre <- c(f$lsens, f$lspec)
  cr_df <- if (ci)   .logit_ellipse(centre, f$vcov_fixed, conf) else NULL
  pr_df <- if (pred) .logit_ellipse(centre, f$vcov_fixed + Psi, conf) else NULL

  summary_pt <- data.frame(fpr = 1 - .inv_logit(f$lspec),
                           tpr = .inv_logit(f$lsens))

  # ---- Numbers for the in-plot summary box ---------------------------------
  sens_ci <- .logit_ci(f$lsens, f$se_lsens, conf)
  spec_ci <- .logit_ci(f$lspec, f$se_lspec, conf)
  i2_biv  <- if (!is.null(fit$heterogeneity)) {
    fit$heterogeneity$I2_biv
  } else NA_real_

  auc_val     <- NA_real_
  auc_ci_pair <- NULL
  auc_text    <- "n/a"
  if (auc) {
    auc_res <- .compute_auc(fit, NULL, auc_method, B, conf, curve,
                            "TP", "FP", "FN", "TN",
                            slope = slope, fpr_grid = fpr_grid)
    auc_val <- auc_res$AUC
    if (!is.null(auc_res$CI)) {
      auc_ci_pair <- auc_res$CI
      auc_text <- sprintf("%.3f (%.3f-%.3f)",
                          auc_val, auc_res$CI[1], auc_res$CI[2])
    } else {
      auc_text <- sprintf("%.3f", auc_val)
    }
  }

  fmt2 <- function(x, lo, hi) sprintf("%.2f (%.2f-%.2f)", x, lo, hi)
  i2_text <- if (is.na(i2_biv)) {
    "n/a"
  } else if (i2_biv >= 0.01) {
    sprintf("%.0f%%", 100 * i2_biv)
  } else {
    sprintf("%.1f%%", 100 * i2_biv)
  }

  # Two columns rendered as separate text annotations so the LEFT edge of
  # every label (Sens / Spec / AUC / I2) lines up exactly, and the "="
  # plus value all start at the same x to the right of the labels.
  rows_label <- c("Sens", "Spec", "AUC", "I2")
  rows_value <- paste(
    "=",
    c(
      fmt2(sens_ci["estimate"], sens_ci["lci"], sens_ci["uci"]),
      fmt2(spec_ci["estimate"], spec_ci["lci"], spec_ci["uci"]),
      auc_text,
      i2_text
    )
  )

  # ---- Box geometry --------------------------------------------------------
  # Both boxes share the same height (matches the legend on the left).
  # The summary box on the right auto-sizes its width to the widest
  # "label value" pair so it stays as tight as possible against the
  # bottom-right corner.
  box_h     <- 0.235
  text_size <- 3.2
  # Approximate character width in plot units for ggplot::annotate text
  # at size = 3.2 on a roughly square panel saved at ~150 dpi.  This is a
  # heuristic, not a measurement -- generous enough to avoid clipping.
  char_w  <- 0.0125
  pad_x   <- 0.012   # horizontal padding inside the box
  gap_lv  <- 0.005   # tiny gap between the label column and the value text

  label_w <- max(nchar(rows_label)) * char_w
  value_w <- max(nchar(rows_value)) * char_w
  box_title <- "Summary"
  # Bold title rendered slightly wider than plain text.
  title_w <- nchar(box_title) * char_w * 1.10

  inner_w <- max(label_w + gap_lv + value_w, title_w)
  box_w   <- inner_w + 2 * pad_x

  bx <- list(xmin = 1 - box_w, xmax = 1.00,
             ymin = 0.00,      ymax = box_h)

  # Vertical layout: title row at the top, then 4 evenly-spaced data rows.
  ys_title <- bx$ymax - 0.030
  row_step <- 0.045
  ys_rows  <- seq(ys_title - row_step, by = -row_step, length.out = 4)

  # Two left-anchored columns so every label (Sens / Spec / AUC / I2)
  # starts at the same x and every "= value" string starts at the same x.
  x_label_anchor <- bx$xmin + pad_x
  x_value_anchor <- x_label_anchor + label_w + gap_lv

  # Bottom-left legend box (5 rows: study / sROC / CI / pred / summary).
  lg <- list(xmin = 0.00, xmax = 0.30,
             ymin = 0.00, ymax = box_h)
  lg_rows <- seq(lg$ymax - 0.030, by = -row_step, length.out = 5)
  lg_x_sym   <- lg$xmin + 0.03
  lg_x_label <- lg$xmin + 0.07
  legend_text <- c("Study estimates", "sROC curve",
                   "95% CI region", "95% prediction region",
                   "Summary point")

  title_text <- sprintf("sROC of %s to predict %s in %s",
                        test, outcome, population)

  # ---- Build plot ----------------------------------------------------------
  p <- ggplot2::ggplot()

  if (!is.null(cr_df)) {
    p <- p + ggplot2::geom_path(data = cr_df,
                                ggplot2::aes(x = fpr, y = tpr),
                                colour = "black",
                                linetype = "dashed",
                                linewidth = 0.5)
  }
  if (!is.null(pr_df)) {
    p <- p + ggplot2::geom_path(data = pr_df,
                                ggplot2::aes(x = fpr, y = tpr),
                                colour = "red",
                                linetype = "dotted",
                                linewidth = 0.6)
  }

  p <- p +
    ggplot2::geom_line(data = curve,
                       ggplot2::aes(x = fpr, y = tpr),
                       colour = "black",
                       linetype = "solid",
                       linewidth = 0.6) +
    ggplot2::geom_point(data = pts,
                        ggplot2::aes(x = fpr, y = tpr),
                        shape = 2, size = 2.4, colour = "black",
                        show.legend = FALSE) +
    ggplot2::geom_point(data = summary_pt,
                        ggplot2::aes(x = fpr, y = tpr),
                        shape = 16, size = 3.5, colour = "black")

  if (isTRUE(labels) && nrow(pts) > 0) {
    pts_lab <- pts
    # Left-anchor the label slightly to the right of each point.  Labels
    # are allowed to overlap each other and the sROC line -- they never
    # reposition the underlying study points.
    pts_lab$lab_x <- pts_lab$fpr + 0.012
    p <- p +
      ggplot2::geom_text(data = pts_lab,
                         ggplot2::aes(x = lab_x, y = tpr, label = studlab),
                         size = 3, colour = "black",
                         hjust = 0, vjust = 0.5)
  }

  # Summary box: white-filled rectangle + bold title + bold label column
  # right-anchored at the separator + plain value column left-anchored
  # immediately after, so labels and values sit flush.
  p <- p +
    ggplot2::annotate("rect",
                      xmin = bx$xmin, xmax = bx$xmax,
                      ymin = bx$ymin, ymax = bx$ymax,
                      fill = "white", colour = "black",
                      linewidth = 0.4) +
    ggplot2::annotate("text",
                      x = (bx$xmin + bx$xmax) / 2, y = ys_title,
                      label = box_title, fontface = "bold",
                      size = text_size + 0.4) +
    ggplot2::annotate("text",
                      x = x_label_anchor, y = ys_rows,
                      label = rows_label,
                      fontface = "bold", hjust = 0, size = text_size) +
    ggplot2::annotate("text",
                      x = x_value_anchor, y = ys_rows,
                      label = rows_value,
                      fontface = "plain", hjust = 0, size = text_size)

  # Bottom-left legend box (5 rows: study / sROC / CI / pred / summary).
  seg_half <- 0.018
  p <- p +
    ggplot2::annotate("rect",
                      xmin = lg$xmin, xmax = lg$xmax,
                      ymin = lg$ymin, ymax = lg$ymax,
                      fill = "white", colour = "black",
                      linewidth = 0.4) +
    ggplot2::annotate("point", x = lg_x_sym, y = lg_rows[1],
                      shape = 2, size = 2.4, colour = "black") +
    ggplot2::annotate("segment",
                      x = lg_x_sym - seg_half, xend = lg_x_sym + seg_half,
                      y = lg_rows[2], yend = lg_rows[2],
                      colour = "black", linetype = "solid",
                      linewidth = 0.6) +
    ggplot2::annotate("segment",
                      x = lg_x_sym - seg_half, xend = lg_x_sym + seg_half,
                      y = lg_rows[3], yend = lg_rows[3],
                      colour = "black", linetype = "dashed",
                      linewidth = 0.5) +
    ggplot2::annotate("segment",
                      x = lg_x_sym - seg_half, xend = lg_x_sym + seg_half,
                      y = lg_rows[4], yend = lg_rows[4],
                      colour = "red", linetype = "dotted",
                      linewidth = 0.6) +
    ggplot2::annotate("point", x = lg_x_sym, y = lg_rows[5],
                      shape = 16, size = 3, colour = "black") +
    ggplot2::annotate("text",
                      x = lg_x_label, y = lg_rows,
                      label = legend_text,
                      hjust = 0, size = 3.1)

  p <- p +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::scale_x_continuous(breaks = seq(0, 1, 0.2),
                                limits = c(0, 1),
                                expand = ggplot2::expansion(mult = 0.02)) +
    ggplot2::scale_y_continuous(breaks = seq(0, 1, 0.2),
                                limits = c(0, 1),
                                expand = ggplot2::expansion(mult = 0.02)) +
    ggplot2::labs(x = "False positive rate (1 - Specificity)",
                  y = "Sensitivity",
                  title = title_text) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid       = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   legend.position  = "none",
                   plot.title       = ggplot2::element_text(face = "bold",
                                                            hjust = 0.5,
                                                            size = 13))

  if (auc) {
    attr(p, "AUC") <- auc_val
    if (!is.null(auc_ci_pair)) attr(p, "AUC_CI") <- auc_ci_pair
  }
  p
}

# AUC dispatcher: prefer dmetatools::AUC_boot; fall back to trapezoid.
# Both paths return a CI -- the trapz path uses a parametric MVN bootstrap
# on the cached fixed-effect VCV so users without dmetatools still see one.
.compute_auc <- function(fit, data, auc_method, B, conf,
                         curve, tp, fp, fn, tn,
                         slope = NULL, fpr_grid = NULL) {
  if (auc_method == "boot") {
    if (!requireNamespace("dmetatools", quietly = TRUE) ||
        !requireNamespace("mada", quietly = TRUE)) {
      warning("dmetatools (and its dep 'mada') not installed; falling back ",
              "to auc_method = 'trapz'.\n",
              "  Install with:\n",
              "    install.packages('mada')\n",
              "    remotes::install_github('nomahi/dmetatools')",
              call. = FALSE)
      return(.compute_auc(fit, data, "trapz", B, conf, curve,
                          tp, fp, fn, tn,
                          slope = slope, fpr_grid = fpr_grid))
    }
    # AUC_boot calls mada::reitsma() and MASS::mvrnorm() without namespace
    # qualifiers, so both must be attached (not merely loaded).
    for (pkg in c("mada", "MASS")) {
      if (!paste0("package:", pkg) %in% search()) {
        suppressPackageStartupMessages(attachNamespace(pkg))
      }
    }
    counts <- .extract_counts(fit, data, tp, fp, fn, tn)
    res <- dmetatools::AUC_boot(counts$TP, counts$FP, counts$FN, counts$TN,
                                B = B, alpha = conf)
    return(list(AUC = unname(res$AUC),
                CI  = unname(res$CI),
                method = "dmetatools::AUC_boot"))
  }

  ord <- order(curve$fpr)
  AUC <- pracma::trapz(curve$fpr[ord], curve$tpr[ord])

  CI <- .trapz_auc_ci(fit, slope, fpr_grid, B, conf)
  list(AUC = AUC, CI = CI, method = "trapezoid")
}

# Parametric MVN bootstrap CI for the trapz AUC.  Samples (lsens, lspec)
# from N(centre, vcov_fixed) and integrates each resampled SROC over the
# same FPR grid.  Slope (tau_sens / tau_spec) is held fixed at the point
# estimate -- a reasonable approximation when the random-effect VCV is
# itself uncertain and we just need a defensible interval.
.trapz_auc_ci <- function(fit, slope, fpr_grid, B, conf) {
  if (is.null(slope) || is.null(fpr_grid) || B <= 1) return(NULL)
  f <- .fixed_se_sp(fit$fit)
  centre <- c(f$lsens, f$lspec)
  V <- f$vcov_fixed
  draws <- tryCatch(
    MASS::mvrnorm(n = B, mu = centre, Sigma = V),
    error = function(e) NULL
  )
  if (is.null(draws)) return(NULL)
  logit_fpr <- stats::qlogis(fpr_grid)
  logit_sp  <- -logit_fpr
  ord <- order(fpr_grid)
  aucs <- apply(draws, 1, function(p) {
    logit_se <- p[1] - slope * (logit_sp - p[2])
    tpr <- stats::plogis(logit_se)
    pracma::trapz(fpr_grid[ord], tpr[ord])
  })
  alpha <- 1 - conf
  unname(stats::quantile(aucs, c(alpha / 2, 1 - alpha / 2), na.rm = TRUE))
}

# Get the raw TP/FP/FN/TN vectors: prefer user-supplied wide `data`,
# otherwise reconstruct from the long data carried in the fit.
.extract_counts <- function(fit, data, tp, fp, fn, tn) {
  if (!is.null(data)) {
    req <- c(tp, fp, fn, tn)
    miss <- setdiff(req, names(data))
    if (length(miss)) stop("data is missing columns: ",
                           paste(miss, collapse = ", "))
    return(list(TP = as.integer(data[[tp]]),
                FP = as.integer(data[[fp]]),
                FN = as.integer(data[[fn]]),
                TN = as.integer(data[[tn]])))
  }
  studies <- split(fit$long, fit$long$studlab)
  df <- do.call(rbind.data.frame, lapply(studies, function(s) {
    r_sens <- s[s$sens == 1, ]
    r_spec <- s[s$spec == 1, ]
    data.frame(TP = r_sens$true,
               FP = r_spec$n - r_spec$true,
               FN = r_sens$n - r_sens$true,
               TN = r_spec$true)
  }))
  list(TP = df$TP, FP = df$FP, FN = df$FN, TN = df$TN)
}
