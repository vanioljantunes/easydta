# ============================================================================
# sroc.R  -  Summary ROC curve, confidence region, prediction region, AUC
#
# SROC line is derived from the bivariate glmer fit via the Harbord
# equivalence (no HSROC refit required):
#
#   logit(TPR) = lsens - (tau_sens / tau_spec) * (logit(FPR) - (-lspec))
#              [ because logit(FPR) = -logit(Sp) ]
#
# AUC: the trapezoidal integral of the fitted SROC curve over a fixed FPR
# grid (`pracma::trapz`). The CI is a parametric MVN bootstrap: draw
# (lsens, lspec) from N(centre, vcov_fixed) -- slope held at the point
# estimate -- integrate each resampled SROC, take quantiles. For a pairwise
# comparison the AUC difference (dAUC) reuses the same per-arm draws: the
# arms are bootstrapped independently and differenced (no external deps).
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
#' @param test.label  Test name shown in the main title (default "test").
#' @param outcome     Outcome name shown in the main title (default "outcome").
#' @param population  Population name shown in the main title (default
#'   "population"). The title is rendered as
#'   `sROC of <test.label> to predict <outcome> in <population>`.
#' @param ci          Draw the 95% confidence region?  (default TRUE)
#' @param pred        Draw the 95% prediction region?  (default TRUE; drawn
#'   as a red dotted loop to distinguish it from the dashed CI region).
#' @param labels      Logical. Print study labels next to each triangle?
#'   (default FALSE).  Labels are drawn on top of the panel and are allowed
#'   to overlap each other and the sROC curve -- they never reposition the
#'   underlying study points.
#' @param group       Optional grouping of the study points, as a vector named
#'   by study label (or an unnamed vector in study order).  Each level gets its
#'   own plot symbol and its own legend row.  `NA` entries keep the default
#'   triangle.  Default `NULL` (all studies drawn alike).
#' @param group.suffix Text appended to each group level in the legend
#'   (default `" studies"`, so `"curvilinear"` reads `"curvilinear studies"`).
#' @param shapes      Plot symbols used for the `group` levels, in order
#'   (default open circle, open square, open diamond, open triangle-down,
#'   cross).
#' @param colors      Point colours for the `group` levels, in order or named
#'   by level (like `shapes`).  Defaults to a colour-blind-safe palette so
#'   grouped points differ by colour as well as symbol; ignored when `group`
#'   is `NULL`.
#' @param legend.pos  Where the symbol legend box sits: `"bottom"` (bottom
#'   centre, the default), `"bottomleft"`, or `"topright"`.
#' @param title.size  Font size of the two-line plot title (default 11;
#'   `dta_sroc_pair()` drops it to 10 for the half-width panels).
#' @param auc         Compute AUC and attach as attribute?  (default TRUE)
#' @param auc_ci      Compute a CI for the AUC (parametric MVN bootstrap)?
#'   (default TRUE). If FALSE only the trapezoidal point estimate is returned.
#' @param B           Number of MVN-bootstrap replicates for the AUC CI
#'   (default 2000).
#' @param conf        Confidence level (default 0.95).
#' @param n_grid      Grid size for the SROC curve (default 200).
#'
#' @return A ggplot object; if `auc = TRUE`, `attr(x, "AUC")` holds the AUC
#'   point estimate and `attr(x, "AUC_CI")` holds the CI.
#' @examples
#' data(anti_ccp2)
#' fit <- dta_fit_single(anti_ccp2, wide = TRUE)
#' dta_sroc(fit, test.label = "anti-CCP2",
#'          outcome = "rheumatoid arthritis", population = "adults")
#' @export
dta_sroc <- function(fit,
                     test.label = "test",
                     outcome    = "outcome",
                     population = "population",
                     ci    = TRUE,
                     pred  = TRUE,
                     labels = FALSE,
                     group  = NULL,
                     group.suffix = " studies",
                     shapes = c(1, 0, 5, 6, 4),
                     colors = c("#0072B2", "#D55E00", "#009E73",
                                "#CC79A7", "#E69F00"),
                     legend.pos = c("bottom", "bottomleft", "topright"),
                     title.size = 11,
                     auc   = TRUE,
                     auc_ci = TRUE,
                     B     = 2000,
                     conf  = 0.95,
                     n_grid = 200,
                     auc_override = NULL) {
  stopifnot(inherits(fit, "dta_single"))
  legend.pos <- match.arg(legend.pos)
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

  # Optional grouping of the study points (e.g. measurement plane): one plot
  # symbol and one legend row per level.  `group` is a vector named by study
  # label; unnamed vectors are matched positionally to `pts$studlab`.
  if (is.null(group)) {
    pts$grp   <- NA_character_
    grp_lv    <- character(0)
    grp_shape <- integer(0)
    grp_col   <- character(0)
  } else {
    g <- if (!is.null(names(group))) as.character(group[pts$studlab]) else as.character(group)
    if (length(g) != nrow(pts))
      stop("`group` must have one entry per study in the fit.")
    pts$grp   <- g
    grp_lv    <- sort(unique(g[!is.na(g)]))
    # Named `shapes` pin a symbol to a level; unnamed ones are taken in order.
    grp_shape <- if (!is.null(names(shapes))) shapes[grp_lv] else shapes[seq_along(grp_lv)]
    if (anyNA(grp_shape))
      stop("`shapes` has no symbol for: ",
           paste(grp_lv[is.na(grp_shape)], collapse = ", "))
    names(grp_shape) <- grp_lv
    # Named `colors` pin a colour to a level; unnamed ones are taken in order.
    grp_col <- if (!is.null(names(colors))) colors[grp_lv] else colors[seq_along(grp_lv)]
    if (anyNA(grp_col))
      stop("`colors` has no colour for: ",
           paste(grp_lv[is.na(grp_col)], collapse = ", "))
    names(grp_col) <- grp_lv
  }

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
    if (!is.null(auc_override)) {
      auc_val     <- auc_override$AUC
      auc_ci_pair <- auc_override$CI
    } else if (!isTRUE(auc_ci)) {
      ord <- order(curve$fpr)
      auc_val <- pracma::trapz(curve$fpr[ord], curve$tpr[ord])
    } else {
      auc_res <- .compute_auc(fit, B, conf, curve,
                              slope = slope, fpr_grid = fpr_grid)
      auc_val     <- auc_res$AUC
      auc_ci_pair <- auc_res$CI
    }
    auc_text <- if (!is.null(auc_ci_pair)) {
      sprintf("%.3f (%.3f-%.3f)", auc_val, auc_ci_pair[1], auc_ci_pair[2])
    } else {
      sprintf("%.3f", auc_val)
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

  # Bottom-left legend box: one row per study symbol (one per `group` level,
  # or a single "Study estimates" row when ungrouped), then sROC / CI / pred /
  # summary.  The box grows with the number of rows and the longest label.
  if (length(grp_lv)) {
    pt_pch  <- unname(grp_shape)
    pt_col  <- unname(grp_col)
    pt_text <- paste0(grp_lv, group.suffix)
    if (any(is.na(pts$grp))) {
      pt_pch  <- c(pt_pch, 2)
      pt_col  <- c(pt_col, "black")
      pt_text <- c(pt_text, "Study estimates")
    }
  } else {
    pt_pch  <- 2
    pt_col  <- "black"
    pt_text <- "Study estimates"
  }
  n_pt <- length(pt_pch)

  # The prediction-region row only appears when the region is drawn.
  if (pred) {
    legend_text <- c(pt_text, "sROC curve", "95% CI region",
                     "95% prediction region", "Summary point")
    lg_pch <- c(pt_pch, NA, NA, NA, 16)
    lg_lty <- c(rep(NA, n_pt), "solid", "dashed", "dotted", NA)
    lg_col <- c(pt_col, "black", "black", "red", "black")
    lg_lwd <- c(rep(NA, n_pt), 0.6, 0.5, 0.6, NA)
    lg_size <- c(rep(2.4, n_pt), NA, NA, NA, 3)
  } else {
    legend_text <- c(pt_text, "sROC curve", "95% CI region", "Summary point")
    lg_pch <- c(pt_pch, NA, NA, 16)
    lg_lty <- c(rep(NA, n_pt), "solid", "dashed", NA)
    lg_col <- c(pt_col, "black", "black", "black")
    lg_lwd <- c(rep(NA, n_pt), 0.6, 0.5, NA)
    lg_size <- c(rep(2.4, n_pt), NA, NA, 3)
  }
  n_lg <- length(legend_text)

  # Legend box geometry, then anchored by `legend.pos`: away from the
  # top-left cloud of study points, which the old bottom-left box could hide
  # on high-specificity data.
  lg_w <- max(0.30, 0.07 + max(nchar(legend_text)) * char_w + 0.012)
  lg_h <- max(box_h, 0.055 + (n_lg - 1) * row_step)
  lg <- switch(legend.pos,
    bottomleft = list(xmin = 0.00,           xmax = lg_w,
                      ymin = 0.00,           ymax = lg_h),
    bottom     = list(xmin = (1 - lg_w) / 2, xmax = (1 + lg_w) / 2,
                      ymin = 0.00,           ymax = lg_h),
    topright   = list(xmin = 1 - lg_w,       xmax = 1.00,
                      ymin = 1 - lg_h,       ymax = 1.00))
  lg_x_sym   <- lg$xmin + 0.03
  lg_x_label <- lg$xmin + 0.07
  lg_rows <- seq(lg$ymax - 0.030, by = -row_step, length.out = n_lg)

  # Two lines: the title has to fit inside a half-width panel in dta_sroc_pair.
  title_text <- sprintf("sROC of %s\nto predict %s in %s",
                        test.label, outcome, population)

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
    ggplot2::geom_point(data = if (length(grp_lv)) pts[is.na(pts$grp), ] else pts,
                        ggplot2::aes(x = fpr, y = tpr),
                        shape = 2, size = 2.4, colour = "black",
                        show.legend = FALSE) +
    ggplot2::geom_point(data = summary_pt,
                        ggplot2::aes(x = fpr, y = tpr),
                        shape = 16, size = 3.5, colour = "black")

  for (lv in grp_lv) {
    p <- p + ggplot2::geom_point(data = pts[!is.na(pts$grp) & pts$grp == lv, ],
                                 ggplot2::aes(x = fpr, y = tpr),
                                 shape = grp_shape[[lv]], size = 2.4,
                                 colour = grp_col[[lv]], show.legend = FALSE)
  }

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

  seg_half <- 0.018
  p <- p +
    ggplot2::annotate("rect",
                      xmin = lg$xmin, xmax = lg$xmax,
                      ymin = lg$ymin, ymax = lg$ymax,
                      fill = "white", colour = "black",
                      linewidth = 0.4)

  for (i in seq_len(n_lg)) {
    p <- if (!is.na(lg_pch[i])) {
      p + ggplot2::annotate("point", x = lg_x_sym, y = lg_rows[i],
                            shape = lg_pch[i], size = lg_size[i],
                            colour = lg_col[i])
    } else {
      p + ggplot2::annotate("segment",
                            x = lg_x_sym - seg_half, xend = lg_x_sym + seg_half,
                            y = lg_rows[i], yend = lg_rows[i],
                            colour = lg_col[i], linetype = lg_lty[i],
                            linewidth = lg_lwd[i])
    }
  }

  p <- p +
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
                                                            lineheight = 1.15,
                                                            size = title.size))

  if (auc) {
    attr(p, "AUC") <- auc_val
    if (!is.null(auc_ci_pair)) attr(p, "AUC_CI") <- auc_ci_pair
  }
  p
}

#' Side-by-side SROC plot for a pairwise / covariate comparison
#'
#' Convenience wrapper around two `dta_sroc()` calls that arranges the
#' intervention (`.e`) and control (`.c`) arms next to each other and,
#' below the panels, prints a small Cochrane-style differences table:
#' for each measure (Sens, Spec, AUC) the per-arm estimate (95% CI),
#' the absolute difference `.e - .c` (95% CI), and a p-value.
#'
#' Statistical inference follows the Cochrane Handbook for DTA Reviews
#' (Appendix 12, Takwoingi et al. 2023):
#' * Sens / Spec p-values come from likelihood-ratio tests of the
#'   bivariate model with the relevant fixed effect constrained equal
#'   between arms (rows 2 and 3 of `x$compare$lr_tests`); the absolute
#'   difference (with delta-method CI) comes from `x$compare$differences`.
#' * AUC inference (when `auc_ic = TRUE`) is a parametric MVN bootstrap:
#'   each arm's AUC + 95% CI come from resampling `(lsens, lspec)` and
#'   integrating each SROC (the `dta_sroc()` method). The difference
#'   `dAUC = AUC.e - AUC.c`, its 95% CI, and the p-value for `dAUC = 0`
#'   are formed by differencing the two arms' bootstrap draws. The arms
#'   are drawn independently, so the dAUC CI ignores within-study
#'   correlation and is mildly conservative. No external packages needed.
#'
#' @param x         A `dta_pairwise_result` (from `dta_pairwise()` /
#'   `dta_compare_tests()`) *or* a named list of two `dta_single` fits.
#' @param arm.e,arm.c Names of the intervention (`.e`) and control (`.c`)
#'   arms in `x$arms`.  For a `dta_pairwise_result` they default to the
#'   `intervention` / `control` arms recorded in the result, so you can omit
#'   them; required when `x` is a plain named list.
#' @param test.label.e,test.label.c Test labels shown in each panel title and
#'   in the differences table header.  Default to `arm.e` / `arm.c`.
#' @param outcome,population Shared title fields passed through to
#'   `dta_sroc()`.
#' @param table     Logical. Show the differences table below the panels?
#'   (default `TRUE`).
#' @param auc_ic    Logical. Compute AUC point estimates + confidence
#'   intervals + dAUC inference?  (default `TRUE`.)  When `FALSE`:
#'   per-panel summary boxes show just the AUC point estimate (no CI),
#'   and the differences table omits the AUC row entirely -- skipping
#'   the bootstrap is much faster and useful for previewing.
#' @param B         MVN-bootstrap replicates for the AUC / dAUC CIs (default
#'   2000); ignored when `auc_ic = FALSE`.
#' @param conf      Confidence level (default 0.95).
#' @param ncol      Number of columns for the SROC panel row (default 2).
#' @param ...       Extra arguments forwarded to both `dta_sroc()` calls
#'   (e.g. `labels`, `pred`).  `auc_ci` is set automatically from `auc_ic`
#'   and should not be overridden.
#'
#' @return A `gtable` of class `"dta_sroc_pair"`; its `print()` method draws
#'   it on the current device, so it renders automatically when returned at
#'   the top level (no `grid::grid.draw()` needed).
#'   `attr(., "panels")` holds the underlying ggplot list and
#'   `attr(., "diff_table")` holds the rendered data.frame.
#' @examples
#' \donttest{
#' data(schuetz)
#' res <- dta_pairwise(schuetz, studlab = "studlab",
#'                     intervention.label = "CT", control.label = "MRI")
#' dta_sroc_pair(res, auc_ic = FALSE)   # arms taken from the result
#' }
#' @export
dta_sroc_pair <- function(x,
                          arm.e = NULL, arm.c = NULL,
                          test.label.e = NULL,
                          test.label.c = NULL,
                          outcome    = "outcome",
                          population = "population",
                          table  = TRUE,
                          auc_ic = TRUE,
                          B      = 2000,
                          conf   = 0.95,
                          ncol   = 2,
                          ...) {
  arms <- if (inherits(x, c("dta_pairwise_result", "dta_compare"))) {
    x$arms
  } else x
  if (!is.list(arms) || is.null(names(arms)))
    stop("`x` must be a dta_pairwise_result or a named list of dta_single fits.")

  # Default the arms to the intervention / control roles recorded on the
  # result, so `dta_sroc_pair(res)` works without naming the arms.
  if (is.null(arm.e)) arm.e <- x$labels$intervention
  if (is.null(arm.c)) arm.c <- x$labels$control
  if (is.null(arm.e) || is.null(arm.c))
    stop("Provide `arm.e` and `arm.c` (no intervention/control roles found).")
  if (is.null(test.label.e)) test.label.e <- arm.e
  if (is.null(test.label.c)) test.label.c <- arm.c

  miss <- setdiff(c(arm.e, arm.c), names(arms))
  if (length(miss))
    stop("Arms not found in x$arms: ", paste(miss, collapse = ", "))

  fit_e <- arms[[arm.e]]
  fit_c <- arms[[arm.c]]

  # Per-arm AUC + dAUC inference up-front so each panel can reuse the same
  # numbers via `auc_override` and the differences table gets one consistent
  # dAUC + p-value.  Skipped entirely when auc_ic = FALSE.
  auc_pair <- if (isTRUE(auc_ic)) {
    .compute_auc_pair(fit_e, fit_c, B = B, conf = conf)
  } else NULL

  panel_args <- list(...)
  panel_args$auc    <- isTRUE(auc_ic)
  panel_args$auc_ci <- isTRUE(auc_ic)
  # Each panel is half the device width, so the two-line title is set smaller
  # here than in a standalone dta_sroc().
  if (is.null(panel_args$title.size)) panel_args$title.size <- 10

  p_e <- do.call(dta_sroc, c(list(fit_e,
                                  test.label = test.label.e,
                                  outcome    = outcome,
                                  population = population,
                                  auc_override = auc_pair$arm.e),
                             panel_args))
  p_c <- do.call(dta_sroc, c(list(fit_c,
                                  test.label = test.label.c,
                                  outcome    = outcome,
                                  population = population,
                                  auc_override = auc_pair$arm.c),
                             panel_args))

  panels <- gridExtra::arrangeGrob(p_e, p_c, ncol = ncol)

  diff_df <- .sroc_pair_diff_table(fit_e, fit_c, p_e, p_c,
                                   test.label.e, test.label.c, arm.e, arm.c,
                                   x, auc_pair, auc_ic, conf)

  if (isTRUE(table)) {
    tbl_grob <- gridExtra::tableGrob(
      diff_df,
      rows = NULL,
      theme = gridExtra::ttheme_minimal(
        core    = list(fg_params = list(cex = 0.85)),
        colhead = list(fg_params = list(cex = 0.85, fontface = "bold"))
      )
    )
    # Wrap the differences table in a single bordered container box that
    # hugs the table extent (header + body), so it reads as one panel.
    tbl_grob <- gtable::gtable_add_grob(
      tbl_grob,
      grobs = grid::rectGrob(gp = grid::gpar(fill = NA, col = "black",
                                             lwd = 1.2)),
      t = 1, b = nrow(tbl_grob), l = 1, r = ncol(tbl_grob),
      name = "container-border"
    )
    # Size the table row from the table itself so extra rows never clip.
    tbl_h <- grid::grobHeight(tbl_grob) + grid::unit(40, "pt")
    g <- gridExtra::arrangeGrob(panels, tbl_grob, nrow = 2,
                                heights = grid::unit.c(
                                  grid::unit(1, "npc") - tbl_h, tbl_h))
  } else {
    g <- panels
  }

  attr(g, "panels")     <- list(.e = p_e, .c = p_c)
  attr(g, "diff_table") <- diff_df
  class(g) <- c("dta_sroc_pair", class(g))
  g
}

#' @export
print.dta_sroc_pair <- function(x, ...) {
  grid::grid.newpage()
  grid::grid.draw(x)
  invisible(x)
}

# Per-arm SROC geometry + AUC point estimate, matching dta_sroc()'s
# construction (Harbord slope = tau_sens / tau_spec; FPR grid 0.001..0.999).
.arm_auc_geom <- function(fit, n_grid = 200) {
  f   <- .fixed_se_sp(fit$fit)
  Psi <- fit$Psi
  tau_sens <- sqrt(Psi[1, 1])
  tau_spec <- sqrt(Psi[2, 2])
  slope    <- if (tau_spec > 0) tau_sens / tau_spec else 0
  fpr_grid <- seq(0.001, 0.999, length.out = n_grid)
  logit_sp <- -stats::qlogis(fpr_grid)
  tpr      <- stats::plogis(f$lsens - slope * (logit_sp - f$lspec))
  ord      <- order(fpr_grid)
  list(slope = slope, fpr_grid = fpr_grid,
       AUC = pracma::trapz(fpr_grid[ord], tpr[ord]))
}

# Per-arm AUC + dAUC inference, mada-free.  Each arm's AUC CI is the
# parametric MVN bootstrap from .trapz_auc_draws(); the dAUC estimate,
# CI, and p-value come from differencing the two arms' independent draws
# (so the dAUC CI ignores within-study correlation -- mildly conservative).
# Returns list($arm.e, $arm.c each list(AUC, CI), $diff list(est, ci, p)).
.compute_auc_pair <- function(fit_e, fit_c, B = 2000, conf = 0.95,
                              n_grid = 200) {
  ge <- .arm_auc_geom(fit_e, n_grid)
  gc <- .arm_auc_geom(fit_c, n_grid)
  draws_e <- .trapz_auc_draws(fit_e, ge$slope, ge$fpr_grid, B)
  draws_c <- .trapz_auc_draws(fit_c, gc$slope, gc$fpr_grid, B)

  alpha <- 1 - conf
  qci <- function(v) unname(stats::quantile(v, c(alpha / 2, 1 - alpha / 2),
                                            na.rm = TRUE))

  arm.e <- list(AUC = ge$AUC, CI = if (!is.null(draws_e)) qci(draws_e) else NULL)
  arm.c <- list(AUC = gc$AUC, CI = if (!is.null(draws_c)) qci(draws_c) else NULL)

  diff <- if (!is.null(draws_e) && !is.null(draws_c)) {
    d <- draws_e - draws_c
    p <- min(1, 2 * min(mean(d <= 0), mean(d >= 0)))
    list(est = ge$AUC - gc$AUC, ci = qci(d), p = p)
  } else {
    list(est = ge$AUC - gc$AUC, ci = c(NA_real_, NA_real_), p = NA_real_)
  }
  list(arm.e = arm.e, arm.c = arm.c, diff = diff)
}

# Build the differences table shown below dta_sroc_pair() panels.
# Columns: Measure, Effect Size (.e - .c), P-value; a final I2 row carries
# the residual Cochran Q heterogeneity of the paired data.
# Sens/Spec p-values come from x$compare$lr_tests (Cochrane LR tests).
# Sens/Spec diff CIs come from x$compare$differences (delta method),
# sign-flipped if the comparison was estimated as (.c - .e) by the
# pairwise model's level ordering.
# AUC row uses auc_pair (per-arm AUC + CI and the bootstrap dAUC est/CI/p);
# omitted when auc_ic = FALSE.
.sroc_pair_diff_table <- function(fit_e, fit_c, p_e, p_c,
                                  test.label.e, test.label.c, arm.e, arm.c,
                                  x, auc_pair, auc_ic, conf) {
  # Result figures rounded to 2 decimals; p-values always to 3.
  fmt <- function(est, lo, hi) {
    if (is.na(est)) return("n/a")
    if (is.na(lo) || is.na(hi)) return(sprintf("%.2f", est))
    sprintf("%.2f (%.2f, %.2f)", est, lo, hi)
  }
  fmt_p <- function(p) {
    if (is.null(p) || is.na(p)) return("n/a")
    if (p < 0.001) "<0.001" else sprintf("%.3f", p)
  }
  arm_row <- function(fit) {
    f <- .fixed_se_sp(fit$fit)
    list(sens = .logit_ci(f$lsens, f$se_lsens, conf),
         spec = .logit_ci(f$lspec, f$se_lspec, conf))
  }
  re <- arm_row(fit_e); rc <- arm_row(fit_c)

  diff_sens <- diff_spec <- list(est = NA_real_, lo = NA_real_, hi = NA_real_)
  p_sens <- p_spec <- NA_real_

  cmp  <- if (inherits(x, "dta_pairwise_result")) x$compare else NULL
  pair <- if (inherits(x, "dta_pairwise_result")) x$pair    else NULL

  if (!is.null(cmp) && !is.null(cmp$differences) && !is.null(pair)) {
    levs <- pair$levels
    s <- if (identical(levs[1], arm.e)) 1 else -1
    d <- cmp$differences
    rs <- d[d$measure == "Absolute diff Sens", ][1, ]
    rp <- d[d$measure == "Absolute diff Spec", ][1, ]
    diff_sens <- list(est = s * rs$estimate,
                      lo  = if (s == 1) rs$lci else -rs$uci,
                      hi  = if (s == 1) rs$uci else -rs$lci)
    diff_spec <- list(est = s * rp$estimate,
                      lo  = if (s == 1) rp$lci else -rp$uci,
                      hi  = if (s == 1) rp$uci else -rp$lci)
    lr <- cmp$lr_tests
    p_sens <- lr$p_value[grepl("Se", lr$comparison)][1]
    p_spec <- lr$p_value[grepl("Sp", lr$comparison)][1]
  } else {
    diff_sens$est <- unname(re$sens["estimate"] - rc$sens["estimate"])
    diff_spec$est <- unname(re$spec["estimate"] - rc$spec["estimate"])
  }

  out <- data.frame(
    Measure = c("Sensitivity", "Specificity"),
    Diff = c(fmt(diff_sens$est, diff_sens$lo, diff_sens$hi),
             fmt(diff_spec$est, diff_spec$lo, diff_spec$hi)),
    P = c(fmt_p(p_sens), fmt_p(p_spec)),
    stringsAsFactors = FALSE
  )

  if (isTRUE(auc_ic)) {
    auc_e    <- if (!is.null(auc_pair$arm.e)) auc_pair$arm.e$AUC else attr(p_e, "AUC")
    auc_e_ci <- if (!is.null(auc_pair$arm.e)) auc_pair$arm.e$CI  else attr(p_e, "AUC_CI")
    auc_c    <- if (!is.null(auc_pair$arm.c)) auc_pair$arm.c$AUC else attr(p_c, "AUC")
    auc_c_ci <- if (!is.null(auc_pair$arm.c)) auc_pair$arm.c$CI  else attr(p_c, "AUC_CI")

    if (!is.null(auc_pair$diff)) {
      d_est <- auc_pair$diff$est; d_lo <- auc_pair$diff$ci[1]
      d_hi  <- auc_pair$diff$ci[2]; p_auc <- auc_pair$diff$p
    } else {
      d_est <- if (!is.na(auc_e) && !is.na(auc_c)) auc_e - auc_c else NA_real_
      d_lo  <- NA_real_; d_hi <- NA_real_; p_auc <- NA_real_
    }

    out <- rbind(out, data.frame(
      Measure = "AUC",
      Diff = fmt(d_est, d_lo, d_hi),
      P = fmt_p(p_auc),
      stringsAsFactors = FALSE
    ))
  }

  # Residual Cochran Q across both arms: fixed-effect logit pooling per
  # measure within each arm, so Q tests between-study heterogeneity beyond
  # the test difference itself (df = sum of k - 1 over the four strata).
  # I2 = (Q - df) / Q (Higgins 2003); the naive form, reported here because
  # the paired table has no threshold axis to condition on.
  q_stat <- 0; q_df <- 0
  for (fit in list(fit_e, fit_c)) {
    long <- fit$long
    for (msr in c("sens", "spec")) {
      r <- long[long[[msr]] == 1, , drop = FALSE]
      k <- r$true; n <- r$n
      cc <- k == 0 | k == n            # continuity correction on the boundary
      k[cc] <- k[cc] + 0.5
      n[cc] <- n[cc] + 1
      th <- stats::qlogis(k / n)
      w  <- k * (n - k) / n            # 1 / var(logit p)
      th_bar <- sum(w * th) / sum(w)
      q_stat <- q_stat + sum(w * (th - th_bar)^2)
      q_df   <- q_df + (length(th) - 1)
    }
  }
  i2  <- if (q_stat > 0) max(0, (q_stat - q_df) / q_stat) else NA_real_
  p_q <- if (q_df > 0) stats::pchisq(q_stat, df = q_df, lower.tail = FALSE)
         else NA_real_
  out <- rbind(out, data.frame(
    Measure = "I2",
    Diff = if (is.na(i2)) "n/a" else sprintf("%.1f%%", 100 * i2),
    P = fmt_p(p_q),
    stringsAsFactors = FALSE
  ))

  # The arm names already head the two panels above the table, so the
  # effect-size column stays short, without "(A - B)".
  names(out) <- c("Measure", "Effect Size", "P-value")
  out
}

# AUC = trapezoidal integral of the fitted SROC curve; CI = parametric MVN
# bootstrap (.trapz_auc_ci).
.compute_auc <- function(fit, B, conf, curve, slope = NULL, fpr_grid = NULL) {
  ord <- order(curve$fpr)
  AUC <- pracma::trapz(curve$fpr[ord], curve$tpr[ord])
  CI  <- .trapz_auc_ci(fit, slope, fpr_grid, B, conf)
  list(AUC = AUC, CI = CI)
}

# B parametric-bootstrap AUC draws: sample (lsens, lspec) ~ N(centre,
# vcov_fixed) and integrate each resampled SROC over the FPR grid.  Slope
# (tau_sens / tau_spec) is held fixed at the point estimate.  Returns a
# numeric vector of length B, or NULL if the draw is not possible.
.trapz_auc_draws <- function(fit, slope, fpr_grid, B) {
  if (is.null(slope) || is.null(fpr_grid) || B <= 1) return(NULL)
  f <- .fixed_se_sp(fit$fit)
  centre <- c(f$lsens, f$lspec)
  draws <- tryCatch(
    MASS::mvrnorm(n = B, mu = centre, Sigma = f$vcov_fixed),
    error = function(e) NULL
  )
  if (is.null(draws)) return(NULL)
  logit_sp <- -stats::qlogis(fpr_grid)
  ord <- order(fpr_grid)
  apply(draws, 1, function(p) {
    tpr <- stats::plogis(p[1] - slope * (logit_sp - p[2]))
    pracma::trapz(fpr_grid[ord], tpr[ord])
  })
}

# Parametric MVN bootstrap CI for the trapz AUC (quantiles of the draws).
.trapz_auc_ci <- function(fit, slope, fpr_grid, B, conf) {
  aucs <- .trapz_auc_draws(fit, slope, fpr_grid, B)
  if (is.null(aucs)) return(NULL)
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
