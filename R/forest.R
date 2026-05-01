# ============================================================================
# forest.R  -  Coupled sensitivity / specificity forest plot via grid
#
# Five ggplot panels composed side-by-side with gridExtra::grid.arrange:
#   (1) study label + raw 2x2 counts (studlab, TP, FN, TN, FP)
#   (2) Sensitivity point + exact 95% binomial CI per study, summary diamond
#   (3) Sensitivity numeric value: "est (lci-uci)" per row
#   (4) Specificity point + exact 95% binomial CI per study, summary diamond
#   (5) Specificity numeric value: "est (lci-uci)" per row
#
# Counts and study labels are reconstructed from fit$long, so the user
# passes only the fit object.  Per-study CIs are exact (Clopper-Pearson);
# the summary row uses the bivariate glmer fixed-effect Se/Sp with
# Wald CI on logit scale, then back-transformed by plogis().
#
# Layout:
#   - Each panel uses an identical y range [0.5, summary_idx + 1.5] so
#     all five panels share the same vertical layout.
#   - The header row (column titles) is rendered as in-data geom_text at
#     y = summary_idx + 1, in bold.  Using in-data labels rather than
#     plot.title means the column titles align across panels regardless
#     of axis-area heights.
#   - Non-CI panels render an invisible (colour = NA) x-axis at the
#     bottom so their inner panel-area heights match the CI panels'
#     (which carry visible axis ticks/title), and study rows line up.
#   - The label panel places the studlab and each integer count via
#     separate geom_text layers at explicit x positions, so column
#     content cannot overrun the panel width.
#
# Zebra shading alternates row backgrounds across all five panels for
# readability; the summary row is never striped.
# ============================================================================

#' Coupled sensitivity / specificity forest plot
#'
#' @param fit    Either a `dta_single` object, or a `dta_pairwise_result`
#'   (output of `dta_compare_tests()`) -- in which case `arm` selects
#'   which test arm to plot.  The 2x2 counts and the study labels are
#'   pulled from `fit$long`, so no separate data frame is needed.
#' @param arm    Required when `fit` is a `dta_pairwise_result`: name of
#'   the arm in `fit$arms` to plot (e.g. `"CCP1"`).  Ignored otherwise.
#' @param conf   Confidence level (default 0.95).
#' @param digits Numeric display digits for the value columns (default 2).
#'
#' @return A gtable object drawn on the current device.
#' @export
dta_forest <- function(fit, arm = NULL, conf = 0.95, digits = 2) {
  if (inherits(fit, "dta_pairwise_result")) {
    if (is.null(arm))
      stop("`fit` is a dta_pairwise_result; pass `arm = \"<test-name>\"` ",
           "(one of: ", paste(shQuote(names(fit$arms)), collapse = ", "), ").")
    if (!arm %in% names(fit$arms))
      stop("Arm '", arm, "' not found in fit$arms (available: ",
           paste(shQuote(names(fit$arms)), collapse = ", "), ").")
    fit <- fit$arms[[arm]]
  } else if (!is.null(arm)) {
    warning("`arm` is ignored when `fit` is a dta_single object.",
            call. = FALSE)
  }
  stopifnot(inherits(fit, "dta_single"))

  long <- fit$long
  studs <- unique(as.character(long$studlab))
  nrow_d <- length(studs)

  TP_v <- integer(nrow_d); FN_v <- integer(nrow_d)
  TN_v <- integer(nrow_d); FP_v <- integer(nrow_d)
  se_df <- data.frame(idx = seq_len(nrow_d),
                      study = studs,
                      est = NA_real_, lci = NA_real_, uci = NA_real_,
                      stringsAsFactors = FALSE)
  sp_df <- se_df

  for (i in seq_len(nrow_d)) {
    s  <- studs[i]
    rs <- long[long$studlab == s & long$sens == 1, , drop = FALSE]
    rp <- long[long$studlab == s & long$spec == 1, , drop = FALSE]
    TP <- as.integer(rs$true); n1 <- as.integer(rs$n); FN <- n1 - TP
    TN <- as.integer(rp$true); n0 <- as.integer(rp$n); FP <- n0 - TN
    TP_v[i] <- TP; FN_v[i] <- FN; TN_v[i] <- TN; FP_v[i] <- FP

    se_ci <- .exact_binom_ci(TP, n1, conf)
    sp_ci <- .exact_binom_ci(TN, n0, conf)
    se_df[i, c("est", "lci", "uci")] <- se_ci
    sp_df[i, c("est", "lci", "uci")] <- sp_ci
  }

  f <- .fixed_se_sp(fit$fit)
  se_sum <- .logit_ci(f$lsens, f$se_lsens, conf)
  sp_sum <- .logit_ci(f$lspec, f$se_lspec, conf)

  summary_idx <- nrow_d + 1
  se_df <- rbind(se_df,
                 data.frame(idx = summary_idx, study = "Summary",
                            est = se_sum["estimate"],
                            lci = se_sum["lci"],
                            uci = se_sum["uci"],
                            stringsAsFactors = FALSE))
  sp_df <- rbind(sp_df,
                 data.frame(idx = summary_idx, study = "Summary",
                            est = sp_sum["estimate"],
                            lci = sp_sum["lci"],
                            uci = sp_sum["uci"],
                            stringsAsFactors = FALSE))

  se_df$ypos <- (summary_idx + 1) - se_df$idx
  sp_df$ypos <- (summary_idx + 1) - sp_df$idx

  fmt <- function(e, l, u, d) {
    fstr <- sprintf("%%.%df (%%.%df-%%.%df)", d, d, d)
    sprintf(fstr, e, l, u)
  }
  se_df$txt <- fmt(se_df$est, se_df$lci, se_df$uci, digits)
  sp_df$txt <- fmt(sp_df$est, sp_df$lci, sp_df$uci, digits)

  # Column totals shown on the summary row.
  TP_sum <- sum(TP_v); FN_sum <- sum(FN_v)
  TN_sum <- sum(TN_v); FP_sum <- sum(FP_v)

  all_studlab <- c(studs, "Summary")
  label_df <- data.frame(
    ypos    = se_df$ypos,
    studlab = all_studlab[se_df$idx],
    TP      = ifelse(se_df$idx == summary_idx, as.character(TP_sum),
                     as.character(c(TP_v, NA)[se_df$idx])),
    FN      = ifelse(se_df$idx == summary_idx, as.character(FN_sum),
                     as.character(c(FN_v, NA)[se_df$idx])),
    TN      = ifelse(se_df$idx == summary_idx, as.character(TN_sum),
                     as.character(c(TN_v, NA)[se_df$idx])),
    FP      = ifelse(se_df$idx == summary_idx, as.character(FP_sum),
                     as.character(c(FP_v, NA)[se_df$idx])),
    is_sum  = se_df$idx == summary_idx,
    stringsAsFactors = FALSE
  )
  label_df$face <- ifelse(label_df$is_sum, "bold", "plain")
  se_df$face    <- ifelse(se_df$study == "Summary", "bold", "plain")
  sp_df$face    <- ifelse(sp_df$study == "Summary", "bold", "plain")

  # y range: a tight slot below the summary row holds the I^2 annotation
  # (ypos = 0.4); the panel bottom (ylim_full[1] = 0) is just below that
  # so the visible x-axis on the CI panels sits close to the diamond.
  header_y  <- summary_idx + 1
  i2_y      <- 0.4
  ylim_full <- c(0.05, header_y + 0.5)

  # Zebra shading: stripe every other study row; never stripe the summary
  # or the header.
  stripe_idx <- seq_len(nrow_d)
  stripe_idx <- stripe_idx[stripe_idx %% 2L == 1L]
  zebra_df <- data.frame(
    ypos = (summary_idx + 1) - stripe_idx
  )

  zebra_layer <- function() {
    ggplot2::geom_rect(
      data = zebra_df, inherit.aes = FALSE,
      ggplot2::aes(xmin = -Inf, xmax = Inf,
                   ymin = ypos - 0.5, ymax = ypos + 0.5),
      fill = "grey92"
    )
  }

  # I^2 split across three panels: biv under summary label, sens under sens
  # diamond, spec under spec diamond.
  het <- fit$heterogeneity
  if (!is.null(het)) {
    i2_biv_text  <- sprintf("I²(biv) = %.1f%%",  100 * het$I2_biv)
    i2_sens_text <- sprintf("I²(sens) = %.1f%%", 100 * het$I2_sens)
    i2_spec_text <- sprintf("I²(spec) = %.1f%%", 100 * het$I2_spec)
  } else {
    i2_biv_text <- i2_sens_text <- i2_spec_text <- "I²: NA"
  }
  i2_size <- 2.6

  # Shared geometry: every panel reserves identical bottom space for the
  # x-axis ink (line/ticks/text). Non-CI panels render that ink with NA
  # colour so the layout remains identical and zebra rows align across
  # every panel. Axis title is dropped on every panel so panel areas
  # match exactly.
  axis_text_size <- 8
  axis_line_col  <- "black"

  base_theme <- ggplot2::theme_bw() +
    ggplot2::theme(
      axis.title.y     = ggplot2::element_blank(),
      axis.title.x     = ggplot2::element_blank(),
      axis.text.y      = ggplot2::element_blank(),
      axis.ticks.y     = ggplot2::element_blank(),
      panel.grid       = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border     = ggplot2::element_blank(),
      plot.margin      = ggplot2::margin(t = 5.5, r = 0, b = 2, l = 0),
      legend.position  = "none"
    )

  invisible_axis_theme <- ggplot2::theme(
    axis.text.x  = ggplot2::element_text(colour = NA, size = axis_text_size),
    axis.ticks.x = ggplot2::element_line(colour = NA),
    axis.line.x  = ggplot2::element_line(colour = NA)
  )

  visible_axis_theme <- ggplot2::theme(
    axis.text.x  = ggplot2::element_text(colour = "black",
                                         size = axis_text_size),
    axis.ticks.x = ggplot2::element_line(colour = axis_line_col,
                                         linewidth = 0.4),
    axis.line.x  = ggplot2::element_line(colour = axis_line_col,
                                         linewidth = 0.4),
    axis.ticks.length = grid::unit(3, "pt")
  )

  # ----- Label panel: 5 columns at fixed x positions ------------------------
  x_studlab <- 0.00
  x_TP      <- 0.62
  x_FN      <- 0.74
  x_TN      <- 0.87
  x_FP      <- 1.00
  hdr_size  <- 3.2
  cell_size <- 3.0

  p_labels <- ggplot2::ggplot(label_df, ggplot2::aes(y = ypos)) +
    zebra_layer() +
    ggplot2::geom_text(ggplot2::aes(label = studlab, fontface = face),
                       x = x_studlab, hjust = 0, size = cell_size) +
    ggplot2::geom_text(ggplot2::aes(label = TP, fontface = face),
                       x = x_TP, hjust = 1, size = cell_size) +
    ggplot2::geom_text(ggplot2::aes(label = FN, fontface = face),
                       x = x_FN, hjust = 1, size = cell_size) +
    ggplot2::geom_text(ggplot2::aes(label = TN, fontface = face),
                       x = x_TN, hjust = 1, size = cell_size) +
    ggplot2::geom_text(ggplot2::aes(label = FP, fontface = face),
                       x = x_FP, hjust = 1, size = cell_size) +
    ggplot2::annotate("text", x = x_studlab, y = header_y,
                      label = "Study", hjust = 0,
                      fontface = "bold", size = hdr_size) +
    ggplot2::annotate("text", x = x_TP, y = header_y,
                      label = "TP", hjust = 1,
                      fontface = "bold", size = hdr_size) +
    ggplot2::annotate("text", x = x_FN, y = header_y,
                      label = "FN", hjust = 1,
                      fontface = "bold", size = hdr_size) +
    ggplot2::annotate("text", x = x_TN, y = header_y,
                      label = "TN", hjust = 1,
                      fontface = "bold", size = hdr_size) +
    ggplot2::annotate("text", x = x_FP, y = header_y,
                      label = "FP", hjust = 1,
                      fontface = "bold", size = hdr_size) +
    ggplot2::annotate("text", x = x_studlab, y = i2_y,
                      label = i2_biv_text, hjust = 0,
                      fontface = "bold", size = i2_size) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = ylim_full) +
    ggplot2::scale_x_continuous(breaks = seq(0, 1, 0.2)) +
    base_theme + invisible_axis_theme

  # ----- Sens / Spec value text panels --------------------------------------
  text_value_panel <- function(df, title, i2_label) {
    ggplot2::ggplot(df, ggplot2::aes(y = ypos)) +
      zebra_layer() +
      ggplot2::geom_text(ggplot2::aes(label = txt, fontface = face),
                         x = 0.05, hjust = 0, size = cell_size) +
      ggplot2::annotate("text", x = 0.05, y = header_y,
                        label = title, hjust = 0,
                        fontface = "bold", size = hdr_size) +
      ggplot2::annotate("text", x = 0.05, y = i2_y,
                        label = i2_label, hjust = 0,
                        fontface = "bold", size = i2_size) +
      ggplot2::coord_cartesian(xlim = c(0, 1), ylim = ylim_full) +
      ggplot2::scale_x_continuous(breaks = seq(0, 1, 0.2)) +
      base_theme + invisible_axis_theme
  }

  p_se_txt <- text_value_panel(se_df, "Sens (95% CI)", i2_sens_text)
  p_sp_txt <- text_value_panel(sp_df, "Spec (95% CI)", i2_spec_text)

  # ----- Sensitivity / Specificity CI panels --------------------------------
  # No box (panel.border) and no title above or below; the x-axis scale
  # (breaks/ticks/text) is preserved. A dashed reference line passes
  # through the centre of the summary diamond.
  make_panel <- function(df) {
    df$is_sum <- df$study == "Summary"
    df$point_shape <- ifelse(df$is_sum, 18L, 15L)  # diamond vs filled square
    df$point_size  <- ifelse(df$is_sum, 4, 2)
    sum_est <- df$est[df$is_sum]
    ggplot2::ggplot(df, ggplot2::aes(x = est, y = ypos)) +
      zebra_layer() +
      ggplot2::geom_vline(xintercept = sum_est,
                          linetype = "dashed", colour = "grey40",
                          linewidth = 0.4) +
      ggplot2::geom_errorbar(ggplot2::aes(xmin = lci, xmax = uci),
                             width = 0.25, orientation = "y") +
      ggplot2::geom_point(ggplot2::aes(shape = point_shape,
                                       size  = point_size)) +
      ggplot2::scale_shape_identity() +
      ggplot2::scale_size_identity() +
      ggplot2::scale_x_continuous(limits = c(0, 1),
                                  breaks = seq(0, 1, 0.2)) +
      ggplot2::coord_cartesian(ylim = ylim_full) +
      base_theme + visible_axis_theme
  }

  p_se <- make_panel(se_df)
  p_sp <- make_panel(sp_df)

  # Fixed per-y-unit height in inches so plots with different study
  # counts render with identical row spacing. Content height scales
  # linearly with the y-axis extent (ylim_full); the plot is top-aligned
  # in the device, leaving unused space blank below the axis.
  row_in        <- 0.20
  ylim_extent   <- ylim_full[2] - ylim_full[1]
  axis_pad_in   <- 0.35  # axis ticks + tick text + bottom margin
  content_h     <- grid::unit(ylim_extent * row_in, "inches") +
                   grid::unit(axis_pad_in, "inches")

  g <- gridExtra::arrangeGrob(
    p_labels, p_se, p_se_txt, p_sp, p_sp_txt,
    ncol = 5,
    widths = c(3.0, 1.8, 1.6, 1.8, 1.6),
    padding = grid::unit(0, "line")
  )

  grid::grid.newpage()
  grid::pushViewport(grid::viewport(
    y = grid::unit(1, "npc"), just = "top",
    height = content_h, width = grid::unit(1, "npc")
  ))
  grid::grid.draw(g)
  grid::popViewport()

  invisible(g)
}
