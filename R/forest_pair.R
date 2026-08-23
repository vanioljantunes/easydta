# Paired forest plot: per-study difference between two tests, for sensitivity
# and/or specificity.  Layout mirrors dta_forest(): a label column, one text
# column per arm, the CI panel, and the difference text column.

# Wilson score interval for a single proportion.
.wilson_ci <- function(k, n, conf) {
  if (is.na(k) || is.na(n) || n <= 0) return(c(NA_real_, NA_real_))
  z   <- stats::qnorm(1 - (1 - conf) / 2)
  p   <- k / n
  d   <- 1 + z^2 / n
  ctr <- (p + z^2 / (2 * n)) / d
  hw  <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d
  c(max(0, ctr - hw), min(1, ctr + hw))
}

# Newcombe (1998) hybrid-score interval for the difference of two INDEPENDENT
# proportions (his method 10).  The two arms of a paired DTA review are not
# independent, but the discordance counts needed for a paired interval are not
# recoverable from two marginal 2x2 tables, so this is used per study and is
# mildly conservative.  The summary row uses the model-based paired difference.
.newcombe_diff <- function(k1, n1, k2, n2, conf) {
  if (any(is.na(c(k1, n1, k2, n2))) || n1 <= 0 || n2 <= 0)
    return(c(estimate = NA_real_, lci = NA_real_, uci = NA_real_))
  p1 <- k1 / n1
  p2 <- k2 / n2
  w1 <- .wilson_ci(k1, n1, conf)
  w2 <- .wilson_ci(k2, n2, conf)
  d  <- p1 - p2
  c(estimate = d,
    lci = max(-1, d - sqrt((p1 - w1[1])^2 + (w2[2] - p2)^2)),
    uci = min( 1, d + sqrt((w1[2] - p1)^2 + (p2 - w2[1])^2)))
}

# Per-study numerator / denominator for one measure from a dta_single fit.
.arm_counts <- function(fit, measure) {
  long  <- fit$long
  studs <- unique(as.character(long$studlab))
  keep  <- if (measure == "sens") long$sens == 1 else long$spec == 1
  r     <- long[keep, , drop = FALSE]
  data.frame(studlab = studs,
             k = as.integer(r$true[match(studs, as.character(r$studlab))]),
             n = as.integer(r$n[match(studs, as.character(r$studlab))]),
             stringsAsFactors = FALSE)
}

#' Paired forest plot of the difference between two tests
#'
#' For every study that contributes both tests, plots the difference in
#' sensitivity (and/or specificity) between the intervention arm (`.e`) and
#' the control arm (`.c`), with each arm's own estimate shown as a text
#' column.  By default sensitivity and specificity are drawn as two blocks in
#' a single figure, sensitivity on top.
#'
#' Per-study differences carry a Newcombe hybrid-score interval computed as if
#' the two arms were independent samples; the discordance counts needed for a
#' properly paired interval are not recoverable from two marginal 2x2 tables,
#' so these intervals are mildly conservative.  The summary row is the
#' model-based paired difference from the bivariate meta-regression
#' (`x$compare$differences`) with the likelihood-ratio p-value from
#' `x$compare$lr_tests`, exactly as reported by `dta_sroc_pair()`.
#'
#' @param x       A `dta_pairwise_result` (from `dta_pairwise()` or
#'   `dta_compare_tests()`).
#' @param arm.e,arm.c Arm names; default to the intervention / control roles
#'   recorded on `x`.
#' @param test.label.e,test.label.c Column headings for the two arms; default
#'   to the arm names.
#' @param which   Which measures to draw, in order: `"sens"`, `"spec"`, or
#'   both (the default, sensitivity above specificity).
#' @param conf    Confidence level (default 0.95).
#' @param digits  Display digits for the value columns (default 2).
#' @param just    Horizontal justification of the numeric value columns
#'   (default `"left"`, which keeps each column tight against the one before it).
#' @param counts Draw the two count columns of each arm (`TP` and `FN` in the
#'   sensitivity block, `TN` and `FP` in the specificity block). Default `TRUE`.
#' @param widths Relative widths of the eleven columns: study labels, then for
#'   each arm its two counts, its estimate and its confidence interval, then the
#'   difference plot and the difference values. Defaults to
#'   `c(1.6, 0.4, 0.4, 0.5, 1.15, 0.4, 0.4, 0.5, 1.15, 1.8, 1.4)`. When
#'   `counts = FALSE` the four count columns are dropped.
#' @param title   Optional plot title, rendered left-aligned above the blocks.
#' @param legend  Optional legend text shown left-aligned below the blocks.
#'   When `NULL` a one-line note on the interval method is written for you;
#'   pass `NA` to suppress it.
#'
#' @return A gtable object, drawn on the current device.
#' @examples
#' \donttest{
#' data(schuetz)
#' res <- dta_pairwise(schuetz, studlab = "studlab",
#'                     intervention.label = "CT", control.label = "MRI")
#' dta_forest_pair(res)
#' }
#' @export
dta_forest_pair <- function(x,
                            arm.e = NULL, arm.c = NULL,
                            test.label.e = NULL, test.label.c = NULL,
                            which = c("sens", "spec"),
                            conf = 0.95, digits = 2,
                            just = c("left", "center", "right"),
                            counts = TRUE,
                            widths = c(1.6, 0.4, 0.4, 0.5, 1.15,
                                       0.4, 0.4, 0.5, 1.15, 1.8, 1.4),
                            title = NULL, legend = NULL) {

  if (!inherits(x, c("dta_pairwise_result", "dta_compare")))
    stop("`x` must be a dta_pairwise_result (dta_pairwise / dta_compare_tests).")

  which     <- match.arg(which, c("sens", "spec"), several.ok = TRUE)
  just      <- match.arg(just)
  if (length(widths) != 11 || !is.numeric(widths) || any(!is.finite(widths)) ||
      any(widths <= 0))
    stop("`widths` must be eleven positive numbers.")
  keep <- if (counts) seq_len(11) else c(1, 4, 5, 8, 9, 10, 11)
  val_x     <- switch(just, center = 0.5, left = 0.02, right = 0.98)
  val_hjust <- switch(just, center = 0.5, left = 0,    right = 1)

  if (is.null(arm.e)) arm.e <- x$labels$intervention
  if (is.null(arm.c)) arm.c <- x$labels$control
  if (is.null(arm.e) || is.null(arm.c))
    stop("Provide `arm.e` and `arm.c` (no intervention/control roles found).")
  if (is.null(test.label.e)) test.label.e <- arm.e
  if (is.null(test.label.c)) test.label.c <- arm.c

  miss <- setdiff(c(arm.e, arm.c), names(x$arms))
  if (length(miss))
    stop("Arms not found in x$arms: ", paste(miss, collapse = ", "))

  fit_e <- x$arms[[arm.e]]
  fit_c <- x$arms[[arm.c]]

  studs <- intersect(unique(as.character(fit_e$long$studlab)),
                     unique(as.character(fit_c$long$studlab)))
  if (!length(studs))
    stop("No study contributes both arms; a paired forest needs paired data.")

  fmt <- function(e, l, u) {
    f <- sprintf("%%.%df (%%.%df to %%.%df)", digits, digits, digits)
    ifelse(is.na(e), "", sprintf(f, e, l, u))
  }
  fmt_est <- function(e) ifelse(is.na(e), "", sprintf(sprintf("%%.%df", digits), e))
  fmt_ci  <- function(l, u) {
    f <- sprintf("(%%.%df to %%.%df)", digits, digits)
    ifelse(is.na(l), "", sprintf(f, l, u))
  }
  fmt_p <- function(p) {
    if (length(p) == 0 || is.na(p)) "p = NA"
    else if (p < 0.001) "p < 0.001"
    else sprintf("p = %.3f", p)
  }

  measure_title <- c(sens = "Sensitivity", spec = "Specificity")
  measure_short <- c(sens = "Sens", spec = "Spec")
  count_titles  <- list(sens = c("TP", "FN"), spec = c("TN", "FP"))
  diff_row      <- c(sens = "Absolute diff Sens", spec = "Absolute diff Spec")
  lr_row        <- c(sens = 2L, spec = 3L)

  axis_text_size <- 8
  hdr_size  <- 3.2
  cell_size <- 3.0
  row_in      <- 0.20
  axis_pad_in <- 0.35

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
    axis.text.x  = ggplot2::element_text(colour = "black", size = axis_text_size),
    axis.ticks.x = ggplot2::element_line(colour = "black", linewidth = 0.4),
    axis.line.x  = ggplot2::element_line(colour = "black", linewidth = 0.4),
    axis.ticks.length = grid::unit(3, "pt")
  )

  blocks  <- list()
  heights <- list()

  for (m in which) {
    ce <- .arm_counts(fit_e, m)
    cc <- .arm_counts(fit_c, m)
    ce <- ce[match(studs, ce$studlab), ]
    cc <- cc[match(studs, cc$studlab), ]

    nrow_d      <- length(studs)
    summary_idx <- nrow_d + 1
    header_y    <- summary_idx + 1
    note_y      <- 0
    ylim_full   <- c(-0.5, header_y + 0.5)

    d <- data.frame(idx = seq_len(nrow_d), study = studs,
                    est = NA_real_, lci = NA_real_, uci = NA_real_,
                    est_e = "", ci_e = "", est_c = "", ci_c = "",
                    n1_e = "", n2_e = "", n1_c = "", n2_c = "",
                    stringsAsFactors = FALSE)

    for (i in seq_len(nrow_d)) {
      ci_e <- .wilson_ci(ce$k[i], ce$n[i], conf)
      ci_c <- .wilson_ci(cc$k[i], cc$n[i], conf)
      # k is TP in the sensitivity block and TN in the specificity block, so
      # n - k is the complementary count (FN or FP).
      d$n1_e[i] <- format(ce$k[i]);          d$n2_e[i] <- format(ce$n[i] - ce$k[i])
      d$n1_c[i] <- format(cc$k[i]);          d$n2_c[i] <- format(cc$n[i] - cc$k[i])
      d$est_e[i] <- fmt_est(ce$k[i] / ce$n[i]); d$ci_e[i] <- fmt_ci(ci_e[1], ci_e[2])
      d$est_c[i] <- fmt_est(cc$k[i] / cc$n[i]); d$ci_c[i] <- fmt_ci(ci_c[1], ci_c[2])
      dd <- .newcombe_diff(ce$k[i], ce$n[i], cc$k[i], cc$n[i], conf)
      d[i, c("est", "lci", "uci")] <- dd
    }

    # Summary row: pooled per-arm estimates and the model-based difference.
    fe <- .fixed_se_sp(fit_e$fit)
    fc <- .fixed_se_sp(fit_c$fit)
    pe <- if (m == "sens") .logit_ci(fe$lsens, fe$se_lsens, conf) else
                           .logit_ci(fe$lspec, fe$se_lspec, conf)
    pc <- if (m == "sens") .logit_ci(fc$lsens, fc$se_lsens, conf) else
                           .logit_ci(fc$lspec, fc$se_lspec, conf)

    dif <- x$compare$differences
    dr  <- dif[dif$measure == diff_row[[m]], ][1, ]
    d <- rbind(d, data.frame(idx = summary_idx, study = "Summary",
                             est = dr$estimate, lci = dr$lci, uci = dr$uci,
                             est_e = fmt_est(pe["estimate"]),
                             ci_e  = fmt_ci(pe["lci"], pe["uci"]),
                             est_c = fmt_est(pc["estimate"]),
                             ci_c  = fmt_ci(pc["lci"], pc["uci"]),
                             n1_e = "", n2_e = "", n1_c = "", n2_c = "",
                             stringsAsFactors = FALSE))

    d$ypos    <- (summary_idx + 1) - d$idx
    d$is_sum  <- d$study == "Summary"
    d$face    <- ifelse(d$is_sum, "bold", "plain")
    d$txt_d   <- fmt(d$est, d$lci, d$uci)
    d$studlab <- d$study

    p_note <- fmt_p(x$compare$lr_tests$p_value[lr_row[[m]]])

    stripe_idx <- seq_len(nrow_d)
    zebra_df <- data.frame(ypos = (summary_idx + 1) -
                                  stripe_idx[stripe_idx %% 2L == 1L])
    zebra <- ggplot2::geom_rect(
      data = zebra_df, inherit.aes = FALSE,
      ggplot2::aes(xmin = -Inf, xmax = Inf,
                   ymin = ypos - 0.5, ymax = ypos + 0.5),
      fill = "grey92")

    p_labels <- ggplot2::ggplot(d, ggplot2::aes(y = ypos)) +
      zebra +
      ggplot2::geom_text(ggplot2::aes(label = studlab, fontface = face),
                         x = 0, hjust = 0, size = cell_size) +
      ggplot2::annotate("text", x = 0, y = header_y, label = "Study",
                        hjust = 0, fontface = "bold", size = hdr_size) +
      ggplot2::annotate("text", x = 0, y = note_y, label = p_note,
                        hjust = 0, size = cell_size) +
      ggplot2::coord_cartesian(xlim = c(0, 1), ylim = ylim_full) +
      base_theme + invisible_axis_theme

    text_panel <- function(lab_col, col_title, x = val_x, hjust = val_hjust) {
      dd <- d
      dd$val <- dd[[lab_col]]
      ggplot2::ggplot(dd, ggplot2::aes(y = ypos)) +
        zebra +
        ggplot2::geom_text(ggplot2::aes(label = val, fontface = face),
                           x = x, hjust = hjust, size = cell_size) +
        ggplot2::annotate("text", x = x, y = header_y, label = col_title,
                          hjust = hjust, fontface = "bold",
                          size = hdr_size) +
        ggplot2::coord_cartesian(xlim = c(0, 1), ylim = ylim_full) +
        base_theme + invisible_axis_theme
    }

    lim  <- max(abs(c(d$lci, d$uci)), na.rm = TRUE)
    lim  <- min(1, ceiling(lim * 10) / 10)
    brks <- pretty(c(-lim, lim), n = 5)

    d$point_shape <- ifelse(d$is_sum, 18L, 15L)
    d$point_size  <- ifelse(d$is_sum, 4, 2)

    p_diff <- ggplot2::ggplot(d, ggplot2::aes(x = est, y = ypos)) +
      zebra +
      ggplot2::geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.4) +
      ggplot2::geom_errorbar(ggplot2::aes(xmin = lci, xmax = uci),
                             width = 0.25, orientation = "y") +
      ggplot2::geom_point(ggplot2::aes(shape = point_shape, size = point_size)) +
      ggplot2::scale_shape_identity() +
      ggplot2::scale_size_identity() +
      ggplot2::scale_x_continuous(limits = c(-lim, lim),
                                  breaks = brks[brks >= -lim & brks <= lim]) +
      ggplot2::coord_cartesian(ylim = ylim_full) +
      base_theme + visible_axis_theme

    ct  <- count_titles[[m]]
    est_title <- measure_short[[m]]

    # Single-line headers: a wrapped header overlaps the first study row,
    # because row height is fixed in inches while the font is not.
    cols <- list(
      p_labels,
      text_panel("n1_e", ct[1], x = 0.5, hjust = 0.5),
      text_panel("n2_e", ct[2], x = 0.5, hjust = 0.5),
      text_panel("est_e", est_title, x = 0.5, hjust = 0.5),
      text_panel("ci_e", paste0(100 * conf, "% CI")),
      text_panel("n1_c", ct[1], x = 0.5, hjust = 0.5),
      text_panel("n2_c", ct[2], x = 0.5, hjust = 0.5),
      text_panel("est_c", est_title, x = 0.5, hjust = 0.5),
      text_panel("ci_c", paste0(100 * conf, "% CI")),
      p_diff,
      text_panel("txt_d", "Diff (95% CI)"))

    w <- widths[keep]
    body <- gridExtra::arrangeGrob(
      grobs = cols[keep], ncol = length(keep),
      widths = w,
      padding = grid::unit(0, "line"))

    # Arm labels span their own block of columns, as in the intervention and
    # control headers of a meta forest plot.
    span_fontsize <- 10
    span_e <- if (counts) sum(widths[2:5]) else sum(widths[c(4, 5)])
    span_c <- if (counts) sum(widths[6:9]) else sum(widths[c(8, 9)])
    arm_label <- function(lab)
      grid::textGrob(lab, gp = grid::gpar(fontface = "bold",
                                          fontsize = span_fontsize))
    spanner <- gridExtra::arrangeGrob(
      grobs = list(grid::nullGrob(), arm_label(test.label.e),
                   arm_label(test.label.c), grid::nullGrob(),
                   grid::nullGrob()),
      ncol = 5,
      widths = c(widths[1], span_e, span_c, widths[10], widths[11]))
    span_h <- grid::unit(1.5 * span_fontsize, "points")

    m_fontsize <- 11
    m_grob <- grid::textGrob(measure_title[[m]],
                             x = grid::unit(2, "pt"), hjust = 0,
                             gp = grid::gpar(fontface = "bold",
                                             fontsize = m_fontsize))
    m_h <- grid::unit(1.8 * m_fontsize, "points")
    body_h <- grid::unit((ylim_full[2] - ylim_full[1]) * row_in, "inches") +
              grid::unit(axis_pad_in, "inches")

    blocks  <- c(blocks,  list(m_grob, spanner, body))
    heights <- c(heights, list(m_h,    span_h,  body_h))
  }

  if (is.null(legend))
    legend <- paste("Study differences use Newcombe hybrid-score intervals",
                    "computed as if the arms were independent;",
                    "the summary row is the model-based paired difference.")
  if (length(legend) == 1 && is.na(legend)) legend <- NULL

  grobs <- blocks
  if (!is.null(title)) {
    title_fontsize <- 12
    grobs   <- c(list(grid::textGrob(title, x = grid::unit(0.5, "npc"),
                                     hjust = 0.5,
                                     gp = grid::gpar(fontface = "bold",
                                                     fontsize = title_fontsize))),
                 grobs)
    heights <- c(list(grid::unit(1.8 * title_fontsize, "points")), heights)
  }
  if (!is.null(legend)) {
    leg_txt      <- paste(legend, collapse = "\n")
    n_lines      <- length(strsplit(leg_txt, "\n", fixed = TRUE)[[1]])
    leg_fontsize <- 9
    grobs   <- c(grobs, list(grid::textGrob(leg_txt, x = grid::unit(2, "pt"),
                                            y = grid::unit(1, "npc"),
                                            hjust = 0, vjust = 1,
                                            gp = grid::gpar(fontsize = leg_fontsize,
                                                            lineheight = 1.2))))
    heights <- c(heights,
                 list(grid::unit(1.3 * leg_fontsize * n_lines + 6, "points")))
  }

  heights_vec <- do.call(grid::unit.c, heights)
  total_h     <- Reduce(`+`, heights)
  g <- gridExtra::arrangeGrob(grobs = grobs, ncol = 1, heights = heights_vec)

  grid::grid.newpage()
  grid::pushViewport(grid::viewport(
    x = grid::unit(0.5, "npc"), y = grid::unit(0.5, "npc"), just = "centre",
    height = total_h, width = grid::unit(1, "npc")))
  grid::grid.draw(g)
  grid::popViewport()

  invisible(g)
}
