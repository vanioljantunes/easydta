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
#'   (from `dta_pairwise()` / `dta_compare_tests()`) -- in which case `test`
#'   selects which arm to plot.  The 2x2 counts and the study labels are
#'   pulled from `fit$long`, so no separate data frame is needed.
#' @param test   Required when `fit` is a `dta_pairwise_result`: which arm to
#'   plot, given as the role `intervention` or `control` (a bare word or a
#'   string); the arm name is looked up from the result.  Ignored otherwise.
#' @param conf   Confidence level (default 0.95).
#' @param digits Numeric display digits for the value columns (default 2).
#' @param just   Horizontal justification of the numeric value columns
#'   (Sens/Spec "est (lci-uci)" and their headers): `"center"` (default),
#'   `"left"`, or `"right"`. Centring places each value in the middle of
#'   its own column.
#' @param title  Optional plot title, rendered left-aligned above the plot.
#' @param legend Optional legend text shown left-aligned below the plot.
#'   A character vector is rendered one line per element (similar to an
#'   "add row" of free text).
#'
#' @param subgroup Optional study-level grouping. Either a named vector
#'   (names = study labels, values = subgroup) or an unnamed vector in the
#'   order of the studies in `fit`. Studies are listed under a heading per
#'   subgroup, in order of first appearance, each subgroup with at least two
#'   studies closed by its own pooled bivariate estimate (italic diamond row),
#'   and the overall estimate last.
#' @param test.subgroup Logical. With exactly two subgroups of at least two
#'   studies each, print the likelihood-ratio test for a difference in
#'   sensitivity and specificity between them (bivariate meta-regression,
#'   via [dta_compare_tests()]) under the I-squared row. Default `TRUE`.
#'
#' @return A gtable object drawn on the current device.
#' @examples
#' data(anti_ccp2)
#' fit <- dta_fit_single(anti_ccp2, wide = TRUE)
#' dta_forest(fit)
#' @export
dta_forest <- function(fit, test, conf = 0.95, digits = 2,
                       just = c("center", "left", "right"),
                       title = NULL, legend = NULL,
                       subgroup = NULL, test.subgroup = TRUE) {
  test_sym <- substitute(test)
  just     <- match.arg(just)
  val_x    <- switch(just, center = 0.5, left = 0.05, right = 0.95)
  val_hjust<- switch(just, center = 0.5, left = 0,    right = 1)
  if (inherits(fit, "dta_pairwise_result")) {
    if (missing(test))
      stop("`fit` is a dta_pairwise_result; pass `test = intervention` or ",
           "`test = control`.")
    fit <- fit$arms[[.role_to_arm(.role_from(test_sym, test), fit)]]
  } else if (!missing(test)) {
    warning("`test` is ignored when `fit` is a dta_single object.",
            call. = FALSE)
  }
  stopifnot(inherits(fit, "dta_single"))

  long   <- fit$long
  studs  <- unique(as.character(long$studlab))
  nrow_d <- length(studs)

  cnt <- data.frame(studlab = studs, TP = NA_integer_, FN = NA_integer_,
                    TN = NA_integer_, FP = NA_integer_,
                    stringsAsFactors = FALSE)
  for (i in seq_len(nrow_d)) {
    rs <- long[long$studlab == studs[i] & long$sens == 1, , drop = FALSE]
    rp <- long[long$studlab == studs[i] & long$spec == 1, , drop = FALSE]
    cnt$TP[i] <- as.integer(rs$true); cnt$FN[i] <- as.integer(rs$n) - cnt$TP[i]
    cnt$TN[i] <- as.integer(rp$true); cnt$FP[i] <- as.integer(rp$n) - cnt$TN[i]
  }

  # Row table, top to bottom.  type: "group" (subgroup heading), "study",
  # "subsum" (pooled estimate of one subgroup), "sum" (overall pooled).
  mk_row <- function(type, label, js, se, sp) {
    data.frame(type = type, label = label,
               TP = if (length(js)) sum(cnt$TP[js]) else NA_integer_,
               FN = if (length(js)) sum(cnt$FN[js]) else NA_integer_,
               TN = if (length(js)) sum(cnt$TN[js]) else NA_integer_,
               FP = if (length(js)) sum(cnt$FP[js]) else NA_integer_,
               se_est = se[1], se_lci = se[2], se_uci = se[3],
               sp_est = sp[1], sp_lci = sp[2], sp_uci = sp[3],
               stringsAsFactors = FALSE, row.names = NULL)
  }
  study_rows <- function(js) {
    do.call(rbind, lapply(js, function(j) mk_row(
      "study", cnt$studlab[j], j,
      .exact_binom_ci(cnt$TP[j], cnt$TP[j] + cnt$FN[j], conf),
      .exact_binom_ci(cnt$TN[j], cnt$TN[j] + cnt$FP[j], conf))))
  }
  pooled <- function(f) {
    list(se = .logit_ci(f$lsens, f$se_lsens, conf),
         sp = .logit_ci(f$lspec, f$se_lspec, conf))
  }
  na3 <- c(NA_real_, NA_real_, NA_real_)

  sub_p <- NULL
  if (is.null(subgroup)) {
    rows <- study_rows(seq_len(nrow_d))
  } else {
    grp  <- .forest_subgroup(subgroup, studs)
    levs <- unique(grp)
    rows <- NULL
    for (lv in levs) {
      js   <- which(grp == lv)
      rows <- rbind(rows, mk_row("group", lv, integer(0), na3, na3),
                    study_rows(js))
      # a subgroup summary needs at least two studies for the bivariate GLMM
      if (length(js) >= 2) {
        sfit <- dta_fit_single(long[long$studlab %in% studs[js], , drop = FALSE],
                               wide = FALSE, conf = conf)
        ps   <- pooled(.fixed_se_sp(sfit$fit))
        rows <- rbind(rows, mk_row("subsum",
                                   sprintf("Subtotal (k = %d)", length(js)),
                                   js, ps$se, ps$sp))
      }
    }
    # likelihood-ratio test (joint Se and Sp) between two subgroups, each
    # with at least two studies; Cochrane Handbook DTA ch. 10 meta-regression
    k_lev <- table(grp)[levs]
    if (isTRUE(test.subgroup) && length(levs) == 2 && all(k_lev >= 2)) {
      wd <- cbind(cnt, subgroup = grp)
      sub_p <- tryCatch(
        dta_compare_tests(wd, test_var = "subgroup",
                          conf = conf)$compare$lr_tests$p_value[1],
        error = function(e) NA_real_)
    }
  }
  po   <- pooled(.fixed_se_sp(fit$fit))
  rows <- rbind(rows, mk_row("sum", if (is.null(subgroup)) "Summary" else "Overall",
                             seq_len(nrow_d), po$se, po$sp))

  n_rows    <- nrow(rows)
  rows$ypos <- (n_rows + 1) - seq_len(n_rows)

  fmt <- function(e, l, u, d) {
    fstr <- sprintf("%%.%df (%%.%df-%%.%df)", d, d, d)
    ifelse(is.na(e), "", sprintf(fstr, e, l, u))
  }
  face <- ifelse(rows$type %in% c("group", "sum"), "bold",
                 ifelse(rows$type == "subsum", "italic", "plain"))

  se_df <- data.frame(ypos = rows$ypos, type = rows$type,
                      est = rows$se_est, lci = rows$se_lci, uci = rows$se_uci,
                      txt = fmt(rows$se_est, rows$se_lci, rows$se_uci, digits),
                      face = face, stringsAsFactors = FALSE)
  sp_df <- data.frame(ypos = rows$ypos, type = rows$type,
                      est = rows$sp_est, lci = rows$sp_lci, uci = rows$sp_uci,
                      txt = fmt(rows$sp_est, rows$sp_lci, rows$sp_uci, digits),
                      face = face, stringsAsFactors = FALSE)

  cnt_chr <- function(x) ifelse(is.na(x), "", as.character(x))
  label_df <- data.frame(ypos = rows$ypos, studlab = rows$label,
                         TP = cnt_chr(rows$TP), FN = cnt_chr(rows$FN),
                         TN = cnt_chr(rows$TN), FP = cnt_chr(rows$FP),
                         face = face, stringsAsFactors = FALSE)

  # y range: the I^2 annotation gets its own full row one unit below the
  # summary (ypos = 0), the subgroup test one more row below (ypos = -1);
  # the panel bottom sits half a row below the last text row so the
  # visible x-axis on the CI panels clears it.
  header_y  <- n_rows + 1
  i2_y      <- 0
  test_y    <- -1
  ylim_full <- c(if (is.null(sub_p)) -0.5 else -1.5, header_y + 0.5)

  # Zebra shading: stripe every other study row; never stripe summaries,
  # subgroup headings or the header.
  study_y  <- rows$ypos[rows$type == "study"]
  zebra_df <- data.frame(ypos = study_y[seq_along(study_y) %% 2L == 1L])

  zebra_layer <- function() {
    ggplot2::geom_rect(
      data = zebra_df, inherit.aes = FALSE,
      ggplot2::aes(xmin = -Inf, xmax = Inf,
                   ymin = ypos - 0.5, ymax = ypos + 0.5),
      fill = "grey92"
    )
  }

  # Single bivariate I^2 annotation, shown under the summary label only.
  het <- fit$heterogeneity
  if (!is.null(het)) {
    i2_biv_text <- sprintf("I² = %.1f%%", 100 * het$I2_biv)
  } else {
    i2_biv_text <- "I²: NA"
  }
  i2_size <- 3.0  # standard cell font, own row (matches cell_size)

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
                      size = i2_size) +
    (if (!is.null(sub_p))
      ggplot2::annotate("text", x = x_studlab, y = test_y,
                        label = if (is.na(sub_p))
                          "Test for subgroup differences: not estimable"
                        else sprintf("Test for subgroup differences: p %s",
                                     if (sub_p < 0.001) "< 0.001"
                                     else sprintf("= %.3f", sub_p)),
                        hjust = 0, size = i2_size)) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = ylim_full) +
    ggplot2::scale_x_continuous(breaks = seq(0, 1, 0.2)) +
    base_theme + invisible_axis_theme

  # ----- Sens / Spec value text panels --------------------------------------
  # No per-panel I^2 annotation; the single bivariate I^2 lives under the
  # summary label in the label panel.
  text_value_panel <- function(df, col_title) {
    ggplot2::ggplot(df, ggplot2::aes(y = ypos)) +
      zebra_layer() +
      ggplot2::geom_text(ggplot2::aes(label = txt, fontface = face),
                         x = val_x, hjust = val_hjust, size = cell_size) +
      ggplot2::annotate("text", x = val_x, y = header_y,
                        label = col_title, hjust = val_hjust,
                        fontface = "bold", size = hdr_size) +
      ggplot2::coord_cartesian(xlim = c(0, 1), ylim = ylim_full) +
      ggplot2::scale_x_continuous(breaks = seq(0, 1, 0.2)) +
      base_theme + invisible_axis_theme
  }

  p_se_txt <- text_value_panel(se_df, "Sens (95% CI)")
  p_sp_txt <- text_value_panel(sp_df, "Spec (95% CI)")

  # ----- Sensitivity / Specificity CI panels --------------------------------
  # No box (panel.border) and no title above or below; the x-axis scale
  # (breaks/ticks/text) is preserved. A dashed reference line passes
  # through the centre of the summary diamond.
  make_panel <- function(df) {
    sum_est <- df$est[df$type == "sum"]
    df <- df[!is.na(df$est), , drop = FALSE]
    df$point_shape <- ifelse(df$type == "study", 15L, 18L)  # square vs diamond
    df$point_size  <- ifelse(df$type == "sum", 4,
                             ifelse(df$type == "subsum", 3.2, 2))
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

  main <- gridExtra::arrangeGrob(
    p_labels, p_se, p_se_txt, p_sp, p_sp_txt,
    ncol = 5,
    widths = c(3.0, 1.8, 1.6, 1.8, 1.6),
    padding = grid::unit(0, "line")
  )

  # Stack optional left-aligned title (above) and legend (below) with the
  # main five-panel grob in a single column. Each band gets a fixed height
  # in points so the overall viewport height stays deterministic.
  grobs   <- list(main)
  heights <- list(content_h)

  if (!is.null(title)) {
    title_fontsize <- 12
    title_grob <- grid::textGrob(
      title, x = grid::unit(2, "pt"), hjust = 0,
      gp = grid::gpar(fontface = "bold", fontsize = title_fontsize)
    )
    grobs   <- c(list(title_grob), grobs)
    heights <- c(list(grid::unit(1.8 * title_fontsize, "points")), heights)
  }

  if (!is.null(legend)) {
    leg_txt      <- paste(legend, collapse = "\n")
    n_lines      <- length(strsplit(leg_txt, "\n", fixed = TRUE)[[1]])
    leg_fontsize <- 9
    legend_grob  <- grid::textGrob(
      leg_txt, x = grid::unit(2, "pt"), y = grid::unit(1, "npc"),
      hjust = 0, vjust = 1,
      gp = grid::gpar(fontsize = leg_fontsize, lineheight = 1.2)
    )
    grobs   <- c(grobs, list(legend_grob))
    heights <- c(heights,
                 list(grid::unit(1.3 * leg_fontsize * n_lines + 6, "points")))
  }

  heights_vec <- do.call(grid::unit.c, heights)
  total_h     <- Reduce(`+`, heights)

  g <- gridExtra::arrangeGrob(grobs = grobs, ncol = 1, heights = heights_vec)

  grid::grid.newpage()
  # Centre the fixed-height composite in the device, both axes.
  grid::pushViewport(grid::viewport(
    x = grid::unit(0.5, "npc"), y = grid::unit(0.5, "npc"),
    just = "centre",
    height = total_h, width = grid::unit(1, "npc")
  ))
  grid::grid.draw(g)
  grid::popViewport()

  invisible(g)
}

# Resolve the `subgroup` argument of dta_forest() to one value per study.
.forest_subgroup <- function(subgroup, studs) {
  if (!is.null(names(subgroup))) {
    miss <- setdiff(studs, names(subgroup))
    if (length(miss))
      stop("`subgroup` has no value for: ", paste(miss, collapse = ", "))
    return(as.character(subgroup[studs]))
  }
  if (length(subgroup) != length(studs))
    stop("unnamed `subgroup` must have one value per study (", length(studs), ").")
  as.character(subgroup)
}
