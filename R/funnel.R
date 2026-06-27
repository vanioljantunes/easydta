# ============================================================================
# funnel.R  -  Deeks funnel plot + Deeks test for asymmetry
#
# Cochrane Handbook for DTA Reviews v2.0 (2023), Chapter 10.6.4 / Deeks et al.
# (2005) "The performance of tests of publication bias and other sample size
# effects in systematic reviews of diagnostic test accuracy was assessed".
#
# Per-study quantities (with continuity correction `add` for cells with zeros):
#   ln(DOR_i)     = log( (TP+a)*(TN+a) / ((FP+a)*(FN+a)) )
#   var(ln(DOR))  = 1/(TP+a) + 1/(FP+a) + 1/(FN+a) + 1/(TN+a)
#   ESS_i         = 4 * n1 * n0 / (n1 + n0)        # effective sample size
#                                                 # n1 = TP+FN, n0 = FP+TN
#
# Funnel: ln(DOR) on x; 1/sqrt(ESS) on y (reversed -- large studies on top).
#         Vertical reference line = REML pooled ln(DOR).
#
# Deeks test for asymmetry:
#   weighted linear regression of ln(DOR_i) on 1/sqrt(ESS_i),
#   with weights = ESS_i.  Test of the slope coefficient = 0.
# ============================================================================

#' Deeks funnel plot + Deeks test for asymmetry (DOR scale)
#'
#' Cochrane-recommended publication-bias diagnostic for DTA meta-analysis
#' (Handbook v2.0 ch. 10.6.4; Deeks, Macaskill & Irwig 2005).  Produces a
#' funnel plot of `ln(DOR)` against `1/sqrt(ESS)` with the REML pooled
#' `ln(DOR)` as the reference line, and below it a small results table
#' from the Deeks weighted linear regression test.
#'
#' @param fit   Either a `dta_single` object or a `dta_pairwise_result`
#'   (output of `dta_compare_tests()`); when pairwise, pass `arm` to
#'   choose which test arm to plot.
#' @param arm   Required when `fit` is a `dta_pairwise_result`: name of
#'   the arm in `fit$arms` to plot.  Ignored otherwise.
#' @param test,outcome,population Title fields ("Deeks funnel of <test>
#'   to predict <outcome> in <population>").
#' @param continuity Continuity correction added to all four 2x2 cells of
#'   any study with a zero cell (default `0.5`).
#' @param conf  Confidence level for the pooled DOR CI (default 0.95).
#' @param table Logical -- show the Deeks results table below the funnel?
#'   (default `TRUE`).
#' @param ...   Reserved.
#'
#' @return A `gtable` of class `"dta_funnel"`; its `print()` method draws it
#'   on the current device, so it renders automatically when returned at the
#'   top level (no `grid::grid.draw()` needed).  The Deeks regression result
#'   is on `attr(., "deeks")` and the per-study data on `attr(., "study_data")`.
#' @examples
#' data(anti_ccp2)
#' fit <- dta_fit_single(anti_ccp2, wide = TRUE)
#' dta_funnel(fit, test = "anti-CCP2",
#'            outcome = "rheumatoid arthritis", population = "adults")
#' @export
dta_funnel <- function(fit,
                       arm = NULL,
                       test       = "test",
                       outcome    = "outcome",
                       population = "population",
                       continuity = 0.5,
                       conf  = 0.95,
                       table = TRUE,
                       ...) {
  if (!requireNamespace("metafor", quietly = TRUE)) {
    stop("Package 'metafor' is required for dta_funnel().\n",
         "  install.packages('metafor')")
  }

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

  counts <- .extract_counts(fit, NULL, "TP", "FP", "FN", "TN")
  studs  <- unique(as.character(fit$long$studlab))
  d <- data.frame(study = studs,
                  TP = counts$TP, FP = counts$FP,
                  FN = counts$FN, TN = counts$TN,
                  stringsAsFactors = FALSE)

  add0 <- ifelse(d$TP == 0 | d$FP == 0 | d$FN == 0 | d$TN == 0,
                 continuity, 0)
  TPa <- d$TP + add0; FPa <- d$FP + add0
  FNa <- d$FN + add0; TNa <- d$TN + add0
  d$lnDOR    <- log(TPa * TNa / (FPa * FNa))
  d$varlnDOR <- 1 / TPa + 1 / FPa + 1 / FNa + 1 / TNa
  d$selnDOR  <- sqrt(d$varlnDOR)

  n1 <- d$TP + d$FN
  n0 <- d$FP + d$TN
  d$ESS    <- 4 * n1 * n0 / (n1 + n0)
  d$inv_re <- 1 / sqrt(d$ESS)

  rma_fit <- metafor::rma(yi = d$lnDOR, vi = d$varlnDOR,
                          method = "REML")
  pooled_lnDOR <- as.numeric(rma_fit$b)
  z <- stats::qnorm(1 - (1 - conf) / 2)
  pooled_lnDOR_ci <- c(pooled_lnDOR - z * rma_fit$se,
                       pooled_lnDOR + z * rma_fit$se)

  deeks_fit <- stats::lm(lnDOR ~ inv_re, data = d, weights = d$ESS)
  ds <- summary(deeks_fit)
  slope <- ds$coefficients["inv_re", ]
  deeks <- list(
    slope    = unname(slope["Estimate"]),
    se       = unname(slope["Std. Error"]),
    t        = unname(slope["t value"]),
    p_value  = unname(slope["Pr(>|t|)"]),
    df       = deeks_fit$df.residual,
    n_studies = nrow(d),
    asymmetric = unname(slope["Pr(>|t|)"]) < (1 - conf)
  )

  title_text <- sprintf("Deeks funnel of %s to predict %s in %s",
                        test, outcome, population)

  # Contour-enhanced pseudo-confidence funnel.  At each y = 1/sqrt(ESS)
  # the half-width of the level-l region is z(l) * k * y, where k is the
  # data-driven scale fit (per study) selnDOR_i ~ 0 + 1/sqrt(ESS_i).
  # Three nested triangles are drawn (99 / 95 / 90 %), darker shading
  # toward the centre, so the layered "multiple-triangles" funnel shape
  # is visible behind the study points.
  k <- tryCatch(
    unname(stats::coef(stats::lm(selnDOR ~ 0 + inv_re, data = d))[1]),
    error = function(e) NA_real_
  )
  funnel_layers <- NULL
  if (is.finite(k) && k > 0) {
    y_apex <- 0
    y_base <- max(d$inv_re) * 1.05
    levs   <- c(0.99, 0.95, 0.90)
    fills  <- c("grey88", "grey78", "grey68")
    funnel_layers <- do.call(rbind, lapply(seq_along(levs), function(i) {
      z_l  <- stats::qnorm(1 - (1 - levs[i]) / 2)
      half <- z_l * k * y_base
      data.frame(
        x = c(pooled_lnDOR,
              pooled_lnDOR + half,
              pooled_lnDOR - half),
        y = c(y_apex, y_base, y_base),
        level = factor(sprintf("%d%%", round(100 * levs[i])),
                       levels = sprintf("%d%%", round(100 * levs))),
        stringsAsFactors = FALSE
      )
    }))
  }

  pooled_df <- data.frame(x = pooled_lnDOR)
  p <- ggplot2::ggplot()
  if (!is.null(funnel_layers)) {
    p <- p +
      ggplot2::geom_polygon(
        data = funnel_layers,
        ggplot2::aes(x = x, y = y, group = level, fill = level),
        colour = "black", linewidth = 0.3) +
      ggplot2::scale_fill_manual(
        name = "Pseudo-CI",
        values = stats::setNames(fills, sprintf("%d%%", round(100 * levs))))
  }
  p <- p +
    ggplot2::geom_vline(data = pooled_df,
                        ggplot2::aes(xintercept = x),
                        linetype = "dashed", colour = "black",
                        linewidth = 0.5) +
    ggplot2::geom_point(data = d,
                        ggplot2::aes(x = lnDOR, y = inv_re),
                        shape = 2, size = 2.4, colour = "black") +
    ggplot2::scale_y_reverse(
      name = expression(1 / sqrt("ESS"))) +
    ggplot2::scale_x_continuous(
      name = "ln(Diagnostic Odds Ratio)") +
    ggplot2::labs(title = title_text) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid       = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position  = "bottom",
      legend.key.size  = ggplot2::unit(0.4, "cm"),
      legend.margin    = ggplot2::margin(t = 0, b = 2, unit = "pt"),
      legend.box.spacing = ggplot2::unit(2, "pt"),
      plot.margin      = ggplot2::margin(t = 5, r = 5, b = 8, l = 5,
                                         unit = "pt"),
      plot.title       = ggplot2::element_text(face = "bold",
                                               hjust = 0.5,
                                               size = 13))

  pooled_DOR    <- exp(pooled_lnDOR)
  pooled_DOR_ci <- exp(pooled_lnDOR_ci)
  fmt_p <- function(x) if (is.na(x)) "n/a" else
                       if (x < 0.001) "<0.001" else sprintf("%.3f", x)

  deeks_df <- data.frame(
    Statistic = c("Number of studies", "Deeks p-value"),
    Value = c(sprintf("%d", deeks$n_studies),
              fmt_p(deeks$p_value)),
    stringsAsFactors = FALSE
  )

  if (isTRUE(table)) {
    tbl <- gridExtra::tableGrob(
      deeks_df, rows = NULL,
      theme = gridExtra::ttheme_default(
        base_size = 10,
        core    = list(
          fg_params = list(cex = 0.9),
          bg_params = list(fill = "white",
                           col  = "grey55", lwd = 0.6)
        ),
        colhead = list(
          fg_params = list(cex = 0.9, fontface = "bold"),
          bg_params = list(fill = "grey85",
                           col  = "grey55", lwd = 0.6)
        )
      )
    )
    # Align the table horizontally with the funnel's vertical reference
    # line (= pooled lnDOR).  We borrow the funnel's column widths from
    # its ggplotGrob -- the panel column carries a `null` unit, so when
    # the composite is drawn the table-row's panel column scales with
    # the funnel's panel column.  Inside that shared panel column we
    # place the table at the data-coord fraction frac_x = (pooled - xmin)
    # / (xmax - xmin), so the table's centre lands directly under the
    # dashed reference line at any device size.
    funnel_gt <- ggplot2::ggplotGrob(p)
    panel_layout <- funnel_gt$layout[funnel_gt$layout$name == "panel", ]
    panel_l <- panel_layout$l[1]
    panel_r <- panel_layout$r[1]
    pb <- ggplot2::ggplot_build(p)
    x_rng <- pb$layout$panel_params[[1]]$x.range
    if (is.null(x_rng)) {
      x_rng <- pb$layout$panel_scales_x[[1]]$get_limits()
    }
    frac_x <- (pooled_lnDOR - x_rng[1]) / (x_rng[2] - x_rng[1])

    tbl_aligned <- grid::grobTree(tbl, vp = grid::viewport(
      x = grid::unit(frac_x, "npc"), y = 0.5, just = "center",
      width  = grid::grobWidth(tbl),
      height = grid::grobHeight(tbl)
    ))

    # Build the table-row gtable with the funnel's column widths and
    # drop the table into the panel column.  Spacer ~20% of table-row
    # height sits between funnel and table.
    table_row <- gtable::gtable(widths  = funnel_gt$widths,
                                heights = grid::unit(1, "null"))
    table_row <- gtable::gtable_add_grob(table_row, tbl_aligned,
                                         t = 1, l = panel_l, r = panel_r,
                                         b = 1)
    spacer <- grid::nullGrob()
    g <- gridExtra::arrangeGrob(funnel_gt, spacer, table_row,
                                nrow = 3, heights = c(8, 0.2, 1))
  } else {
    g <- gridExtra::arrangeGrob(p, ncol = 1)
  }

  attr(g, "deeks")      <- deeks
  attr(g, "study_data") <- d
  attr(g, "pooled")     <- list(lnDOR = pooled_lnDOR,
                                lnDOR_CI = pooled_lnDOR_ci,
                                DOR = pooled_DOR,
                                DOR_CI = pooled_DOR_ci)
  class(g) <- c("dta_funnel", class(g))
  g
}

#' @export
print.dta_funnel <- function(x, ...) {
  grid::grid.newpage()
  grid::grid.draw(x)
  invisible(x)
}
