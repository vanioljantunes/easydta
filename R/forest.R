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
# Row alignment: ALL five panels share the same x-axis space at the bottom
# (the non-CI panels render an invisible axis with `colour = NA` to keep
# the same vertical space), so the inner panel areas have identical
# heights and study rows line up exactly across all panels.
#
# Zebra shading alternates row backgrounds across all five panels for
# readability; the summary row is never striped.
# ============================================================================

#' Coupled sensitivity / specificity forest plot
#'
#' @param fit    A `dta_single` object. The 2x2 counts and the study labels
#'   are pulled from `fit$long`, so no separate data frame is needed.
#' @param conf   Confidence level (default 0.95).
#' @param digits Numeric display digits for the value columns (default 2).
#'
#' @return A gtable object drawn on the current device.
#' @export
dta_forest <- function(fit, conf = 0.95, digits = 2) {
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

  fmt <- function(e, l, u, d) sprintf(sprintf("%%.%df (%%.%df-%%.%df)", d, d, d),
                                       e, l, u)
  se_df$txt <- fmt(se_df$est, se_df$lci, se_df$uci, digits)
  sp_df$txt <- fmt(sp_df$est, sp_df$lci, sp_df$uci, digits)

  # Zebra shading: stripe every other study row; never stripe the summary.
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

  # ----- Label panel: studlab + TP + FN + TN + FP --------------------------
  fmt_int <- function(x) formatC(x, width = 4)              # right-aligned
  fmt_lab <- function(x, w = 22) {
    s <- substr(x, 1, w)
    formatC(s, width = -w, flag = "-")                       # left-aligned
  }

  study_rows <- sprintf("%s %s %s %s %s",
                        fmt_lab(studs),
                        fmt_int(TP_v),
                        fmt_int(FN_v),
                        fmt_int(TN_v),
                        fmt_int(FP_v))
  summary_row <- sprintf("%-22s %4s %4s %4s %4s",
                         "Summary", "", "", "", "")

  label_text <- ifelse(se_df$study == "Summary",
                       summary_row,
                       study_rows[match(se_df$idx, seq_len(nrow_d))])

  label_df <- data.frame(
    ypos  = se_df$ypos,
    label = label_text,
    stringsAsFactors = FALSE
  )

  header_label <- sprintf("%-22s %4s %4s %4s %4s",
                          "Study", "TP", "FN", "TN", "FP")

  # Common theme for non-CI panels: keeps the bottom axis space (so the
  # inner panel rectangles have the same height as the CI panels, which
  # have visible x-axis ticks/title/labels), but renders that area in NA
  # colour so it is invisible.
  invisible_axis <- ggplot2::theme(
    axis.text.x  = ggplot2::element_text(colour = NA),
    axis.text.y  = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_line(colour = NA),
    axis.ticks.y = ggplot2::element_blank(),
    axis.title.x = ggplot2::element_text(colour = NA),
    axis.title.y = ggplot2::element_blank(),
    panel.grid   = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(colour = NA, fill = NA),
    plot.title   = ggplot2::element_text(family = "mono",
                                         size = 9, hjust = 0)
  )

  text_panel <- function(df, lab_col, title) {
    ggplot2::ggplot(df, ggplot2::aes(x = 0, y = ypos,
                                     label = .data[[lab_col]])) +
      zebra_layer() +
      ggplot2::geom_text(hjust = 0, size = 3, family = "mono") +
      ggplot2::coord_cartesian(xlim = c(0, 1)) +
      ggplot2::scale_x_continuous(breaks = seq(0, 1, 0.2)) +
      ggplot2::theme_bw() +
      invisible_axis +
      ggplot2::labs(title = title, x = " ")
  }

  p_labels <- text_panel(label_df, "label", header_label)
  p_se_txt <- text_panel(se_df,    "txt",   "Sens (95% CI)")
  p_sp_txt <- text_panel(sp_df,    "txt",   "Spec (95% CI)")

  # ----- Sensitivity / Specificity CI panels --------------------------------
  make_panel <- function(df, title) {
    df$is_sum <- df$study == "Summary"
    df$point_shape <- ifelse(df$is_sum, 18L, 15L)  # diamond vs filled square
    df$point_size  <- ifelse(df$is_sum, 4, 2)
    ggplot2::ggplot(df, ggplot2::aes(x = est, y = ypos)) +
      zebra_layer() +
      ggplot2::geom_vline(xintercept = c(0, 0.5, 1),
                          colour = "grey80", linetype = "dotted") +
      ggplot2::geom_errorbar(ggplot2::aes(xmin = lci, xmax = uci),
                             width = 0.25, orientation = "y") +
      ggplot2::geom_point(ggplot2::aes(shape = point_shape,
                                       size  = point_size)) +
      ggplot2::scale_shape_identity() +
      ggplot2::scale_size_identity() +
      ggplot2::scale_x_continuous(limits = c(0, 1),
                                  breaks = seq(0, 1, 0.2)) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "none",
                     axis.text.y     = ggplot2::element_blank(),
                     axis.ticks.y    = ggplot2::element_blank(),
                     axis.title.y    = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank(),
                     plot.title = ggplot2::element_text(size = 9, hjust = 0)) +
      ggplot2::labs(x = title, title = title)
  }

  p_se <- make_panel(se_df, "Sensitivity (95% CI)")
  p_sp <- make_panel(sp_df, "Specificity (95% CI)")

  gridExtra::grid.arrange(
    p_labels, p_se, p_se_txt, p_sp, p_sp_txt,
    ncol = 5,
    widths = c(2.4, 1.8, 1.5, 1.8, 1.5)
  )
}
