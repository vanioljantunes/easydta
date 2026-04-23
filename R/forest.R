# ============================================================================
# forest.R  -  Coupled sensitivity / specificity forest plot via grid
#
# Three ggplot panels composed side-by-side with gridExtra::grid.arrange:
#   (1) study label + raw counts (TP/FN and FP/TN)
#   (2) Sensitivity point + exact 95% binomial CI per study, summary diamond
#   (3) Specificity point + exact 95% binomial CI per study, summary diamond
#
# Summary (diamond) row uses the bivariate glmer fixed-effect Se and Sp with
# Wald CI on logit scale then back-transformed by plogis().
# ============================================================================

#' Coupled sensitivity / specificity forest plot
#'
#' @param fit    A `dta_single` object.
#' @param data   The original WIDE data frame (one row per study) carrying
#'   `Study_ID`, `TP`, `FP`, `FN`, `TN`.
#' @param tp,fp,fn,tn,study  Column names.
#' @param conf   Confidence level (default 0.95).
#' @param digits Numeric display digits (default 2).
#'
#' @return A gtable object drawn on the current device.
#' @export
dta_forest <- function(fit,
                       data,
                       tp    = "TP",
                       fp    = "FP",
                       fn    = "FN",
                       tn    = "TN",
                       studlab = "studlab",
                       conf  = 0.95,
                       digits = 2) {
  stopifnot(inherits(fit, "dta_single"))
  req <- c(tp, fp, fn, tn, studlab)
  miss <- setdiff(req, names(data))
  if (length(miss)) stop("Missing columns in data: ",
                         paste(miss, collapse = ", "))

  nrow_d <- nrow(data)
  se_df <- data.frame(idx = seq_len(nrow_d),
                      study = as.character(data[[studlab]]),
                      est = NA_real_, lci = NA_real_, uci = NA_real_)
  sp_df <- se_df
  for (i in seq_len(nrow_d)) {
    TP <- data[[tp]][i]; FN <- data[[fn]][i]
    FP <- data[[fp]][i]; TN <- data[[tn]][i]
    se_ci <- .exact_binom_ci(TP, TP + FN, conf)
    sp_ci <- .exact_binom_ci(TN, TN + FP, conf)
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
                            uci = se_sum["uci"]))
  sp_df <- rbind(sp_df,
                 data.frame(idx = summary_idx, study = "Summary",
                            est = sp_sum["estimate"],
                            lci = sp_sum["lci"],
                            uci = sp_sum["uci"]))

  se_df$ypos <- (summary_idx + 1) - se_df$idx
  sp_df$ypos <- (summary_idx + 1) - sp_df$idx

  label_df <- data.frame(
    ypos = se_df$ypos,
    label = ifelse(se_df$study == "Summary",
                   "Summary",
                   paste0(se_df$study, "   ",
                          data[[tp]][se_df$idx], "/",
                          data[[tp]][se_df$idx] + data[[fn]][se_df$idx],
                          "   ",
                          data[[fp]][se_df$idx], "/",
                          data[[fp]][se_df$idx] + data[[tn]][se_df$idx]))
  )
  label_df$label[is.na(label_df$label)] <- "Summary"

  p_labels <- ggplot2::ggplot(label_df, ggplot2::aes(x = 0, y = ypos,
                                                     label = label)) +
    ggplot2::geom_text(hjust = 0, size = 3, family = "mono") +
    ggplot2::coord_cartesian(xlim = c(0, 1)) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text  = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   axis.title = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank()) +
    ggplot2::labs(title = "Study   TP/n1   FP/n0")

  make_panel <- function(df, title) {
    df$is_sum <- df$study == "Summary"
    df$point_shape <- ifelse(df$is_sum, 18L, 15L)  # diamond vs filled square
    df$point_size  <- ifelse(df$is_sum, 4, 2)
    ggplot2::ggplot(df, ggplot2::aes(x = est, y = ypos)) +
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
                     axis.text.y = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank(),
                     axis.title.y = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank()) +
      ggplot2::labs(x = title, title = title)
  }

  p_se <- make_panel(se_df, "Sensitivity (95% CI)")
  p_sp <- make_panel(sp_df, "Specificity (95% CI)")

  gridExtra::grid.arrange(
    p_labels, p_se, p_sp,
    ncol = 3, widths = c(2, 2, 2)
  )
}
