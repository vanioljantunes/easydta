# ============================================================================
# roc_points.R  -  Raw study points on ROC space (sanity-check scatter)
# ============================================================================

#' Raw study-level ROC scatter
#'
#' @param data  Wide data frame (one row per study) with TP/FP/FN/TN columns.
#' @param tp,fp,fn,tn Column names.
#' @param title Plot title.
#'
#' @return A ggplot object.
#' @examples
#' data(anti_ccp2)
#' dta_roc_points(anti_ccp2)
#' @export
dta_roc_points <- function(data,
                           tp = "TP",
                           fp = "FP",
                           fn = "FN",
                           tn = "TN",
                           title = "Study-level ROC scatter") {
  stopifnot(is.data.frame(data))
  req <- c(tp, fp, fn, tn)
  miss <- setdiff(req, names(data))
  if (length(miss)) stop("Missing columns in data: ",
                         paste(miss, collapse = ", "))

  n1 <- as.integer(data[[tp]]) + as.integer(data[[fn]])
  n0 <- as.integer(data[[fp]]) + as.integer(data[[tn]])
  df <- data.frame(
    fpr = as.integer(data[[fp]]) / pmax(n0, 1L),
    tpr = as.integer(data[[tp]]) / pmax(n1, 1L),
    n_total = n1 + n0
  )
  ggplot2::ggplot(df, ggplot2::aes(x = fpr, y = tpr, size = n_total)) +
    ggplot2::geom_point(shape = 21, fill = "grey70", alpha = 0.7) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::scale_x_continuous(breaks = seq(0, 1, 0.2)) +
    ggplot2::scale_y_continuous(breaks = seq(0, 1, 0.2)) +
    ggplot2::labs(x = "False positive rate (1 - Specificity)",
                  y = "Sensitivity",
                  title = title,
                  size  = "Total n") +
    ggplot2::theme_bw()
}
