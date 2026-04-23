# ============================================================================
# reshape.R  -  Wide 2x2 -> long sens/spec format for glmer
#
# Two reshape entry points:
#
#   dta_reshape()           -- single-test input: one row per study with
#                              columns TP, FP, FN, TN plus a `studlab` col.
#
#   dta_reshape_pairwise()  -- paired-design input: one row per study with
#                              columns author, year, TP.e/FP.e/FN.e/TN.e
#                              (intervention) and TP.c/FP.c/FN.c/TN.c
#                              (control), plus user-supplied intervention
#                              and control names.  Returns the same long
#                              format as dta_reshape() with an added test
#                              column carrying the two test labels.
#
# Both end in the canonical Cochrane long layout:
#   studlab | [test] | sens (0/1) | spec (0/1) | true | n
#
# Cochrane ref: Appendix 5, p.12 (single); Appendix 12, p.55-56 (pairwise).
# ============================================================================

#' Reshape wide 2x2 data (single-test) into long sens/spec format
#'
#' @param data    Data frame with one row per study.
#' @param tp,fp,fn,tn  Column names (strings) for TP/FP/FN/TN counts.
#' @param studlab Column name of the study label.
#' @param extra   Character vector of additional columns to carry through.
#'
#' @return A long-format data frame with columns: studlab, sens, spec,
#'   true, n, plus any `extra` columns.
#' @export
dta_reshape <- function(data,
                        tp      = "TP",
                        fp      = "FP",
                        fn      = "FN",
                        tn      = "TN",
                        studlab = "studlab",
                        extra   = NULL) {
  stopifnot(is.data.frame(data))
  req <- c(tp, fp, fn, tn, studlab)
  miss <- setdiff(req, names(data))
  if (length(miss)) {
    stop("Missing required columns: ", paste(miss, collapse = ", "))
  }

  X <- data
  X$.n1     <- as.integer(X[[tp]]) + as.integer(X[[fn]])
  X$.n0     <- as.integer(X[[fp]]) + as.integer(X[[tn]])
  X$.true1  <- as.integer(X[[tp]])
  X$.true0  <- as.integer(X[[tn]])
  X$.rec_id <- seq_len(nrow(X))

  keep_cols <- c(studlab, extra)

  Y <- stats::reshape(
    X,
    direction = "long",
    varying   = list(c(".n1",    ".n0"),
                     c(".true1", ".true0")),
    timevar   = "sens",
    times     = c(1L, 0L),
    v.names   = c("n", "true"),
    idvar     = ".rec_id"
  )

  Y$spec     <- 1L - Y$sens
  Y$studlab  <- Y[[studlab]]
  Y <- Y[order(Y$.rec_id, -Y$sens), , drop = FALSE]
  rownames(Y) <- NULL

  out_cols <- c("studlab", setdiff(keep_cols, "studlab"),
                "sens", "spec", "true", "n")
  out_cols <- out_cols[!duplicated(out_cols)]
  Y <- Y[, out_cols, drop = FALSE]
  Y
}

#' Reshape paired-design wide data (intervention vs control) into long form
#'
#' Accepts one row per study with both arms of a paired comparison side by
#' side (the `.e` suffix for the experimental / intervention arm, `.c` for
#' the control arm) and stacks them into the Cochrane long format with
#' a `test` column carrying the two arm labels.
#'
#' @param data          Data frame with one row per study.
#' @param author        Column name of the first-author label.
#' @param year          Column name of the publication year.
#' @param intervention  Label to assign to the `.e` rows (e.g. "CT").
#' @param control       Label to assign to the `.c` rows (e.g. "MRI").
#' @param tp.e,fp.e,fn.e,tn.e  Experimental-arm 2x2 column names.
#' @param tp.c,fp.c,fn.c,tn.c  Control-arm 2x2 column names.
#' @param studlab       Optional: column to use as the study label instead
#'   of `paste(author, year)`.  If NULL (default), the studlab is built as
#'   `paste(data[[author]], data[[year]])`.
#' @param test_var      Name for the stacked test-type column (default
#'   `"test"`).  Pass this same name to `dta_fit_pairwise(test_var = ...)`.
#'
#' @return Long-format data frame with columns: studlab, <test_var>,
#'   sens, spec, true, n.
#' @export
dta_reshape_pairwise <- function(data,
                                 author       = "author",
                                 year         = "year",
                                 intervention = "Intervention",
                                 control      = "Control",
                                 tp.e         = "TP.e",
                                 fp.e         = "FP.e",
                                 fn.e         = "FN.e",
                                 tn.e         = "TN.e",
                                 tp.c         = "TP.c",
                                 fp.c         = "FP.c",
                                 fn.c         = "FN.c",
                                 tn.c         = "TN.c",
                                 studlab      = NULL,
                                 test_var     = "test") {
  stopifnot(is.data.frame(data))
  need <- c(tp.e, fp.e, fn.e, tn.e, tp.c, fp.c, fn.c, tn.c)
  if (is.null(studlab)) need <- c(need, author, year) else need <- c(need, studlab)
  miss <- setdiff(need, names(data))
  if (length(miss)) {
    stop("Missing columns: ", paste(miss, collapse = ", "))
  }

  studlab_vec <- if (is.null(studlab)) {
    paste(data[[author]], data[[year]])
  } else {
    as.character(data[[studlab]])
  }

  # Stack the two arms into one wide-per-arm frame, then call dta_reshape().
  E <- data.frame(
    studlab = studlab_vec,
    test    = intervention,
    TP      = as.integer(data[[tp.e]]),
    FP      = as.integer(data[[fp.e]]),
    FN      = as.integer(data[[fn.e]]),
    TN      = as.integer(data[[tn.e]]),
    stringsAsFactors = FALSE
  )
  Cc <- data.frame(
    studlab = studlab_vec,
    test    = control,
    TP      = as.integer(data[[tp.c]]),
    FP      = as.integer(data[[fp.c]]),
    FN      = as.integer(data[[fn.c]]),
    TN      = as.integer(data[[tn.c]]),
    stringsAsFactors = FALSE
  )
  stacked <- rbind(E, Cc)
  # Rename `test` to the user-requested column name, if different
  if (test_var != "test") {
    stacked[[test_var]] <- stacked$test
    stacked$test <- NULL
  }

  dta_reshape(stacked,
              tp = "TP", fp = "FP", fn = "FN", tn = "TN",
              studlab = "studlab",
              extra   = test_var)
}
