# ============================================================================
# guides.R  -  Documentation-only "guide" help pages (no exported objects).
#
# Each block below generates a standalone Rd topic. The numeric prefix in the
# @name keeps them grouped and ordered in the package Help Pages index:
#   1 data preparation -> 2 single-arm -> 3 pairwise -> 4 forests -> 5 funnels
# ============================================================================

#' easydta guide 1: Data preparation
#'
#' How to shape diagnostic test accuracy (DTA) data for `easydta` and which
#' bundled datasets to use.
#'
#' @details
#' `easydta` accepts two wide layouts, one row per study:
#'
#' \strong{Single-arm wide} -- one test, many studies. Columns
#' `studlab, TP, FP, FN, TN` (plus optional covariates). Feed to
#' [dta_fit_single()] (with `wide = TRUE`) or stack two arms and feed to
#' [dta_compare_tests()] as a between-study covariate.
#'
#' \strong{Paired wide} -- two tests in the \emph{same} studies. One row per
#' study with `.e` (index) and `.c` (comparator) suffixed counts
#' (`TP.e, FP.e, ..., TP.c, FP.c, ...`) and a shared `studlab`. Feed to
#' [dta_pairwise()].
#'
#' Bundled datasets: [anti_ccp1] / [anti_ccp2] (single-arm anti-CCP subsets)
#' and [schuetz] (paired CT vs MRI). [dta_reshape()] and
#' [dta_reshape_pairwise()] convert wide to the long format the fitters use,
#' but the one-call wrappers reshape internally.
#'
#' @examples
#' data(anti_ccp2)   # single-arm: studlab, TP, FP, FN, TN, test
#' head(anti_ccp2)
#'
#' data(schuetz)     # paired wide: studlab + .e (CT) / .c (MRI) counts
#' head(schuetz)
#'
#' @seealso [dta_reshape()], [dta_reshape_pairwise()],
#'   the next guide [easydta-2-single-arm].
#' @name easydta-1-data-preparation
NULL

#' easydta guide 2: Single-arm analysis
#'
#' Meta-analyse one diagnostic test across many studies with the bivariate
#' binomial GLMM.
#'
#' @details
#' [dta_fit_single()] fits the Cochrane Appendix-5 model
#' `cbind(true, n - true) ~ 0 + sens + spec + (0 + sens + spec | studlab)`
#' and attaches heterogeneity (\eqn{\tau}, \eqn{\rho}, Zhou-Dendukuri bivariate
#' \eqn{I^2}, prediction-region ellipse). Summaries come from [dta_summary()]
#' and [dta_derived()] (Se, Sp, DOR, LR+/LR- with CIs). Plot it with
#' [dta_forest()] (guide 4), [dta_sroc()], and [dta_funnel()] (guide 5).
#'
#' @examples
#' data(anti_ccp2)
#' fit <- dta_fit_single(anti_ccp2, wide = TRUE)
#' print(fit)
#' dta_summary(fit)
#'
#' @seealso [dta_fit_single()], [dta_summary()], [dta_derived()], [dta_sroc()],
#'   the next guide [easydta-3-pairwise].
#' @name easydta-2-single-arm
NULL

#' easydta guide 3: Pairwise analysis
#'
#' Compare two diagnostic tests head-to-head (Cochrane Handbook Appendix 12).
#'
#' @details
#' Two entry points return the same `dta_pairwise_result`:
#' \itemize{
#'   \item [dta_pairwise()] for a \strong{paired} wide `.e`/`.c` frame (each
#'     study did both tests, e.g. [schuetz]).
#'   \item [dta_compare_tests()] for a \strong{between-study covariate} frame
#'     (each study did one test, e.g. `rbind(anti_ccp1, anti_ccp2)`).
#' }
#' Both fit the nested models A/B/C/D, run likelihood-ratio tests (overall,
#' Sens-differs, Spec-differs) and delta-method Se/Sp differences. `variance =
#' "equal"` is Cochrane model B (default); `"unequal"` is model E (separate
#' between-study variances per test). Visualise with [dta_sroc_pair()] and
#' per-arm [dta_forest()] / [dta_funnel()].
#'
#' @examples
#' data(schuetz)
#' res <- dta_pairwise(schuetz, studlab = "studlab",
#'                     intervention = "CT", control = "MRI")
#' print(res)
#'
#' @seealso [dta_pairwise()], [dta_compare_tests()], [dta_sroc_pair()],
#'   the next guide [easydta-4-forests].
#' @name easydta-3-pairwise
NULL

#' easydta guide 4: Forest plots
#'
#' Coupled sensitivity / specificity forest plots.
#'
#' @details
#' [dta_forest()] draws a single composite of aligned panels (study labels,
#' 2x2 counts, sensitivity forest, specificity forest, numeric estimates) and
#' draws itself on the current device. For a `dta_single` fit, pass it
#' directly. For a pairwise result, select an arm with `arm = "<test>"`
#' (e.g. `arm = "CT"`).
#'
#' @examples
#' data(anti_ccp2)
#' fit <- dta_fit_single(anti_ccp2, wide = TRUE)
#' dta_forest(fit)
#'
#' \donttest{
#' data(schuetz)
#' res <- dta_pairwise(schuetz, studlab = "studlab",
#'                     intervention = "CT", control = "MRI")
#' dta_forest(res, arm = "CT")
#' }
#'
#' @seealso [dta_forest()], the next guide [easydta-5-funnels].
#' @name easydta-4-forests
NULL

#' easydta guide 5: Funnel plots
#'
#' Deeks funnel plot and asymmetry test (publication-bias diagnostic).
#'
#' @details
#' [dta_funnel()] (Cochrane Handbook v2.0 sec. 10.6.4; Deeks, Macaskill &
#' Irwig 2005) plots `ln(DOR)` against `1 / sqrt(ESS)` (y-axis reversed so
#' large studies sit at the top) with contour pseudo-CI bands, and reports the
#' weighted-regression Deeks asymmetry test beneath. A slope p-value < 0.05
#' suggests asymmetry consistent with publication bias. The returned object
#' draws itself when printed. For a pairwise result, select an arm with `arm`.
#'
#' @examples
#' data(anti_ccp2)
#' fit <- dta_fit_single(anti_ccp2, wide = TRUE)
#' dta_funnel(fit, test = "anti-CCP2",
#'            outcome = "rheumatoid arthritis", population = "adults")
#'
#' @seealso [dta_funnel()], the package overview [easydta-package].
#' @name easydta-5-funnels
NULL
