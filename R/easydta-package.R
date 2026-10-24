#' easydta: Cochrane-compliant diagnostic test accuracy meta-analysis
#'
#' Bivariate binomial GLMM meta-analysis of diagnostic test accuracy (DTA)
#' studies on top of \pkg{lme4}, following the Cochrane Handbook for
#' Systematic Reviews of Diagnostic Test Accuracy (v2.0, 2023), Chapter 10
#' and Appendix 12 (Takwoingi et al.).  One call per deliverable, no other
#' meta-analysis dependency.
#'
#' @section Start here:
#' \tabular{ll}{
#'   [easydta-1-data-preparation] \tab Import an Excel sheet, shape it, check it \cr
#'   [easydta-2-single-arm] \tab One test: fit, read the summary, plot it \cr
#'   [easydta-3-pairwise] \tab Two tests: paired or between-study, and the variance choice \cr
#'   [easydta-4-forests] \tab Coupled sensitivity / specificity forests \cr
#'   [easydta-5-funnels] \tab Deeks funnel and the asymmetry test \cr
#'   [easydta-package-datasets] \tab The three bundled datasets
#' }
#'
#' @section What you have, and what to call:
#' \tabular{ll}{
#'   \strong{You have} \tab \strong{Call} \cr
#'   One test, one row per study \tab [dta_fit_single()] \cr
#'   Two tests, the same studies ran both \tab [dta_pairwise()] \cr
#'   Two tests, each study ran one of them \tab [dta_compare_tests()] \cr
#'   Counts already stacked in long format \tab [dta_fit_pairwise()] after [dta_reshape()]
#' }
#'
#' @section One test, step by step:
#' \tabular{ll}{
#'   [dta_fit_single()] \tab Fit the bivariate binomial GLMM \cr
#'   [dta_summary()] \tab Pooled sensitivity and specificity with CIs \cr
#'   [dta_derived()] \tab DOR, LR+ and LR- with delta-method CIs \cr
#'   [dta_forest()] \tab Coupled forest of both measures \cr
#'   [dta_sroc()] \tab SROC curve, confidence and prediction regions, AUC \cr
#'   [dta_roc_points()] \tab Raw study points, no model \cr
#'   [dta_funnel()] \tab Deeks funnel plus asymmetry test \cr
#'   [dta_loo()], [dta_loo_plot()] \tab Leave-one-out table and sROC panel
#' }
#'
#' @section Two tests, step by step:
#' \tabular{ll}{
#'   [dta_pairwise()] \tab Whole pipeline for a paired `.e` / `.c` frame \cr
#'   [dta_compare_tests()] \tab Whole pipeline for a between-study covariate frame \cr
#'   [dta_fit_pairwise()] \tab The nested models A / B / C / D on their own \cr
#'   [dta_compare()] \tab Likelihood-ratio tests and Se / Sp differences \cr
#'   [dta_sroc_pair()] \tab Both arms side by side with a differences table \cr
#'   [dta_forest_pair()] \tab Forest of the difference between the two tests
#' }
#'
#' @section Plumbing you rarely call yourself:
#' \tabular{ll}{
#'   [dta_reshape()] \tab Single-arm wide to long \cr
#'   [dta_reshape_pairwise()] \tab Paired wide to long
#' }
#'
#' @examples
#' data(anti_ccp2)
#' fit <- dta_fit_single(anti_ccp2, wide = TRUE)
#' dta_summary(fit)
#'
#' data(schuetz)
#' res <- dta_pairwise(schuetz, studlab = "studlab",
#'                     intervention.label = "CT", control.label = "MRI")
#' print(res)
#'
"_PACKAGE"
