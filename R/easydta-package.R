#' easydta: Cochrane-compliant diagnostic test accuracy meta-analysis
#'
#' Pedagogical, \pkg{lme4}-based workflow for bivariate binomial GLMM
#' meta-analysis of diagnostic test accuracy (DTA) studies, following the
#' Cochrane Handbook for Systematic Reviews of Diagnostic Test Accuracy
#' (v2.0, 2023), Chapter 10 and Appendix 12 (Takwoingi et al.).
#'
#' @section Guides (read in order):
#' \enumerate{
#'   \item [easydta-1-data-preparation] -- shaping data and the bundled datasets.
#'   \item [easydta-2-single-arm] -- one-test meta-analysis.
#'   \item [easydta-3-pairwise] -- two-test comparison.
#'   \item [easydta-4-forests] -- coupled forest plots.
#'   \item [easydta-5-funnels] -- Deeks funnel + asymmetry test.
#' }
#'
#' @section Single-test workflow:
#' \describe{
#'   \item{[dta_fit_single()]}{Fit the bivariate GLMM for one test.}
#'   \item{[dta_summary()], [dta_derived()]}{Se/Sp, DOR, LR+/LR- with CIs.}
#'   \item{[dta_forest()], [dta_sroc()], [dta_funnel()],
#'     [dta_roc_points()]}{Coupled forest, SROC, Deeks funnel, raw ROC scatter.}
#' }
#'
#' @section Pairwise / comparison workflow:
#' \describe{
#'   \item{[dta_pairwise()]}{One call for a paired (wide `.e`/`.c`) design.}
#'   \item{[dta_compare_tests()]}{One call for a between-study covariate design.}
#'   \item{[dta_fit_pairwise()], [dta_compare()]}{Lower-level nested-model fit
#'     (models A/B/C/D, equal/unequal variance) and the LR tests + Se/Sp
#'     differences.}
#'   \item{[dta_sroc_pair()]}{Side-by-side SROC with a Cochrane differences
#'     table.}
#' }
#'
#' @section Reshaping helpers:
#' \describe{
#'   \item{[dta_reshape()], [dta_reshape_pairwise()]}{Wide -> long; usually
#'     called internally.}
#' }
#'
#' @section Example data:
#' \describe{
#'   \item{[anti_ccp1], [anti_ccp2]}{Single-arm anti-CCP subsets (CCP1, CCP2).}
#'   \item{[schuetz]}{Paired CT vs MRI (direct comparison).}
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
#' @keywords internal
"_PACKAGE"
