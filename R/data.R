#' Anti-CCP antibody for diagnosing rheumatoid arthritis (single-arm subsets)
#'
#' Diagnostic test accuracy data for anti-cyclic citrullinated peptide (anti-CCP)
#' antibody testing for rheumatoid arthritis, split by assay generation into two
#' single-arm datasets: `anti_ccp1` (first-generation, CCP1) and `anti_ccp2`
#' (second-generation, CCP2).  Each study used exactly one assay generation, so
#' the two subsets are disjoint at the study level (a between-study design).
#'
#' @format Data frames with one row per study and 6 variables:
#' \describe{
#'   \item{studlab}{character. Study label (first author and year).}
#'   \item{TP}{integer. True positives.}
#'   \item{FP}{integer. False positives.}
#'   \item{FN}{integer. False negatives.}
#'   \item{TN}{integer. True negatives.}
#'   \item{test}{factor with levels \code{"CCP1"} and \code{"CCP2"} (constant
#'     within each subset; kept so the two can be row-bound into a covariate
#'     comparison via [dta_compare_tests()]).}
#' }
#' `anti_ccp1` has 8 studies; `anti_ccp2` has 29.
#' @details Distributed as `inst/extdata/anti_ccp.xlsx`, one sheet per arm
#'   (`CCP1`, `CCP2`).  Each subset suits a single-arm fit via
#'   [dta_fit_single()]; `rbind(anti_ccp1, anti_ccp2)` reconstructs the full
#'   covariate dataset.
#' @source Nishimura K, Sugiyama D, Kogata Y, et al. Meta-analysis: diagnostic
#'   accuracy of anti-cyclic citrullinated peptide antibody and rheumatoid
#'   factor for rheumatoid arthritis. Ann Intern Med. 2007;146(11):797-808.
#'   Reproduced via the Cochrane Handbook for DTA Reviews v2.0 (2023),
#'   Chapter 10 Supplementary Material 1.
#' @examples
#' data(anti_ccp2)
#' head(anti_ccp2)
#' @name anti_ccp
"anti_ccp1"

#' @rdname anti_ccp
#' @format NULL
"anti_ccp2"

#' CT vs MRI for coronary artery disease (paired, direct comparison)
#'
#' Diagnostic test accuracy data comparing computed tomography (CT) and magnetic
#' resonance imaging (MRI) for coronary artery disease.  Five studies, each of
#' which evaluated \emph{both} tests -- a paired within-study (direct)
#' comparison.  Stored wide: one row per study with `.e` (CT) and `.c` (MRI)
#' suffixed 2x2 counts.
#'
#' @format A data frame with 5 rows (one per study) and 9 variables:
#' \describe{
#'   \item{studlab}{character. Study label (first author and year).}
#'   \item{TP.e, FP.e, FN.e, TN.e}{integer. 2x2 counts for CT (index arm).}
#'   \item{TP.c, FP.c, FN.c, TN.c}{integer. 2x2 counts for MRI (comparator).}
#' }
#' @details Distributed as `inst/extdata/schuetz.xlsx` (single sheet,
#'   `schuetz`).  Use with [dta_pairwise()], naming the arms via
#'   `intervention.label = "CT"`, `control.label = "MRI"`.
#' @source Cochrane Handbook for DTA Reviews v2.0 (2023), Appendix 12
#'   (`schuetz.csv`), direct-comparison subset (studies with both tests).
#' @examples
#' data(schuetz)
#' head(schuetz)
"schuetz"
