# ============================================================================
# data.R  -  One help page for every bundled dataset.
#
# The topic name carries "-package", which R's index writer hoists to the top
# of the package Help Pages list, so the datasets are the first thing a reader
# meets.  `?anti_ccp1`, `?anti_ccp2` and `?schuetz` all land here.
# ============================================================================

#' Datasets: the three example datasets shipped with easydta
#'
#' Three datasets from the Cochrane Handbook for DTA Reviews, one per study
#' design.  `anti_ccp1` and `anti_ccp2` are single-arm (one test per study),
#' `schuetz` is paired (every study ran both tests).  Which one you imitate
#' decides which function you call.
#'
#' @section Which dataset matches which design:
#' \tabular{llll}{
#'   \strong{Dataset} \tab \strong{Design} \tab \strong{Shape} \tab \strong{Fit with} \cr
#'   `anti_ccp1` \tab one test, 8 studies \tab single-arm wide \tab [dta_fit_single()] \cr
#'   `anti_ccp2` \tab one test, 29 studies \tab single-arm wide \tab [dta_fit_single()] \cr
#'   `rbind(anti_ccp1, anti_ccp2)` \tab two tests, different studies \tab stacked, `test` column \tab [dta_compare_tests()] \cr
#'   `schuetz` \tab two tests, same 5 studies \tab paired wide, `.e` / `.c` \tab [dta_pairwise()]
#' }
#'
#' @section anti_ccp1, anti_ccp2:
#' Anti-cyclic citrullinated peptide (anti-CCP) antibody for rheumatoid
#' arthritis, split by assay generation: `anti_ccp1` is the first generation
#' (CCP1, 8 studies), `anti_ccp2` the second (CCP2, 29 studies).  Each study
#' used exactly one generation, so the two subsets share no study and a
#' comparison between them is a between-study (indirect) covariate comparison,
#' never a paired one.
#'
#' One row per study, six columns:
#' \describe{
#'   \item{studlab}{character. Study label (first author and year).}
#'   \item{TP}{integer. True positives.}
#'   \item{FP}{integer. False positives.}
#'   \item{FN}{integer. False negatives.}
#'   \item{TN}{integer. True negatives.}
#'   \item{test}{factor, `"CCP1"` / `"CCP2"`.  Constant inside each subset;
#'     it is there so the two can be row-bound into one covariate frame.}
#' }
#'
#' Also shipped as `inst/extdata/anti_ccp.xlsx`, one sheet per generation.
#'
#' @section schuetz:
#' Computed tomography (CT) against magnetic resonance imaging (MRI) for
#' coronary artery disease.  Five studies, each of which ran \emph{both}
#' tests, so the two arms are correlated within study and the comparison is
#' paired (direct).
#'
#' One row per study, nine columns:
#' \describe{
#'   \item{studlab}{character. Study label (first author and year).}
#'   \item{TP.e, FP.e, FN.e, TN.e}{integer. 2x2 counts for CT, the index arm.}
#'   \item{TP.c, FP.c, FN.c, TN.c}{integer. 2x2 counts for MRI, the comparator.}
#' }
#'
#' Also shipped as `inst/extdata/schuetz.xlsx`.
#'
#' @source Anti-CCP: Nishimura K, Sugiyama D, Kogata Y, et al. Meta-analysis:
#'   diagnostic accuracy of anti-cyclic citrullinated peptide antibody and
#'   rheumatoid factor for rheumatoid arthritis. Ann Intern Med.
#'   2007;146(11):797-808, as reproduced in the Cochrane Handbook for DTA
#'   Reviews v2.0 (2023), Chapter 10 Supplementary Material 1.
#'
#'   schuetz: Cochrane Handbook for DTA Reviews v2.0 (2023), Appendix 12
#'   (`schuetz.csv`), direct-comparison subset.
#'
#' @format See the two sections above: `anti_ccp1` (8 x 6), `anti_ccp2`
#'   (29 x 6), `schuetz` (5 x 9).
#'
#' @examples
#' data(anti_ccp2)
#' head(anti_ccp2)
#'
#' data(schuetz)
#' head(schuetz)
#'
#' @seealso The reading order starts at [easydta-1-data-preparation].
#' @docType data
#' @keywords datasets
#' @name easydta-package-datasets
#' @aliases anti_ccp anti_ccp1 anti_ccp2 schuetz
NULL
