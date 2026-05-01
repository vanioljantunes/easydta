#' Anti-CCP antibody for diagnosing rheumatoid arthritis
#'
#' Diagnostic test accuracy data for anti-cyclic citrullinated peptide (anti-CCP)
#' antibody testing for rheumatoid arthritis, including both first-generation
#' (CCP1) and second-generation (CCP2) assays. A canonical example dataset used
#' throughout the Cochrane Handbook for DTA Reviews.
#'
#' @format A data frame with 37 rows and 6 variables:
#' \describe{
#'   \item{studlab}{character. Study label (first author and year).}
#'   \item{TP}{integer. True positives.}
#'   \item{FP}{integer. False positives.}
#'   \item{FN}{integer. False negatives.}
#'   \item{TN}{integer. True negatives.}
#'   \item{generation}{factor with levels \code{"CCP1"} and \code{"CCP2"}.
#'     Assay generation, used as the covariate in pairwise meta-regression.}
#' }
#' @source Nishimura K, Sugiyama D, Kogata Y, et al. Meta-analysis: diagnostic
#'   accuracy of anti-cyclic citrullinated peptide antibody and rheumatoid
#'   factor for rheumatoid arthritis. Ann Intern Med. 2007;146(11):797-808.
#'   Reproduced via the Cochrane Handbook for DTA Reviews v2.0 (2023),
#'   Chapter 10 Supplementary Material 1.
#' @examples
#' data(anti_ccp)
#' head(anti_ccp)
#' table(anti_ccp$generation)
"anti_ccp"
