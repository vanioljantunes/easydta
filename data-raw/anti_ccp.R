# Build data/anti_ccp.rda from inst/extdata/anti_ccp.csv
# Run from package root: source("data-raw/anti_ccp.R")

anti_ccp <- utils::read.csv(
  system.file("extdata", "anti_ccp.csv", package = "easydta"),
  stringsAsFactors = FALSE
)

if (nrow(anti_ccp) == 0L) {
  anti_ccp <- utils::read.csv("inst/extdata/anti_ccp.csv", stringsAsFactors = FALSE)
}

anti_ccp$generation <- factor(anti_ccp$generation, levels = c("CCP1", "CCP2"))

usethis::use_data(anti_ccp, overwrite = TRUE)
