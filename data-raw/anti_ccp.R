# Build data/anti_ccp1.rda + data/anti_ccp2.rda from
# inst/extdata/anti_ccp.xlsx (one sheet per arm: "CCP1", "CCP2").
# Each is a single-arm subset; the `test` column is kept (constant) so the
# two can be row-bound into a covariate comparison.
# Run from package root: source("data-raw/anti_ccp.R")

path <- system.file("extdata", "anti_ccp.xlsx", package = "easydta")
if (!nzchar(path)) {
  path <- "inst/extdata/anti_ccp.xlsx"
}

lev <- c("CCP1", "CCP2")
anti_ccp1 <- as.data.frame(readxl::read_excel(path, sheet = "CCP1"))
anti_ccp2 <- as.data.frame(readxl::read_excel(path, sheet = "CCP2"))
anti_ccp1$test <- factor(anti_ccp1$test, levels = lev)
anti_ccp2$test <- factor(anti_ccp2$test, levels = lev)

usethis::use_data(anti_ccp1, overwrite = TRUE)
usethis::use_data(anti_ccp2, overwrite = TRUE)
