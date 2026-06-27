# Build data/schuetz.rda + inst/extdata/schuetz.xlsx
# Paired CT vs MRI for coronary artery disease (Cochrane Handbook for DTA
# Reviews v2.0 (2023), Appendix 12 schuetz.csv -- direct-comparison subset:
# the 5 studies that evaluated BOTH tests).
# Wide layout: one row per study; .e = CT (index), .c = MRI (comparator).
# Run from package root: source("data-raw/schuetz.R")

schuetz <- data.frame(
  studlab = c("Dewey 2006", "Kefer 2005", "Langer 2009",
              "Maintz 2007", "Pouleur 2008"),
  TP.e = c(62, 32, 25, 15, 16), FP.e = c(5, 6, 2, 2, 7),
  FN.e = c(4, 2, 1, 1, 1),      TN.e = c(46, 12, 40, 2, 53),
  TP.c = c(42, 30, 18, 15, 17), FP.c = c(2, 9, 15, 1, 17),
  FN.c = c(7, 4, 8, 1, 0),      TN.c = c(39, 9, 27, 3, 43),
  stringsAsFactors = FALSE
)

openxlsx::write.xlsx(list(schuetz = schuetz),
                     file = "inst/extdata/schuetz.xlsx")

usethis::use_data(schuetz, overwrite = TRUE)
