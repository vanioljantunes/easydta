# easydta: the Cochrane Appendix 12 analysis (CT vs MRI, bivariate
# meta-regression with likelihood-ratio tests) on the same 5 paired studies.
# Paste the whole file into R.
# install.packages("remotes"); remotes::install_github("vanioljantunes/easydta")

library(easydta)

# Same data as cochrane_appendix12.R, one row per study (also: data(schuetz)).
# .e = CT (index), .c = MRI (comparator)
d <- data.frame(
  studlab = c("Dewey 2006", "Kefer 2005", "Langer 2009", "Maintz 2007",
              "Pouleur 2008"),
  TP.e = c(62, 32, 25, 15, 16), FP.e = c(5, 6, 2, 2, 7),
  FN.e = c(4, 2, 1, 1, 1),      TN.e = c(46, 12, 40, 2, 53),
  TP.c = c(42, 30, 18, 15, 17), FP.c = c(2, 9, 15, 1, 17),
  FN.c = c(7, 4, 8, 1, 0),      TN.c = c(39, 9, 27, 3, 43)
)

# Equal variances: Cochrane models A, B, C, D in one call
res <- dta_pairwise(d, studlab = "studlab",
                    intervention.label = "CT", control.label = "MRI")

res$compare$lr_tests     # Cochrane: lrtest(A,B), lrtest(B,C), lrtest(B,D)
dta_summary(res$pair)    # Cochrane: sespci(B), plus DOR / LR+ / LR- per test
res$compare$differences  # Cochrane: absolute and relative differences, model B
res$arms$CT              # Cochrane: ma_CT  (+ Zhou-Dendukuri I^2)
res$arms$MRI             # Cochrane: ma_MRI (+ Zhou-Dendukuri I^2)

# Unequal variances: Cochrane model E
resE <- dta_pairwise(d, studlab = "studlab",
                     intervention.label = "CT", control.label = "MRI",
                     variance = "unequal")
dta_summary(resE$pair)    # Cochrane: sespci(E)
resE$compare$differences  # Cochrane: absolute and relative differences, model E

# Cochrane: lrtest(B,E), do the tests need separate variances?
lmtest::lrtest(res$pair$models$B, resE$pair$models$B)

# easydta adds: side-by-side SROC with AUC and the differences table
dta_sroc_pair(res, outcome = "coronary artery disease",
              population = "adults with suspected CAD",
              legend.pos = "bottomleft")
