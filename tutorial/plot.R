# Renders the example plots shown in the tutorial image.
# Run from the package root: Rscript tutorial/plot.R

suppressMessages(devtools::load_all(quiet = TRUE))
dir.create("man/figures", showWarnings = FALSE, recursive = TRUE)

# Single-arm: anti-CCP2, SROC with confidence + prediction regions and AUC
data(anti_ccp2)
fit <- dta_fit_single(anti_ccp2, tp = "TP", fp = "FP", fn = "FN", tn = "TN",
                      studlab = "studlab")
png("man/figures/tutorial-sroc.png", width = 1200, height = 1100, res = 170)
print(dta_sroc(fit, test.label = "anti-CCP2",
               outcome = "rheumatoid arthritis", population = "adults",
               legend.pos = "bottomleft"))
dev.off()

# Pairwise: CT vs MRI (schuetz), side-by-side SROC + differences table
data(schuetz)
res <- dta_pairwise(schuetz, studlab = "studlab",
                    intervention.label = "CT", control.label = "MRI")
png("man/figures/tutorial-sroc-pair.png", width = 2600, height = 1350, res = 170)
print(dta_sroc_pair(res, outcome = "coronary artery disease",
                    population = "adults with suspected CAD",
                    legend.pos = "bottomleft"))
dev.off()

cat("man/figures/tutorial-sroc.png, tutorial-sroc-pair.png written\n")
