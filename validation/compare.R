# Runs the four copy-paste scripts in validation/ and checks the Cochrane
# Handbook results against easydta, number by number.
# Run from the package root: Rscript validation/compare.R
# Set EASYDTA_DEV=1 to test the working copy (devtools::load_all) instead of
# the installed package.

suppressMessages({
  library(lme4)
  library(msm)
  library(lmtest)
})
if (identical(Sys.getenv("EASYDTA_DEV"), "1")) {
  suppressMessages(devtools::load_all(quiet = TRUE))
}

run <- function(file) {
  env <- new.env(parent = globalenv())
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  utils::capture.output(
    suppressMessages(suppressWarnings(sys.source(file, envir = env)))
  )
  env
}

c5  <- run("validation/cochrane_appendix5.R")
e5  <- run("validation/easydta_appendix5.R")
c12 <- run("validation/cochrane_appendix12.R")
e12 <- run("validation/easydta_appendix12.R")

z   <- qnorm(0.975)
tol <- 1e-6
rows <- list()
add <- function(section, measure, cochrane, easydta, note = "CHECK") {
  labs <- if (length(cochrane) == 3) c("", " lower", " upper") else ""
  for (i in seq_along(cochrane)) {
    d <- abs(cochrane[i] - easydta[i])
    rows[[length(rows) + 1]] <<- data.frame(
      section = section, measure = paste0(measure, labs[i]),
      cochrane = cochrane[i], easydta = easydta[i], abs_diff = d,
      status = if (d < tol) "match" else note
    )
  }
}
est3 <- function(tbl, ...) {
  keep <- Reduce(`&`, Map(function(col, val) tbl[[col]] == val,
                          names(list(...)), list(...)))
  as.numeric(tbl[keep, c("estimate", "lci", "uci")])
}
logit3 <- function(cf, name) {
  plogis(cf[name, 1] + c(0, -z, z) * cf[name, 2])
}
ci_exp <- function(est, se) c(est, exp(log(est) - z * se), exp(log(est) + z * se))
ci_lin <- function(est, se) c(est, est - z * se, est + z * se)

# -- Appendix 5: single test -------------------------------------------------
s5 <- dta_summary(e5$fit)
add("App 5", "Sensitivity", with(c5, plogis(Sens)), est3(s5, measure = "Sensitivity"))
add("App 5", "Specificity", with(c5, plogis(Spec)), est3(s5, measure = "Specificity"))
add("App 5", "DOR", with(c5, ci_exp(DOR, se.logDOR)), est3(s5, measure = "DOR"))
add("App 5", "LR+", with(c5, ci_exp(LRp, se.logLRp)), est3(s5, measure = "LR+"))
add("App 5", "LR-", with(c5, ci_exp(LRn, se.logLRn)), est3(s5, measure = "LR-"))

# -- Appendix 12, equal variances (models A-D) -------------------------------
lr_e <- e12$res$compare$lr_tests
lr_c <- list(lrtest(c12$A, c12$B), lrtest(c12$B, c12$C), lrtest(c12$B, c12$D))
lr_lab <- c("LR A vs B", "LR B vs C", "LR B vs D")
for (i in 1:3) {
  add("App 12 (B)", paste(lr_lab[i], "chi-sq"), abs(lr_c[[i]]$Chisq[2]), lr_e$chisq[i])
  add("App 12 (B)", paste(lr_lab[i], "p"), lr_c[[i]][["Pr(>Chisq)"]][2], lr_e$p_value[i])
}
sB <- dta_summary(e12$res$pair)
cB <- summary(c12$B)$coefficients
add("App 12 (B)", "CT sensitivity",  logit3(cB, "seCT"),  est3(sB, test = "CT",  measure = "Sensitivity"))
add("App 12 (B)", "CT specificity",  logit3(cB, "spCT"),  est3(sB, test = "CT",  measure = "Specificity"))
add("App 12 (B)", "MRI sensitivity", logit3(cB, "seMRI"), est3(sB, test = "MRI", measure = "Sensitivity"))
add("App 12 (B)", "MRI specificity", logit3(cB, "spMRI"), est3(sB, test = "MRI", measure = "Specificity"))
dB <- e12$res$compare$differences
add("App 12 (B)", "Sens difference", with(c12, ci_lin(diff_Se_B, se.diff_Se_B)),
    as.numeric(dB[1, c("estimate", "lci", "uci")]))
add("App 12 (B)", "Spec difference", with(c12, ci_lin(diff_Sp_B, se.diff_Sp_B)),
    as.numeric(dB[2, c("estimate", "lci", "uci")]))
cCT  <- summary(c12$ma_CT)$coefficients
cMRI <- summary(c12$ma_MRI)$coefficients
aCT  <- dta_summary(e12$res$arms$CT)
aMRI <- dta_summary(e12$res$arms$MRI)
add("App 12 (arms)", "CT alone sensitivity",  logit3(cCT, "sens"),  est3(aCT,  measure = "Sensitivity"))
add("App 12 (arms)", "CT alone specificity",  logit3(cCT, "spec"),  est3(aCT,  measure = "Specificity"))
add("App 12 (arms)", "MRI alone sensitivity", logit3(cMRI, "sens"), est3(aMRI, measure = "Sensitivity"))
add("App 12 (arms)", "MRI alone specificity", logit3(cMRI, "spec"), est3(aMRI, measure = "Specificity"))

# -- Appendix 12, unequal variances (model E) --------------------------------
note_E <- "block order"
sE <- dta_summary(e12$resE$pair)
cE <- summary(c12$E)$coefficients
add("App 12 (E)", "CT sensitivity",  logit3(cE, "seCT"),  est3(sE, test = "CT",  measure = "Sensitivity"), note_E)
add("App 12 (E)", "CT specificity",  logit3(cE, "spCT"),  est3(sE, test = "CT",  measure = "Specificity"), note_E)
add("App 12 (E)", "MRI sensitivity", logit3(cE, "seMRI"), est3(sE, test = "MRI", measure = "Sensitivity"), note_E)
add("App 12 (E)", "MRI specificity", logit3(cE, "spMRI"), est3(sE, test = "MRI", measure = "Specificity"), note_E)
dE <- e12$resE$compare$differences
add("App 12 (E)", "Sens difference", with(c12, ci_lin(diff_Se, se.diff_Se)),
    as.numeric(dE[1, c("estimate", "lci", "uci")]), note_E)
add("App 12 (E)", "Spec difference", with(c12, ci_lin(diff_Sp, se.diff_Sp)),
    as.numeric(dE[2, c("estimate", "lci", "uci")]), note_E)
lrBE_c <- lrtest(c12$B, c12$E)
lrBE_e <- lrtest(e12$res$pair$models$B, e12$resE$pair$models$B)
add("App 12 (E)", "LR B vs E chi-sq", lrBE_c$Chisq[2], lrBE_e$Chisq[2], note_E)
add("App 12 (E)", "LR B vs E p", lrBE_c[["Pr(>Chisq)"]][2], lrBE_e[["Pr(>Chisq)"]][2], note_E)
add("App 12 (E)", "logLik", as.numeric(logLik(c12$E)), as.numeric(logLik(e12$resE$pair$models$B)), note_E)

# Model E with the CT block first (easydta's order) must match exactly
cE2 <- summary(c12$E_ctfirst)$coefficients
add("App 12 (E, CT block first)", "CT sensitivity",  logit3(cE2, "seCT"),  est3(sE, test = "CT",  measure = "Sensitivity"))
add("App 12 (E, CT block first)", "MRI sensitivity", logit3(cE2, "seMRI"), est3(sE, test = "MRI", measure = "Sensitivity"))
add("App 12 (E, CT block first)", "CT specificity",  logit3(cE2, "spCT"),  est3(sE, test = "CT",  measure = "Specificity"))
add("App 12 (E, CT block first)", "MRI specificity", logit3(cE2, "spMRI"), est3(sE, test = "MRI", measure = "Specificity"))
add("App 12 (E, CT block first)", "logLik", as.numeric(logLik(c12$E_ctfirst)),
    as.numeric(logLik(e12$resE$pair$models$B)))

tbl <- do.call(rbind, rows)
options(width = 200)
print(format(tbl, digits = 6), row.names = FALSE)
cat("\n", sum(tbl$status == "match"), " match (|diff| < ", tol, "), ",
    sum(tbl$status == note_E), " model E block-order rows, ",
    sum(tbl$status == "CHECK"), " to check\n", sep = "")
if (any(tbl$status == "CHECK")) quit(status = 1)
