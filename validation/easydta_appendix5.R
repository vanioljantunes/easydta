# easydta: the Cochrane Appendix 5 analysis (single test, bivariate GLMM)
# on the same 37 anti-CCP studies. Paste the whole file into R.
# install.packages("remotes"); remotes::install_github("vanioljantunes/easydta")

library(easydta)

# Same data as cochrane_appendix5.R (also: rbind(anti_ccp1, anti_ccp2))
d <- data.frame(
  studlab = c("Bas 2003", "Bizzaro 2001", "Goldbach-Mansky 2000", "Jansen 2003",
              "Saraux 2003", "Schellekens 2000", "Vincent 2002", "Zeng 2003",
              "Aotsuka 2005", "Bombardieri 2004", "Choi 2005", "Correa 2004",
              "De Rycke 2004", "Dubucquoi 2004", "Fernandez-Suarez 2005",
              "Garcia-Berrocal 2005", "Girelli 2004", "Greiner 2005",
              "Grootenboer-Mignot 2004", "Hitchon 2004", "Kamali 2005",
              "Kumagai 2004", "Kwok 2005", "Lee 2003", "Lopez-Hoyos 2004",
              "Nell 2005", "Nielen 2005", "Quinn 2006", "Rantapaa-Dahlqvist 2003",
              "Raza 2005", "Sauerland 2005", "Soderlin 2004", "Suzuki 2003",
              "Vallbracht 2004", "van Gaalen 2005", "van Venrooij 2004",
              "Vittecoq 2004"),
  TP = c(110, 40, 43, 110, 40, 72, 139, 90, 115, 23, 236, 74, 89, 90, 31, 69, 25,
         70, 167, 26, 26, 64, 71, 68, 38, 42, 149, 147, 47, 24, 171, 7, 481, 190,
         82, 865, 69),
  FP = c(24, 5, 1, 3, 11, 14, 7, 7, 17, 0, 20, 11, 4, 2, 0, 8, 2, 5, 8, 8, 1, 14,
         2, 14, 3, 2, 7, 10, 7, 3, 26, 2, 23, 12, 13, 79, 5),
  FN = c(86, 58, 63, 148, 46, 77, 101, 101, 16, 7, 88, 8, 29, 50, 22, 18, 10, 17,
         98, 15, 20, 15, 58, 35, 0, 60, 109, 35, 20, 18, 60, 9, 68, 105, 71, 252,
         107),
  TN = c(215, 227, 120, 118, 146, 298, 464, 313, 73, 39, 231, 130, 142, 129, 75,
         38, 40, 228, 88, 15, 56, 293, 66, 132, 73, 96, 114, 106, 375, 79, 443,
         51, 185, 408, 301, 2218, 133)
)

# Same model: glmer(cbind(true, n - true) ~ 0 + sens + spec +
#                   (0 + sens + spec | studlab), family = binomial, nAGQ = 1)
fit <- dta_fit_single(d, tp = "TP", fp = "FP", fn = "FN", tn = "TN",
                      studlab = "studlab")
print(fit)

# Cochrane: plogis(Sens), plogis(Spec) and the DOR / LR+ / LR- data frame
dta_summary(fit)

# easydta adds: between-study SDs, correlation, Zhou-Dendukuri I^2
fit$heterogeneity[c("tau_sens", "tau_spec", "rho",
                    "I2_sens", "I2_spec", "I2_biv")]

# easydta adds: SROC with confidence and prediction regions, AUC
dta_sroc(fit, test.label = "anti-CCP", outcome = "rheumatoid arthritis",
         population = "adults", legend.pos = "bottomleft")
