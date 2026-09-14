# Cochrane Handbook for Systematic Reviews of Diagnostic Test Accuracy v2.0
# Supplementary material 1 to Chapter 10 (Takwoingi et al. 2023)
# Appendix 5: Meta-analysis of a single index test using glmer in R to fit a
# bivariate model.
#
# The code below is the published code. Lines that differ from the PDF are
# marked "# CHANGED:". Paste the whole file into R; no files are needed.
# install.packages(c("lme4", "msm"))

################## 1. DATA IMPORT ###################################
# CHANGED: the published code runs setwd("U:/Handbook 2020") and
# (X=read.csv("anti-ccp.csv")). The same 37 anti-CCP studies are typed in here.
X = data.frame(
  Study.ID = c("Bas 2003", "Bizzaro 2001", "Goldbach-Mansky 2000", "Jansen 2003",
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

############## 2. META-ANALYSIS OF DTA - USING GLMER ###################
### Install lme4 package if required (to run remove the # and select an
### appropriate CRAN mirror).
# install.packages("lme4")
### Load the package lme4.
library(lme4)
### In order to specify the generalized linear model, first, we need to set
### up the data.
### Generate 5 new variables of type long. We need these before we can
### reshape the data.
# n1 is number diseased
# n0 is number without disease
# true1 is number of true positives
# true0 is the number of true negatives
# recordid is the unique identifier for each observation in the dataset
X$n1 <- X$TP+X$FN
X$n0 <- X$FP+X$TN
X$true1 <- X$TP
X$true0 <- X$TN
X$recordid <- 1:37
### Reshape the data from wide to long format ###
Y = reshape(X, direction="long", varying=list(c("n1", "n0"), c("true1",
  "true0")), timevar="sens", times=c(1,0), v.names=c("n","true"))
### Sort data by study to cluster the 2 records per study together. ###
Y = Y[order(Y$id),]
Y$spec<- 1-Y$sens
#############################
### Perform meta-analysis ###
#############################
## Now, let's examine the model specification and output in more detail.
# The variable true specifies the response while sens and spec are dummy
# variables.
# The fixed effect for logit sensitivity and logit specificity are the
# coefficients of sens and spec.
# The constant term is supressed by adding (0 + ...) to the model formula.
# Adding (0 + sens + spec | study) to the model includes study-level random
# effects.
# family = binomial specifies the data are in binomial form.
# Specified in the form above, the between study covariance matrix is
# unstructured.
# nAGQ controls number of points per axis for evaluating the adaptive
# Gauss-Hermite approximation to the log-likelihood. Defaults to 1,
# corresponding to the Laplace approximation.
# CHANGED: verbose=2 removed (it only prints the optimizer trace).
(MA_Y = glmer(formula=cbind(true, n - true) ~ 0 + sens + spec + (0+sens +
  spec|Study.ID), data=Y, family=binomial, nAGQ=1))
### More detail can be obtained by using the summary command.
(ma_Y = summary(MA_Y))
### To obtain the between study covariance between logit sensitivity and
### specificity for each test use
(summary(MA_Y))$vcov
### For the full list of outputs
labels( ma_Y)
### Therefore, to extract the coefficients
ma_Y$coeff
(lsens = ma_Y$coeff[1,1])
(lspec = ma_Y$coeff[2,1])
se.lsens = ma_Y$coeff[1,2]
se.lspec = ma_Y$coeff[2,2]
### Then we can manually create 95% confidence intervals for logit sens and
### spec
Sens = c(lsens, lsens-qnorm(0.975)*se.lsens, lsens+qnorm(0.975)*se.lsens)
Spec = c(lspec, lspec-qnorm(0.975)*se.lspec, lspec+qnorm(0.975)*se.lspec)
### Or as a dataframe
logit_sesp = data.frame(estimate = c(lsens, lspec),
      lci = c(lsens-qnorm(0.975)*se.lsens, lspec-qnorm(0.975)*se.lspec),
      uci = c(lsens+qnorm(0.975)*se.lsens, lspec+qnorm(0.975)*se.lspec),
      row.names = c("lSens", "lSpec"))
### R has a built in logit and inv.logit function (use qlogis and plogis).
plogis(Sens)
plogis(Spec)
### Obtaining diagnostic odds rato (DOR), positive likelihood ratio (LRp)
### and negative likelihood ratio (LRn)
(DOR = exp(lsens+lspec))
(LRp = plogis(lsens)/(1-plogis(lspec)))
(LRn = ((1-plogis(lsens))/plogis(lspec)))
### Standard errors and confidence intervals of DOR and LRs can be
### calculated using delta method. This requires the package msm.
# install.packages("msm")
library(msm)
se.logDOR = deltamethod (~ (x1+x2), mean=c(lsens,lspec), cov=ma_Y$vcov)
se.logLRp = deltamethod (~ log((exp(x1)/(1+exp(x1)))/(1-
  (exp(x2)/(1+exp(x2))))), mean=c(lsens,lspec), cov=ma_Y$vcov)
se.logLRn = deltamethod (~ log((1-
  (exp(x1)/(1+exp(x1))))/(exp(x2)/(1+exp(x2)))), mean=c(lsens,lspec),
  cov=ma_Y$vcov)
data.frame(estimate = c(DOR, LRp, LRn),
      lci = c(exp(log(DOR)-qnorm(0.975)*se.logDOR), exp(log(LRp)-
        qnorm(0.975)*se.logLRp), exp(log(LRn)-qnorm(0.975)*se.logLRn)),
      uci = c(exp(log(DOR)+qnorm(0.975)*se.logDOR), exp(log(LRp)+
        qnorm(0.975)*se.logLRp), exp(log(LRn)+qnorm(0.975)*se.logLRn)),
      row.names = c("DOR", "LR+", "LR-"))
