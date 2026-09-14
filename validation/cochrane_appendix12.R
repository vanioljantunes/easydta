# Cochrane Handbook for Systematic Reviews of Diagnostic Test Accuracy v2.0
# Supplementary material 1 to Chapter 10 (Takwoingi et al. 2023)
# Appendix 12: Bivariate meta-regression using glmer in R, comparison of CT
# and MRI for coronary artery disease.
#
# The code below is the published code. Lines that differ from the PDF are
# marked "# CHANGED:"; blocks appended at the end are marked "# ADDED:".
# Paste the whole file into R; no files are needed.
# install.packages(c("lme4", "lmtest", "msm"))

#################### 1. DATA IMPORT ############################
# CHANGED: the published code runs setwd("U:/Handbook 2020") and
# (X=read.csv("schuetz.csv")), 108 records. Here: the 5 studies that evaluated
# both CT and MRI (easydta's `schuetz`), one row per study and test.
X = data.frame(
  Study_ID = rep(c("Dewey 2006", "Kefer 2005", "Langer 2009", "Maintz 2007",
                   "Pouleur 2008"), 2),
  Test = rep(c("CT", "MRI"), each = 5),
  TP = c(62, 32, 25, 15, 16,   42, 30, 18, 15, 17),
  FP = c( 5,  6,  2,  2,  7,    2,  9, 15,  1, 17),
  FN = c( 4,  2,  1,  1,  1,    7,  4,  8,  1,  0),
  TN = c(46, 12, 40,  2, 53,   39,  9, 27,  3, 43)
)
################## 2. PREPARE THE DATASETS #####################
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
X$recordid <- 1:10   # CHANGED: 1:108 in the published code (108 records)
### Reshape the data from wide to long format. ###
Y = reshape(X, direction="long", varying=list(c("n1", "n0"),
  c("true1","true0")), timevar="sens", times=c(1,0), v.names=
  c("n","true"))
### Sort data by study to cluster the 2 records per study together. ###
Y = Y[order(Y$id),]
Y$spec<- 1-Y$sens
### Generate a seperate data frame for each test type.
Y.CT = Y[Y$Test=="CT",]
Y.MRI = Y[Y$Test=="MRI",]
########### 3. META-REGRESSION - USING GLMER WITH A COVARIATE ###########
### Meta-analysis of CT ###
(ma_CT = glmer(formula=cbind(true, n - true ) ~ 0 + sens + spec + (0+sens +
  spec|Study_ID), data=Y.CT, family=binomial))
### More detail about the parameter estimates can be obtained by using the
### summary command.
summary(ma_CT)
### Meta-analysis of MRI ###
(ma_MRI = glmer(formula=cbind(true, n - true ) ~ 0 + sens + spec + (0+sens
  + spec|Study_ID), data=Y.MRI, family=binomial))
### More detail about the parameter estimates can be obtained by using the
### summary command.
summary(ma_MRI)

### Fit the model without the covariate. ###
(A = glmer(formula=cbind(true, n - true ) ~ 0 + sens + spec + (0+sens +
  spec|Study_ID), data=Y, family=binomial))
### Add covariate terms to the model for both logit sensitivity and logit
### specificity. This model assumes equal variances for both tests.
Y$CT <- as.numeric((Y$Test == "CT"))
Y$MRI <- as.numeric((Y$Test == "MRI"))
Y$seCT <- (Y$CT)*(Y$sens)
Y$seMRI <- (Y$MRI)*(Y$sens)
Y$spCT <- (Y$CT)*(Y$spec)
Y$spMRI <- (Y$MRI)*(Y$spec)
(B = glmer(formula=cbind(true, n - true) ~ 0 + seCT + seMRI + spCT + spMRI
  + (0+sens + spec|Study_ID), data=Y, family=binomial))
### The models can be formally compared using a LR test.
### Install lmtest package if required (to run remove the #).
# install.packages("lmtest")
### Load the package lmtest.
library(lmtest)
lrtest(A,B)
### Is there a statistically significant difference in sensitivity between
### CT and MRI?
(C = glmer(formula=cbind(true, n - true) ~ 0 + sens + spCT + spMRI +
  (0+sens + spec|Study_ID), data=Y, family=binomial))
lrtest(B,C)
### Is there a statistically significant difference in specificity between
### CT and MRI?
(D = glmer(formula=cbind(true, n - true) ~ 0 + seCT + seMRI + spec +
  (0+sens + spec|Study_ID), data=Y, family=binomial))
lrtest(B,D)
### Different variances for each test
(E = glmer(formula=cbind(true, n - true ) ~ 0 + seCT + seMRI + spCT + spMRI
  +(0 +seMRI + spMRI |Study_ID) +(0 +seCT + spCT |Study_ID), data=Y,
  family=binomial))
### More detail can be obtained by using the summary command.
summary(E)
lrtest(B,E)
lrtest(A,E)
### More detail about the parameter estimates can be obtained by using the
### summary command.
summary(E)
### To obtain the between study covariance between logit sensitivity and
### specificity for each test use
(vcovE = (summary(E))$vcov)
### The overall logit-sensitivity and -specificity are given in
cB = summary(B)$coefficients
cE = summary(E)$coefficients
### Therefore confidence intervals can be extracted with the following
### function
sespci <-function(X) {
  cX = summary(X)$coefficients
  rows = nrow(cX)
  rnames = rownames(summary(X)$coefficients)
  lsesp = matrix(data=NA, nrow=rows, ncol=3)
  for(i in 1:rows) {
    lsesp[i,1] = cX[i,1]
    lsesp[i,2] = cX[i,1] - qnorm(0.975)*cX[i,2]
    lsesp[i,3] = cX[i,1] + qnorm(0.975)*cX[i,2]
  }
  rnames
  sesp = round(100 * plogis(lsesp), digits=1)
  row.names(sesp) <- c(rnames)
  sesp
}
sespci(ma_CT)
sespci(ma_MRI)
sespci(A)
sespci(B)
sespci(C)
sespci(D)
sespci(E)

### Standard errors and confidence intervals for absolute and relative
### differences can be calculated using the delta method.
# This requires the package msm.
# install.packages("msm")
### Load the package msm.
library(msm)
seCT = cE[1,1]
seMRI = cE[2,1]
spCT = cE[3,1]
spMRI = cE[4,1]
se_cov = vcovE[1:2,1:2]
sp_cov = vcovE[3:4,3:4]

### Absolute difference
diff_Se=(exp(seCT)/(1+exp(seCT)))-(exp(seMRI)/(1+exp(seMRI)))
diff_Sp=(exp(spCT)/(1+exp(spCT)))-(exp(spMRI)/(1+exp(spMRI)))
se.diff_Se=deltamethod (~ (exp(x1)/(1+exp(x1)))-(exp(x2)/(1+exp(x2))),
  mean=c(seCT,seMRI), cov=se_cov)
se.diff_Sp=deltamethod (~ (exp(x1)/(1+exp(x1)))-(exp(x2)/(1+exp(x2))),
  mean=c(spCT,spMRI), cov=sp_cov)
data.frame(estimate=c(diff_Se, diff_Sp),
            lci=c(diff_Se-qnorm(0.975)*se.diff_Se, diff_Sp-
              qnorm(0.975)*se.diff_Sp),
            uci=c(diff_Se+qnorm(0.975)*se.diff_Se, diff_Sp+
              qnorm(0.975)*se.diff_Sp),
            row.names=c("Absolute difference Sens", "Absolute difference Spec"))

### Relative difference
rel_Se = (exp(seCT)/(1+exp(seCT)))/(exp(seMRI)/(1+exp(seMRI)))
rel_Sp = (exp(spCT)/(1+exp(spCT)))/(exp(spMRI)/(1+exp(spMRI)))
se.rel_Se = deltamethod (~ log((exp(x1)/(1+exp(x1)))/
  (exp(x2)/(1+exp(x2)))), mean=c(seCT,seMRI), cov=se_cov)
se.rel_Sp = deltamethod (~ log((exp(x1)/(1+exp(x1)))/
  (exp(x2)/(1+exp(x2)))), mean=c(spCT,spMRI), cov=sp_cov)
data.frame(estimate=c(rel_Se, rel_Sp),
            lci=c(exp(log(rel_Se)-qnorm(0.975)*se.rel_Se), exp(
              log(rel_Sp)-qnorm(0.975)*se.rel_Sp)),
            uci=c(exp(log(rel_Se)+qnorm(0.975)*se.rel_Se),
              exp(log(rel_Sp)+qnorm(0.975)*se.rel_Sp)),
            row.names=c("Relative difference Sens", "Relative difference Spec"))

# ADDED: the same absolute differences from model B (equal variances), which
# is what easydta reports by default. Same formulas, B instead of E.
vcovB = (summary(B))$vcov
diff_Se_B = plogis(cB[1,1]) - plogis(cB[2,1])
diff_Sp_B = plogis(cB[3,1]) - plogis(cB[4,1])
se.diff_Se_B = deltamethod (~ (exp(x1)/(1+exp(x1)))-(exp(x2)/(1+exp(x2))),
  mean=c(cB[1,1],cB[2,1]), cov=vcovB[1:2,1:2])
se.diff_Sp_B = deltamethod (~ (exp(x1)/(1+exp(x1)))-(exp(x2)/(1+exp(x2))),
  mean=c(cB[3,1],cB[4,1]), cov=vcovB[3:4,3:4])
data.frame(estimate=c(diff_Se_B, diff_Sp_B),
           lci=c(diff_Se_B-qnorm(0.975)*se.diff_Se_B, diff_Sp_B-qnorm(0.975)*se.diff_Sp_B),
           uci=c(diff_Se_B+qnorm(0.975)*se.diff_Se_B, diff_Sp_B+qnorm(0.975)*se.diff_Sp_B),
           row.names=c("Absolute difference Sens (model B)", "Absolute difference Spec (model B)"))

# ADDED: model E with the two random-effect blocks in the other order (CT
# block first, as easydta writes it). Same model; on these 5 studies the fit
# sits on a boundary and this order reaches a slightly higher log-likelihood.
(E_ctfirst = glmer(formula=cbind(true, n - true ) ~ 0 + seCT + seMRI + spCT + spMRI
  +(0 +seCT + spCT |Study_ID) +(0 +seMRI + spMRI |Study_ID), data=Y,
  family=binomial))
c(logLik_E = as.numeric(logLik(E)), logLik_E_ctfirst = as.numeric(logLik(E_ctfirst)))
sespci(E_ctfirst)
