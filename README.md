# easydta

**Cochrane-compliant Diagnostic Test Accuracy meta-analysis in R.**

A pedagogical, `lme4`-based workflow for bivariate binomial GLMM
meta-analysis of DTA studies, following the
[Cochrane Handbook for Systematic Reviews of Diagnostic Test Accuracy v2.0
(Chapter 10, Supplementary Material 1)](https://training.cochrane.org/handbook-diagnostic-test-accuracy/current).
Replaces `mada` (which uses a normal-approximation Cochrane discourages)
with methods Cochrane explicitly endorses.

Heterogeneity is reported via the **Zhou & Dendukuri (2014) bivariate I²**,
which is defined on the bivariate model and is safe for DTA — unlike the
naive Higgins I² that Cochrane §10.2.5 warns against because of the
threshold effect.

---

## Installation

```r
# lightweight install from a local clone
devtools::load_all("path/to/easydta")
# or
install.packages("path/to/easydta", repos = NULL, type = "source")
```

Dependencies (auto-installed from CRAN if needed via `pacman`):
`lme4`, `msm`, `lmtest`, `ggplot2`, `grid`, `gridExtra`, `MASS`, `pracma`.

**Optional** (recommended for AUC inference):

```r
install.packages("mada")                    # dmetatools runtime dep
remotes::install_github("nomahi/dmetatools")
```

When `dmetatools` is available, `dta_sroc(fit, auc_method = "boot")` is the
default and attaches a bootstrap 95% CI to the AUC. Otherwise it falls
back to trapezoidal integration over the SROC curve.

---

## The three workflows

`easydta` supports three workflows, selected by your input data and
research question.

### Workflow A — Single index test

**When to use:** one diagnostic test, many studies; no arm comparison.

**Input schema:** one row per study with columns
`studlab, TP, FP, FN, TN` (plus any covariates you want to carry through).

**Pipeline:**

```
wide 2x2 (studlab, TP, FP, FN, TN)
      │
      ▼   dta_reshape()        -- wide -> long (sens/spec rows)
long (studlab, sens, spec, true, n)
      │
      ▼   dta_fit_single()     -- bivariate GLMM via lme4::glmer
dta_single object
      │
      ├──► dta_summary()         tidy table: Se, Sp, DOR, LR+, LR- with 95% CI
      ├──► dta_derived()         DOR / LR+ / LR- with delta-method CIs
      ├──► dta_heterogeneity()   tau, rho, Zhou-Dendukuri I^2, prediction region
      ├──► dta_forest()          coupled sens+spec forest as a single plot
      ├──► dta_sroc()            SROC + confidence region + prediction region + AUC
      └──► dta_roc_points()      raw study-level ROC scatter
```

See [examples/example_single.Rmd](examples/example_single.Rmd).

### Workflow B — Pairwise test comparison (paired design, `.e`/`.c` input)

**When to use:** two tests evaluated head-to-head across studies (paired
design), e.g. CT vs MRI, ELISA vs IFA. You want to know whether
sensitivity and/or specificity differ between the two tests.

**Input schema:** one row per study with
`author, year, TP.e, FP.e, FN.e, TN.e, TP.c, FP.c, FN.c, TN.c`,
where `.e` = experimental/intervention arm and `.c` = control arm.

**Pipeline:**

```
wide .e/.c (author, year, TP.e..TN.e, TP.c..TN.c)
      │
      ▼   dta_reshape_pairwise(intervention, control)
long with a `test` column carrying the two arm labels
      │
      ▼   dta_fit_pairwise(test_var = "test")
dta_pairwise object (models A/B/C/D, Cochrane App. 12)
      │
      ├──► dta_summary()    per-test Se, Sp, DOR, LR+, LR- with 95% CI
      ├──► dta_compare()    three LR tests + absolute & relative Se/Sp differences
      └──► (split by arm -> dta_fit_single -> dta_forest / dta_sroc per arm)
```

See [examples/example_pairwise.Rmd](examples/example_pairwise.Rmd).

### Workflow C — Plots only (use the fit you already have)

**When to use:** you've already fit the model and want publication-ready
graphics.

See [examples/example_plots.Rmd](examples/example_plots.Rmd).

---

## Public API

### Data reshaping

| Function | Purpose |
|----------|---------|
| `dta_reshape(data, tp, fp, fn, tn, studlab, extra)` | Single-test wide → long. `studlab` is the column carrying the study label; `extra` lets you carry additional covariates through. |
| `dta_reshape_pairwise(data, author, year, intervention, control, tp.e, fp.e, fn.e, tn.e, tp.c, fp.c, fn.c, tn.c, studlab, test_var)` | Paired-design wide (`.e`/`.c`) → long, stacking the two arms into a single `test` column. `intervention` and `control` are the labels assigned to the `.e` and `.c` rows respectively. If `studlab` is `NULL` (default) the label is built as `paste(author, year)`. |

### Model fitting

| Function | Purpose |
|----------|---------|
| `dta_fit_single(long, nAGQ = 1)` | Bivariate binomial GLMM (Cochrane Appendix 5): `glmer(cbind(true, n - true) ~ 0 + sens + spec + (0 + sens + spec \| studlab))`. Returns an S3 `dta_single` carrying the raw fit, the fixed-effect VCV, the between-study VCV `Psi`, and the long data. |
| `dta_fit_pairwise(long, test_var, nAGQ = 1)` | Fits the four nested models A/B/C/D from Cochrane Appendix 12 using the test-type covariate. Returns a `dta_pairwise` object. |

### Summary & derived measures

| Function | Purpose |
|----------|---------|
| `dta_summary(object)` | Tidy data frame of Se, Sp, DOR, LR+, LR- with 95% CIs. Dispatches on `dta_single` / `dta_pairwise`. |
| `dta_derived(fit)` | DOR / LR+ / LR- with delta-method CIs (`msm::deltamethod`). |

### Heterogeneity

| Function | Purpose |
|----------|---------|
| `dta_heterogeneity(fit)` | τ_sens, τ_spec, correlation ρ, and the **Zhou & Dendukuri (2014) bivariate I²** (per-dimension + joint). Also returns the 95% prediction-region ellipse coordinates used by `dta_sroc()`. |

### Comparison

| Function | Purpose |
|----------|---------|
| `dta_compare(pair_fit)` | Three likelihood-ratio tests (null vs full, full vs Se-common, full vs Sp-common) + absolute and relative Se and Sp differences between the two arms, each with delta-method 95% CIs. |

### Plots (all return a single graphic)

| Function | Purpose |
|----------|---------|
| `dta_forest(fit, data, studlab, ...)` | Coupled sensitivity + specificity forest as **one composite output**: three panels (labels + counts, sens + CI, spec + CI) aligned side-by-side via `gridExtra::grid.arrange(ncol = 3)`. Study-level CIs are exact (Clopper–Pearson). |
| `dta_sroc(fit, ci, pred, auc, auc_method, B, ...)` | SROC curve + 95% confidence region + 95% prediction region + AUC. `auc_method = "boot"` uses `dmetatools::AUC_boot` and attaches a bootstrap CI; `auc_method = "trapz"` integrates the SROC with `pracma::trapz`. |
| `dta_roc_points(data, ...)` | Raw study-level ROC scatter (no model overlay). |

---

## Quick start

```r
library(easydta)

# Workflow A — single test
d    <- read.csv(system.file("extdata/anti_ccp.csv", package = "easydta"))
long <- dta_reshape(d, tp = "TP", fp = "FP", fn = "FN", tn = "TN",
                    studlab = "studlab", extra = "generation")
fit  <- dta_fit_single(long)
print(fit)                       # Se, Sp, DOR, LR+, LR- with CIs
dta_heterogeneity(fit)           # includes Zhou-Dendukuri bivariate I^2
dta_forest(fit, d, studlab = "studlab")
dta_sroc(fit)

# Workflow B — paired-design pairwise
d2    <- read.csv(system.file("extdata/ct_mri.csv", package = "easydta"))
long2 <- dta_reshape_pairwise(d2, intervention = "CT", control = "MRI")
pair  <- dta_fit_pairwise(long2, test_var = "test")
print(pair)                      # per-arm Se, Sp, DOR, LR+/-
dta_compare(pair)                # LR tests + Se/Sp differences
```

---

## Bundled example datasets

| File | Purpose |
|------|---------|
| [`inst/extdata/anti_ccp.csv`](inst/extdata/anti_ccp.csv) | 37-study anti-CCP dataset (Cochrane Handbook Appendix 6). Columns: `studlab, TP, FP, FN, TN, generation`. Used for Workflow A; `generation` (CCP1/CCP2) can also drive a non-paired Workflow B via `extra = "generation"` → `dta_fit_pairwise(test_var = "generation")`. |
| [`inst/extdata/ct_mri.csv`](inst/extdata/ct_mri.csv) | 10 fabricated CT-vs-MRI paired studies. Columns: `author, year, TP.e, FP.e, FN.e, TN.e, TP.c, FP.c, FN.c, TN.c`. Use for Workflow B. |

---

## References

- Takwoingi Y *et al.* *Supplementary material 1 to Chapter 10: Code for
  undertaking meta-analysis.* Cochrane Handbook for Systematic Reviews of
  Diagnostic Test Accuracy v2.0 (July 2023).
- Zhou Y, Dendukuri N. *Statistics for quantifying heterogeneity in
  univariate and bivariate meta-analyses of binary data: the case of
  meta-analyses of diagnostic accuracy.* Stat Med. 2014;33(16):2701-2717.
- Noma H *et al.* [`dmetatools`](https://github.com/nomahi/dmetatools) —
  bootstrap AUC and related diagnostic meta-analysis tools.

## License

MIT.
