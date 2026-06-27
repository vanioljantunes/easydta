# easydta

**Cochrane-compliant Diagnostic Test Accuracy meta-analysis in R.**

---

## 1. What is this package

`easydta` is a pedagogical, opinionated R package for diagnostic test
accuracy (DTA) meta-analysis. It implements the bivariate binomial
GLMM workflow recommended by the
[Cochrane Handbook for Systematic Reviews of Diagnostic Test Accuracy
v2.0 (Chapter 10, Supplementary Material 1, Takwoingi et al. 2023)](https://training.cochrane.org/handbook-diagnostic-test-accuracy/current),
exposes it through a small set of `dta_*()` functions, and replaces
older `mada`-style normal-approximation methods with the ones Cochrane
explicitly endorses.

**What you get:**

- **One model engine** — a bivariate binomial GLMM via
  `lme4::glmer(cbind(true, n - true) ~ 0 + sens + spec + (0 + sens + spec | studlab))`,
  the same parameterisation as Cochrane Appendix 5.
- **Cochrane-style heterogeneity** — the
  Zhou & Dendukuri (2014) bivariate I² (per-dimension and joint),
  plus τ, ρ, and the prediction-region ellipse — computed automatically
  inside `dta_fit_single()`. Avoids the naive Higgins I² that
  Cochrane §10.2.5 warns against because of the threshold effect.
- **Cochrane Appendix 12 inference** — likelihood-ratio tests for
  test-comparison studies (overall test effect, Sens-differs, Spec-differs)
  and delta-method 95% CIs for absolute and relative Se/Sp differences.
- **Plots out of the box** — coupled sens/spec forest, SROC with
  confidence + prediction regions and bootstrap AUC (via
  `dmetatools::AUC_boot` or `AUC_comparison`), and a Deeks funnel with
  contour-enhanced pseudo-CI bands and the Deeks asymmetry test.
- **Sensible defaults; few arguments** — every plot function takes
  `test`/`outcome`/`population` for a Cochrane-style title, and every
  composite plot returns a single self-drawing `gtable` (it renders when
  printed — no `grid::grid.draw()` needed).

---

## 2. How to install it

Install directly from GitHub:

```r
# install.packages("remotes")
remotes::install_github("vanioljantunes/easydta")
```

**Required CRAN dependencies** (auto-installed):
`lme4`, `msm`, `lmtest`, `ggplot2`, `grid`, `gridExtra`, `MASS`,
`metafor`, `pracma`.

**Optional but recommended** — installs paired-bootstrap AUC inference
for `dta_sroc_pair()` (joint AUC + dAUC + p-value via
`dmetatools::AUC_comparison`) and per-test bootstrap CIs for
`dta_sroc(auc_method = "boot")`:

```r
install.packages("mada")                       # dmetatools runtime dep
remotes::install_github("nomahi/dmetatools")
```

When `dmetatools` is unavailable, AUC machinery degrades gracefully:
`dta_sroc()` falls back to trapezoidal integration with a parametric
MVN bootstrap CI on `(lsens, lspec)`, and `dta_sroc_pair(auc_ic = TRUE)`
emits a single warning then drops the dAUC p-value column.

---

## 3. Data preparation

### Input formats

`easydta` accepts two wide layouts depending on the workflow:

**(a) Single-test wide** — one row per study:

| studlab    | TP | FP | FN | TN | (covariates...) |
|------------|----|----|----|----|----------------|
| Smith 2003 | 90 | 10 |  5 | 50 | ...             |

Use this for the single-arm workflow (§4) and also for the pairwise
workflow (§5) when the two arms are encoded as a between-study
covariate (e.g. anti-CCP1 vs anti-CCP2 across different studies).

**(b) Paired wide** — one row per study, `.e` (intervention / index) and
`.c` (control / comparator) suffixes; each study evaluated **both** tests:

| studlab    | TP.e | FP.e | FN.e | TN.e | TP.c | FP.c | FN.c | TN.c |
|------------|------|------|------|------|------|------|------|------|
| Dewey 2006 |  62  |   5  |   4  |  46  |  42  |   2  |   7  |  39  |

Pass this single frame to `dta_pairwise()` (§5). The two arms share a
`studlab`, so the random effect captures the within-study correlation.

### Bundled example datasets

| File / `data()` name | Description |
|----------------------|-------------|
| `data(anti_ccp1)`, `data(anti_ccp2)` (also `inst/extdata/anti_ccp.xlsx`, one sheet per arm: `CCP1`, `CCP2`) | anti-CCP single-arm subsets (Cochrane Handbook ch. 10): 8 CCP1 + 29 CCP2 studies. Columns: `studlab, TP, FP, FN, TN, test`. Each drives the single-arm workflow (§4); `rbind(anti_ccp1, anti_ccp2)` drives the covariate comparison via `dta_compare_tests()`. |
| `data(schuetz)` (also `inst/extdata/schuetz.xlsx`, single sheet `schuetz`) | 5-study CT-vs-MRI dataset for coronary artery disease (Cochrane Appendix 12, direct subset). Wide: `studlab` + `.e` (CT) / `.c` (MRI) counts. Each study did both tests (paired) — drives the wide `dta_pairwise()` workflow (§5). |

### Reshaping (only needed when calling fit functions yourself)

```r
# Single-test wide -> long
long  <- dta_reshape(data,
                     tp = "TP", fp = "FP", fn = "FN", tn = "TN",
                     studlab = "studlab",
                     extra   = "test")   # optional covariate
```

In day-to-day use you can usually skip this step:
`dta_fit_single(wide = TRUE, ...)`, `dta_pairwise(...)`, and
`dta_compare_tests(...)` reshape internally.

---

## 4. Single-arm DTA analysis

**When to use:** one diagnostic test, many studies, no head-to-head
comparison.

```r
library(easydta)
data(anti_ccp2)

fit <- dta_fit_single(anti_ccp2, wide = TRUE,
                      tp = "TP", fp = "FP", fn = "FN", tn = "TN",
                      studlab = "studlab")
print(fit)            # Se, Sp, DOR, LR+, LR- with CIs
                      # plus tau, rho, Zhou-Dendukuri bivariate I^2
```

`fit` is a `dta_single` carrying the raw `glmer` fit, the fixed-effect
VCV, the between-study VCV (`Psi`), the long data, and `$heterogeneity`
(τ_sens, τ_spec, ρ, joint and per-dimension I², the prediction-region
ellipse).

### Plots

```r
# Coupled sens / spec forest -- single composite (5 aligned panels)
dta_forest(fit)

# SROC + 95% confidence region + 95% prediction region + AUC
dta_sroc(fit, test = "anti-CCP2",
         outcome    = "rheumatoid arthritis",
         population = "adults")
```

The vignette: [`examples/example_single.Rmd`](examples/example_single.Rmd).

### Tidy summaries

```r
dta_summary(fit)      # Se, Sp, DOR, LR+, LR- with 95% CIs
dta_derived(fit)      # DOR / LR+ / LR- with delta-method CIs
```

---

## 5. Pairwise DTA analysis

**When to use:** two tests evaluated head-to-head. For a **paired** design
(each study did both tests) pass the wide `.e`/`.c` frame to
`dta_pairwise()`. For a **between-study covariate** design (each study did
one test, e.g. `anti_ccp`) use `dta_compare_tests()`.

### One-call analysis

```r
library(easydta)

# Paired wide frame: one row per study, .e = CT (index), .c = MRI
data(schuetz)
res <- dta_pairwise(schuetz,
                    studlab      = "studlab",
                    intervention = "CT",
                    control      = "MRI")
print(res)
# Prints:
#   - per-arm Se / Sp / DOR / LR+/- summary
#   - LR tests A vs B (overall), C vs B (Sens differs?),
#                D vs B (Spec differs?)
#   - Absolute & relative Se/Sp differences with delta-method 95% CIs

# Variance structure (Cochrane Appendix 12): "equal" = model B (default),
# "unequal" = model E (separate between-study variances per test):
res_E <- dta_pairwise(schuetz, studlab = "studlab",
                      intervention = "CT", control = "MRI",
                      variance = "unequal")
```

For a single frame carrying the test column as a between-study covariate
(e.g. `anti_ccp`), use the sister entry point:

```r
data(anti_ccp1); data(anti_ccp2)
res_cov <- dta_compare_tests(rbind(anti_ccp1, anti_ccp2), test_var = "test")
```

`res` is a `dta_pairwise_result` with `$pair`, `$compare`, and
`$arms` (a named list of per-arm `dta_single` fits), plus `$labels`
identifying which level is intervention and which is control.

### Plots

```r
# Per-arm forest -- pass the test name directly, no $arms[[]] gymnastics
dta_forest(res, arm = "CT")
dta_forest(res, arm = "MRI")

# Side-by-side SROC + Cochrane-style differences table beneath
#   .e = intervention, .c = control
dta_sroc_pair(res,
              arm.e = "CT",  test.e = "CT",
              arm.c = "MRI", test.c = "MRI",
              outcome    = "coronary artery disease",
              population = "adults with suspected CAD",
              auc_ic = TRUE)   # FALSE skips the AUC bootstrap (faster)
```

The differences table reports per-arm Sens, Spec and AUC (each with
95% CI), the absolute difference `.e - .c` (95% CI), and a p-value:

| Measure | CT | MRI | Diff (CT - MRI) | P-value |
|---|---|---|---|---|
| Sensitivity | 0.945 (0.897, 0.971) | 0.861 (0.794, 0.909) | 0.084 (0.017, 0.150) | 0.012 |
| Specificity | 0.861 (0.744, 0.929) | 0.708 (0.546, 0.830) | 0.153 (0.047, 0.259) | 0.001 |

*(AUC row appears when `auc_ic = TRUE`; it requires the
`dmetatools::AUC_comparison()` bootstrap.)*

- Sens / Spec p-values are LR tests (Cochrane Appendix 12, rows 2-3 of
  `res$compare$lr_tests`).
- Sens / Spec diff CIs are delta-method on the joint pairwise model.
- AUC + dAUC + p-value come from a single `dmetatools::AUC_comparison()`
  call (Noma & Matsushima 2020), so the table and per-panel summary
  boxes always show the same numbers.

The vignette: [`examples/example_pairwise.Rmd`](examples/example_pairwise.Rmd).

### Lower-level pairwise API (if you don't want the wrapper)

```r
long  <- dta_reshape_pairwise(schuetz, studlab = "studlab",
                              intervention = "CT", control = "MRI")
pair  <- dta_fit_pairwise(long, test_var = "test")  # models A/B/C/D
cmp   <- dta_compare(pair)                          # LR tests + diffs
```

---

## 6. Small effect analysis

Cochrane-recommended publication-bias diagnostic for DTA reviews
(Handbook v2.0 §10.6.4; Deeks, Macaskill & Irwig 2005).

```r
# Single-arm
dta_funnel(fit,
           test       = "anti-CCP2",
           outcome    = "rheumatoid arthritis",
           population = "adults")

# Per arm of a pairwise comparison
dta_funnel(res, arm = "CT",
           test = "CT",
           outcome    = "coronary artery disease",
           population = "adults with suspected CAD")
```

The plot puts `ln(DOR)` on the x-axis and `1/sqrt(ESS)` on the y-axis
(reversed — large studies on top), with a dashed vertical reference
line at the REML-pooled `ln(DOR)` and three nested pseudo-confidence
triangles (90 / 95 / 99 %, darker toward the centre) so you can read
asymmetry by eye.

Underneath the funnel, horizontally aligned with the reference line, a
small grid table reports:

| Statistic | Value |
|---|---|
| Number of studies | 29 |
| Deeks p-value | 0.246 |

The Deeks p-value comes from a weighted linear regression of
`ln(DOR_i)` on `1/sqrt(ESS_i)` with weights `ESS_i` — slope ≠ 0
indicates funnel asymmetry consistent with publication bias. The full
pooled DOR (REML, 95% CI), regression slope (SE), t and df are
attached for programmatic use:

```r
g <- dta_funnel(fit, ...)
attr(g, "deeks")        # full Deeks regression list
attr(g, "pooled")       # pooled DOR + 95% CI on log and natural scales
attr(g, "study_data")   # per-study TP/FP/FN/TN, lnDOR, SE, ESS
```

`continuity = 0.5` (default) is added to all four cells of any study
with a zero cell — matches `metafor::escalc(to = "only0")`.

---

## References

- Takwoingi Y *et al.* *Supplementary material 1 to Chapter 10: Code
  for undertaking meta-analysis.* Cochrane Handbook for Systematic
  Reviews of Diagnostic Test Accuracy v2.0 (July 2023).
- Zhou Y, Dendukuri N. *Statistics for quantifying heterogeneity in
  univariate and bivariate meta-analyses of binary data: the case of
  meta-analyses of diagnostic accuracy.* Stat Med. 2014;33(16):2701-2717.
- Deeks JJ, Macaskill P, Irwig L. *The performance of tests of
  publication bias and other sample size effects in systematic reviews
  of diagnostic test accuracy was assessed.* J Clin Epidemiol.
  2005;58(9):882-893.
- Noma H, Matsushima Y. *Confidence interval for the AUC of SROC curve
  and some related methods using bootstrap for meta-analysis of
  diagnostic accuracy studies.* arXiv:2004.04339 (2020).
- Noma H *et al.* [`dmetatools`](https://github.com/nomahi/dmetatools) —
  bootstrap AUC and dAUC tools used by `dta_sroc()` / `dta_sroc_pair()`.

## License

MIT.
