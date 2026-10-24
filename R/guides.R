# ============================================================================
# guides.R  -  Documentation-only "guide" help pages (no exported objects).
#
# Each block below generates a standalone Rd topic. The numeric prefix in the
# @name keeps them grouped and ordered in the package Help Pages index:
#   1 data preparation -> 2 single-arm -> 3 pairwise -> 4 forests -> 5 funnels
# ============================================================================

#' Guide 1: from the Excel file to a checked data frame
#'
#' Everything that happens before the first model: installing, importing the
#' sheet, the two layouts `easydta` accepts, and the checks worth running
#' before you fit anything.
#'
#' @section 1. Install once:
#' ```r
#' install.packages(c("remotes", "rstudioapi", "readxl"))
#' remotes::install_github("vanioljantunes/easydta")
#' library(easydta); library(rstudioapi); library(readxl)
#' ```
#'
#' @section 2. Import the workbook:
#' \if{html}{\figure{tutorial.png}{options: width="100\%" alt="easydta tutorial"}}
#'
#' Read every sheet into a named list, then pick the sheet you want. Keeping
#' the whole workbook in `ma` means a second outcome is one line away, and the
#' sheet names stay visible:
#' ```r
#' myfile      <- selectFile()                     # or "C:/path/data.xlsx"
#' sheet_names <- excel_sheets(myfile)
#' ma <- lapply(sheet_names, function(x)
#'   as.data.frame(read_excel(myfile, sheet = x)))
#' names(ma) <- sheet_names
#' d <- ma$CCP2
#' ```
#' A CSV works the same way: `d <- read.csv("data.csv")`.
#'
#' @section 3. The two layouts:
#' \strong{Single-arm wide} -- one test, one row per study. Columns
#' `studlab`, `TP`, `FP`, `FN`, `TN`, plus any covariates you want to keep for
#' subgroups. This is what [dta_fit_single()] eats.
#'
#' \strong{Paired wide} -- two tests read in the \emph{same} studies. One row
#' per study, counts suffixed `.e` for the index test and `.c` for the
#' comparator (`TP.e, FP.e, FN.e, TN.e, TP.c, ...`), and one shared `studlab`.
#' This is what [dta_pairwise()] eats.
#'
#' If each study contributed only one of the two tests, you do not have a
#' paired design. Stack the arms into one single-arm frame with a two-level
#' test column and use [dta_compare_tests()] instead. Treating an unpaired
#' comparison as paired invents a within-study correlation that is not in the
#' data.
#'
#' @section 4. Check before you fit:
#' ```r
#' str(d)                       # counts must be numeric, not character
#' colSums(is.na(d[c("TP","FP","FN","TN")]))
#' d$n <- d$TP + d$FP + d$FN + d$TN
#' d[d$TP + d$FN == 0 | d$TN + d$FP == 0, ]   # no diseased, or no non-diseased
#' anyDuplicated(d$studlab)     # one row per study, or the model double counts
#' ```
#' Four things bite most often:
#' \itemize{
#'   \item \strong{Text counts.} A stray space or a comma in Excel turns the
#'     column into `character` and `lme4` fails with an unhelpful message.
#'   \item \strong{Zero cells.} The binomial GLMM handles a zero cell without
#'     a continuity correction, which is the whole point of fitting counts
#'     rather than transformed proportions, so do not add 0.5 yourself. A
#'     study with an empty \emph{margin} (no diseased patients at all) carries
#'     no information about sensitivity and should be dropped.
#'   \item \strong{Several readers or thresholds per study.} Decide the unit
#'     of analysis before fitting: average the readers into one row, or pick
#'     one threshold. Two rows with the same `studlab` are two studies to the
#'     model.
#'   \item \strong{Study labels.} They are printed verbatim in the forest and
#'     the sROC, so write them as you want them to appear ("Smith 2019").
#' }
#'
#' @section 5. Long format, if you ever need it:
#' The fitters work on a long frame with two rows per study, one for
#' sensitivity and one for specificity. [dta_reshape()] and
#' [dta_reshape_pairwise()] build it, and the one-call wrappers call them for
#' you, so you only do this by hand when you want to fit [dta_fit_pairwise()]
#' directly or carry an extra column:
#' ```r
#' long <- dta_reshape(d, tp = "TP", fp = "FP", fn = "FN", tn = "TN",
#'                     studlab = "studlab", extra = "test")
#' ```
#'
#' @examples
#' data(anti_ccp2)   # single-arm wide: studlab, TP, FP, FN, TN, test
#' head(anti_ccp2)
#'
#' data(schuetz)     # paired wide: studlab + .e (CT) / .c (MRI) counts
#' head(schuetz)
#'
#' @seealso [easydta-package-datasets], [dta_reshape()],
#'   [dta_reshape_pairwise()], and the next guide [easydta-2-single-arm].
#' @name easydta-1-data-preparation
NULL

#' Guide 2: one test, from the fit to the finished plots
#'
#' The whole single-test workflow: fit the object, read what it says, then
#' draw it. Every step below names the arguments you actually have to choose.
#'
#' @section 1. Fit:
#' ```r
#' data(anti_ccp2)
#' fit <- dta_fit_single(anti_ccp2, wide = TRUE)
#' ```
#' [dta_fit_single()] fits the Cochrane Appendix-5 model
#' `cbind(true, n - true) ~ 0 + sens + spec + (0 + sens + spec | studlab)`:
#' two fixed effects, the mean logit sensitivity and logit specificity, and a
#' correlated study-level random effect for each. Arguments worth a thought:
#' \describe{
#'   \item{`wide`}{`NA`, the default, auto-detects the layout. Pass
#'     `wide = TRUE` to be explicit, and `tp` / `fp` / `fn` / `tn` /
#'     `studlab` if your columns are named differently.}
#'   \item{`nAGQ`}{`1` is the Laplace approximation, the Handbook default and
#'     what the published Cochrane examples use. Raising it (`nAGQ = 5`)
#'     refines the likelihood but is slower and will not rescue a model that
#'     is singular for lack of data. Leave it at 1 unless you have a reason.}
#'   \item{`conf`}{Confidence level, 0.95 by default. It drives the
#'     prediction-region ellipse stored on the fit and the default of every
#'     summary and plot that reads the fit.}
#' }
#'
#' @section 2. Read the object:
#' ```r
#' print(fit)          # convergence, k, the two pooled logits
#' dta_summary(fit)    # pooled sensitivity and specificity with CIs
#' dta_derived(fit)    # DOR, LR+ and LR- with delta-method CIs
#' ```
#' The fit also carries the heterogeneity block: \eqn{\tau} for each measure,
#' their correlation \eqn{\rho}, and the Zhou-Dendukuri bivariate \eqn{I^2}.
#' Read \eqn{I^2} as the share of the total variation that is between studies
#' rather than sampling noise; in DTA it is routinely above 50 percent and it
#' is a description, not a verdict. A \eqn{\tau} near zero together with a
#' singular-fit warning means the data cannot identify between-study
#' variation, usually because there are too few studies: the pooled point
#' estimate is still usable, the prediction region is not.
#'
#' @section 3. Plot:
#' ```r
#' dta_forest(fit)
#' dta_sroc(fit, test.label = "anti-CCP2",
#'          outcome = "rheumatoid arthritis", population = "adults")
#' dta_funnel(fit, test.label = "anti-CCP2",
#'            outcome = "rheumatoid arthritis", population = "adults")
#' ```
#' [dta_sroc()] is the one with the most choices:
#' \describe{
#'   \item{`test.label`, `outcome`, `population`}{Build the title sentence.
#'     Write them as they should read in print; they are not cosmetic, the
#'     figure is unreadable without them.}
#'   \item{`ci`, `pred`}{The dashed confidence region around the summary point
#'     and the dotted prediction region. The prediction region is the honest
#'     one to show when heterogeneity is high: it says where the next study
#'     would land. It is on by default for a plain plot and off in the
#'     subgroup figure, where the panel is already crowded.}
#'   \item{`summary`}{The in-plot box with Sens, Spec, AUC and \eqn{I^2}. On
#'     by default for a plain plot, off when subgroups are drawn, and you can
#'     force it either way.}
#'   \item{`auc`, `auc_ci`, `B`}{The AUC and its bootstrap CI. `B = 2000`
#'     replicates is the default; drop it to a few hundred while you iterate
#'     on a figure and put it back for the final render.}
#'   \item{`group`, `group.name`, `group.curves`, `group.table`}{A subgroup,
#'     given as a vector named by study label. Each level gets its own symbol,
#'     its own dashed curve refitted on its own studies, and a table under the
#'     panel; the panel title carries the likelihood-ratio test for a subgroup
#'     difference. A level with fewer than two studies is drawn but not
#'     refitted.}
#'   \item{`labels`, `legend.pos`, `legend.style`, `shapes`, `colors`}{Study
#'     labels next to the points, and where the legend sits.}
#' }
#'
#' @section 4. Check how fragile the result is:
#' ```r
#' dta_loo(fit)        # one refit per omitted study, as a table
#' dta_loo_plot(fit)   # the same as an sROC panel with the table beside it
#' ```
#' `auc_ci = TRUE` bootstraps a CI for every refit, which is the slow part,
#' and `B` controls it. Read the table for a study whose removal moves the
#' pooled estimate outside the CI of the full model: that is influence worth
#' naming in the text.
#'
#' @examples
#' data(anti_ccp2)
#' fit <- dta_fit_single(anti_ccp2, wide = TRUE)
#' print(fit)
#' dta_summary(fit)
#' dta_derived(fit)
#'
#' @seealso [dta_fit_single()], [dta_summary()], [dta_derived()], [dta_sroc()],
#'   [dta_loo()], and the next guide [easydta-3-pairwise].
#' @name easydta-2-single-arm
NULL

#' Guide 3: two tests, the design, the variance choice and the tests
#'
#' Comparing two index tests (Cochrane Handbook Appendix 12). The design
#' decides the entry point, and the random-effects structure is the one
#' argument that changes the answer, so it gets a section of its own.
#'
#' @section 1. Which entry point:
#' \tabular{ll}{
#'   Every study ran both tests \tab [dta_pairwise()], paired wide `.e` / `.c` \cr
#'   Each study ran one test \tab [dta_compare_tests()], one stacked frame plus `test_var`
#' }
#' Both return the same `dta_pairwise_result`, so everything downstream is
#' identical. The paired route keeps the within-study correlation, which is
#' what makes a head-to-head comparison more precise than comparing two
#' separate meta-analyses; the between-study route cannot, and its comparison
#' is confounded by whatever else differs between the two sets of studies.
#' ```r
#' data(schuetz)
#' res <- dta_pairwise(schuetz, studlab = "studlab",
#'                     intervention.label = "CT", control.label = "MRI")
#'
#' data(anti_ccp1); data(anti_ccp2)
#' res2 <- dta_compare_tests(rbind(anti_ccp1, anti_ccp2), test_var = "test")
#' ```
#'
#' @section 2. The variance argument:
#' `variance` sets the random-effects structure carried by all four nested
#' models:
#' \describe{
#'   \item{`"unequal"`, the default, Cochrane model E}{Two per-test blocks,
#'     `(0 + seA + spA | studlab) + (0 + seB + spB | studlab)`. Each test gets
#'     its own between-study variances and its own correlation. This is the
#'     default because two tests seldom scatter across studies to the same
#'     degree: a reader-dependent test varies more than an automated one, and
#'     forcing one shared variance on both understates the uncertainty of the
#'     noisier test.}
#'   \item{`"equal"`, Cochrane model B}{One shared block,
#'     `(0 + sens + spec | studlab)`. Three fewer variance parameters, so it
#'     converges where model E does not. Switch to it when the unequal fit is
#'     singular, throws a convergence warning, or returns a correlation pinned
#'     at -1 or 1, which happens routinely with fewer than about ten studies
#'     per test.}
#' }
#' Say in the methods which one you fitted. The two structures can give
#' different difference CIs and therefore different p-values, so the choice
#' belongs in the protocol, not in the output.
#' ```r
#' res_eq <- dta_pairwise(schuetz, studlab = "studlab",
#'                        intervention.label = "CT", control.label = "MRI",
#'                        variance = "equal")
#' ```
#' The other arguments: `nAGQ` as in guide 2, `conf` for every CI in the
#' result, and `intervention.label` / `control.label`, which name the arms
#' everywhere downstream, so set them to the test names rather than leaving
#' the defaults.
#'
#' @section 3. What comes back:
#' ```r
#' print(res)        # the fit, the nested models, the differences
#' res$compare       # LR tests and Se / Sp differences
#' res$arms          # one dta_single fit per test, named by your labels
#' ```
#' The ladder of nested models is fitted in the Handbook's order: A, neither
#' measure differs; B, both may differ; C, only specificity differs; D, only
#' sensitivity differs. Three likelihood-ratio tests follow: the overall
#' 2-degree-of-freedom test of any difference, then one test per measure.
#' Read the overall test first. If it is not significant, reading the two
#' single-measure tests as if they stood alone is multiple testing. The
#' differences themselves come from the delta method and are reported with
#' their CIs; a difference in sensitivity of 0.05 (-0.02 to 0.12) is a result,
#' "p = 0.18" on its own is not.
#'
#' @section 4. Plot:
#' ```r
#' dta_sroc_pair(res, outcome = "coronary artery disease",
#'               population = "adults with suspected CAD")
#' dta_forest_pair(res)
#' dta_forest(res, test = intervention)   # one arm at a time
#' dta_funnel(res, test = control)
#' ```
#' [dta_sroc_pair()] draws the two arms side by side and, with `table = TRUE`,
#' the Cochrane differences table underneath; `group.e` / `group.c` take a
#' subgroup per arm, exactly like `group` in [dta_sroc()].
#' [dta_forest_pair()] plots the per-study difference between the tests, with
#' `which` choosing sensitivity, specificity or both, and `counts` showing the
#' raw cells behind each row.
#'
#' @examples
#' data(schuetz)
#' res <- dta_pairwise(schuetz, studlab = "studlab",
#'                     intervention.label = "CT", control.label = "MRI")
#' print(res)
#'
#' @seealso [dta_pairwise()], [dta_compare_tests()], [dta_fit_pairwise()],
#'   [dta_compare()], [dta_sroc_pair()], and the next guide
#'   [easydta-4-forests].
#' @name easydta-3-pairwise
NULL

#' Guide 4: coupled forest plots and their columns
#'
#' [dta_forest()] draws one composite figure: study labels, the 2x2 counts, a
#' sensitivity forest and a specificity forest on the same rows, and the
#' numeric estimates. It draws itself on the current device, so wrap it in
#' `tiff()` or `png()` and `dev.off()` to save it.
#'
#' @section 1. From a fit, or from one arm of a comparison:
#' ```r
#' data(anti_ccp2)
#' fit <- dta_fit_single(anti_ccp2, wide = TRUE)
#' dta_forest(fit)
#'
#' data(schuetz)
#' res <- dta_pairwise(schuetz, studlab = "studlab",
#'                     intervention.label = "CT", control.label = "MRI")
#' dta_forest(res, test = intervention)   # or test = control, or test = "CT"
#' ```
#'
#' @section 2. The arguments you will reach for:
#' \describe{
#'   \item{`subgroup`, `test.subgroup`}{A study-level grouping, as a vector
#'     named by study label or the name of a column. Each level is drawn as
#'     its own block with its own pooled row; with exactly two levels of at
#'     least two studies, `test.subgroup = TRUE` adds the likelihood-ratio
#'     test for a difference between them.}
#'   \item{`columns`}{Extra per-study columns (year, country, threshold,
#'     risk-of-bias judgements) as a data frame, drawn to the left of the
#'     forests. This is how a QUADAS-2 column gets into the figure.}
#'   \item{`overall`}{`FALSE` drops the pooled diamond. Use it when the
#'     studies are too clinically different to pool and the figure is there to
#'     display them, not to summarise them.}
#'   \item{`title`, `legend`}{Title above and legend below, both left
#'     aligned. One house style puts the descriptive title in the x-axis label
#'     of the measure instead of above the plot, so check what the target
#'     journal wants.}
#'   \item{`digits`, `just`, `conf`}{Display digits, justification of the
#'     numeric columns, and the confidence level.}
#' }
#'
#' @section 3. Saving at publication size:
#' ```r
#' tiff("figure_1_forest.tif", width = 4000, height = 2200,
#'      res = 300, compression = "lzw")
#' dta_forest(fit)
#' dev.off()
#' ```
#' Size the device in R rather than resizing the image afterwards: enlarging a
#' saved bitmap blurs the text, and most journals ask for 300 dpi TIFF.
#'
#' @examples
#' data(anti_ccp2)
#' fit <- dta_fit_single(anti_ccp2, wide = TRUE)
#' dta_forest(fit)
#'
#' \donttest{
#' data(schuetz)
#' res <- dta_pairwise(schuetz, studlab = "studlab",
#'                     intervention.label = "CT", control.label = "MRI")
#' dta_forest(res, test = intervention)
#' }
#'
#' @seealso [dta_forest()], [dta_forest_pair()], and the next guide
#'   [easydta-5-funnels].
#' @name easydta-4-forests
NULL

#' Guide 5: Deeks funnel and reading the asymmetry test
#'
#' The publication-bias diagnostic Cochrane recommends for DTA (Handbook v2.0
#' section 10.6.4; Deeks, Macaskill and Irwig 2005). An ordinary funnel plot
#' of sensitivity is not valid here; Deeks' version is built for the DOR
#' scale.
#'
#' @section 1. Draw it:
#' ```r
#' data(anti_ccp2)
#' fit <- dta_fit_single(anti_ccp2, wide = TRUE)
#' dta_funnel(fit, test.label = "anti-CCP2",
#'            outcome = "rheumatoid arthritis", population = "adults")
#' ```
#' The plot puts `ln(DOR)` on the x axis against `1 / sqrt(ESS)` on the y
#' axis, with that axis reversed so the large studies sit at the top, and
#' draws contour pseudo-confidence bands. The returned object draws itself
#' when printed, and is saved inside a device as in guide 4. For a comparison,
#' pick the arm with `test = intervention` or `test = control`.
#'
#' @section 2. The arguments:
#' \describe{
#'   \item{`continuity`}{Correction added to all four cells of a study with a
#'     zero, only so that its `ln(DOR)` exists. It affects the funnel and the
#'     Deeks test, never the fitted model.}
#'   \item{`table`}{The results table under the plot: slope, its CI and the
#'     p-value. Keep it on, because the figure alone does not carry the test.}
#'   \item{`test.label`, `outcome`, `population`}{The title sentence, as in
#'     [dta_sroc()].}
#'   \item{`conf`}{Confidence level for the pooled DOR CI.}
#' }
#'
#' @section 3. Reading it honestly:
#' A slope p-value below 0.05 is evidence of \emph{asymmetry}, which is
#' consistent with publication bias but also with heterogeneity, with small
#' studies enrolling easier cases, and with chance. Deeks' test has little
#' power with fewer than ten studies: below that, do not run it, and write
#' that it was not assessable rather than reporting a non-significant p-value
#' as reassurance.
#'
#' @examples
#' data(anti_ccp2)
#' fit <- dta_fit_single(anti_ccp2, wide = TRUE)
#' dta_funnel(fit, test.label = "anti-CCP2",
#'            outcome = "rheumatoid arthritis", population = "adults")
#'
#' @seealso [dta_funnel()], [dta_loo()], and the package overview
#'   [easydta-package].
#' @name easydta-5-funnels
NULL
