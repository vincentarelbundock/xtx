#import "/.calepin/calepin.typ" as calepin
#show: calepin.document

#set document(title: "Linear models without the model matrix")
#metadata((
  title: "lm_xtx",
  pagetitle: "Linear models without the model matrix — lm_xtx",
  aliases: ("lm_xtx", "predict.xtx_lm"),
  "source-files": ("R/fe.R", "R/lm.R", "R/methods.R"),
)) <website-metadata>

#calepin.setup(echo: true, eval: true, results: "verbatim", warning: false, message: false)

#calepin.chunk(echo: false, results: "hide")[
```r
library(xtx)
library(xtx.fit)
```
]

= Linear models without the model matrix

== Description

`lm_xtx()` fits a linear model by least squares, like #link("https://rdrr.io/pkg/stats/man/lm.html")[`stats::lm()`], without forming the dense model matrix\. Given a #link("https://rdrr.io/pkg/xtx/man/design_matrix.html")[`xtx::design_matrix()`] and a response `y` instead of a formula and data, it fits that design, like #link("https://rdrr.io/pkg/stats/man/lm.fit.html")[`stats::lm.fit()`]\.

== Usage

```r
#| eval: false
lm_xtx(
  formula,
  data,
  subset,
  weights,
  na.action,
  offset,
  contrasts = NULL,
  compress = "auto",
  tol = 1e-10,
  model = TRUE,
  y = NULL,
  vcov = "iid",
  fixef.rm = "singleton",
  ...
)

# S3 method for class 'xtx_lm'
predict(
  object,
  newdata = NULL,
  se.fit = FALSE,
  interval = c("none", "confidence", "prediction"),
  level = 0.95,
  na.action = stats::na.pass,
  ...
)
```

== Arguments

#table(
  columns: 2,
  table.header([Argument], [Description]),
  [`formula`], [A model formula, as in #link("https://rdrr.io/pkg/stats/man/lm.html")[`stats::lm()`], or a #link("https://rdrr.io/pkg/xtx/man/design_matrix.html")[`xtx::design_matrix()`] \(then give the response as `y`\)\.],
  [`data, subset, weights, na.action, offset, contrasts`], [As in #link("https://rdrr.io/pkg/stats/man/lm.html")[`stats::lm()`] \(`na.action`: `na.omit`, `na.exclude`, `na.fail` or `na.pass`, those remote data can apply, in memory too; with fixed effects `na.omit` or `na.exclude`\)\. With a design matrix, `weights` and `offset` are numeric vectors with one value per row, and `data`, `subset`, `na.action` and `contrasts` do not apply\. `data` can also be remote \(a DuckDB connection with `table` or `sql`, or Parquet paths; see #link("https://rdrr.io/pkg/xtx/man/design_remote.html")[`xtx::design_remote`]\); `weights`, `offset` and `subset` are then columns, `DBI::SQL()` or R expressions over the columns \(see #link("remote_fits.typ")[`remote_fits`]\)\.],
  [`compress`], [`"auto"`, `TRUE` or `FALSE`: fit on the sufficient statistics of rows that share a design pattern \(see Details\)\.],
  [`tol`], [The tolerance for detecting aliased columns \(see #emph[Aliased columns]\)\.],
  [`model`], [Keep the model frame and its design matrix, for #link("https://rdrr.io/r/base/residuals.html")[`residuals()`], #link("https://rdrr.io/r/base/fitted.html")[`fitted()`] and #link("https://rdrr.io/r/base/model.frame.html")[`model.frame()`]\.],
  [`y`], [With a design matrix as first argument, the numeric response, one value per row\.],
  [`vcov`], [With fixed effects: `"iid"` or a one-sided formula naming a cluster variable \(`~firm`\)\. Without fixed effects only `"iid"`\.],
  [`fixef.rm`], [With fixed effects: `"singleton"` \(drop singletons iteratively, as #link("https://rdrr.io/pkg/fixest/man/feols.html")[`fixest::feols()`]; `"perfect_fit"` is the same for a linear model\) or `"none"`\.],
  [`...`], [For remote data \(a DuckDB connection or Parquet paths as `data`, or as `newdata` in `predict()`\): `table`, `sql`, `engine` and `chunk_size`; see #link("https://rdrr.io/pkg/xtx/man/design_remote.html")[`xtx::design_remote`]\. Otherwise not used\.],
  [`object`], [A fit of `lm_xtx()`\.],
  [`newdata`], [A data frame of new observations, or remote data \(a DuckDB connection with `table` or `sql`, or Parquet paths\)\. Omitted, predictions are those of the fit's rows \(which need `model = TRUE`\)\. From a rank-deficient fit, predictions on new data that are not estimable are `NA`, with a warning\.],
  [`se.fit`], [Return standard errors of the predictions?],
  [`interval`], [Type of interval: `"none"`, `"confidence"` or `"prediction"`\.],
  [`level`], [Confidence level of the intervals\.],
)

== Value

An object of class `c("xtx_lm", "xtx_fit")` \(with fixed effects, `"xtx_lm_fe"` first\): a list with, among others, `coefficients`, `sigma`, `df.residual`, `nobs`, `rank`, `cov.unscaled` \(the inverse of `X'WX` on the columns kept, `NA` for aliased columns\), `rss`, `compressed` \(whether the fit was compressed\) and `patterns` \(the number of rows or patterns the fit was computed on\), `spec` \(the design specification\), `call` and `terms`\.

== Details

#strong[Method\.] The fit solves the normal equations from the moments `X'WX` and `X'Wy`, computed from factor codes and numeric variables, by a pivoted Cholesky factorization, on the design with numeric variables far from zero centered \(an exact transform of the design, mapped back\)\. It is the estimator of fits with fixed effects, with none to absorb, in memory and on remote data alike\. When the design is badly conditioned, or the residual sum of squares would lose digits to cancellation, the coefficients are refined and the sums of squares computed from the residuals \(one more pass over the rows\)\. Results agree with those of `lm()` up to floating point rounding\.

#strong[Aliased columns\.] Columns are taken in order; a column is aliased \(its coefficient is `NA`\) when its residual variance after the columns kept before it is below `tol` times its variance, the rule of fits with fixed effects \(and of #link("https://rdrr.io/pkg/fixest/man/feols.html")[`fixest::feols()`]\)\. An exactly collinear column is aliased as in `lm()`\. The rule differs from that of `lm()` for columns that are collinear only numerically: it is applied to the centered design, so an interaction of variables far from zero, which `lm()` drops because their offsets dominate its norm, is kept\.

#strong[Compression\.] When many rows share a design pattern \(the same factor levels and numeric values in every variable of the design\), the fit only needs, for each pattern, the sum of the weights, the weighted mean of the response \(minus the offset\) and the weighted sum of squares around that mean\. The coefficients, residual sum of squares and covariance matrix computed from these sufficient statistics are those of the fit on every row, so a fit on millions of rows of categorical data costs about as much as one on its distinct patterns\. Patterns are found by the deduplication machinery of the products \(see #link("https://rdrr.io/pkg/xtx/man/design_algebra.html")[`xtx::design_algebra`] and the option `xtx.dedup`\)\. With `compress = "auto"`, the fit is compressed when at most half the rows are distinct patterns, and the search for patterns gives up as soon as that is exceeded; the search is skipped when the option `xtx.dedup` is `FALSE`\. `compress = TRUE` always compresses \(unless every row is distinct\) and `FALSE` never does\.

#strong[What is kept\.] The fit keeps the coefficients, their unscaled covariance, the residual and explained sums of squares, and the design specification, which #link("https://rdrr.io/r/base/predict.html")[`predict()`] applies to `newdata` with `design_matrix(spec, newdata)`\. Per-row residuals and fitted values are not stored: with `model = TRUE` \(the default, as in #link("https://rdrr.io/r/base/lm.html")[`lm()`]\) the fit keeps the model frame and its design matrix, which share memory with the frame, and #link("https://rdrr.io/r/base/residuals.html")[`residuals()`], #link("https://rdrr.io/r/base/fitted.html")[`fitted()`] and `predict()` without `newdata` recompute them \(one product with the design\)\. With `model = FALSE` these need `newdata`\.

#strong[Methods\.] #link("https://rdrr.io/r/base/coef.html")[`coef()`], #link("https://rdrr.io/r/base/vcov.html")[`vcov()`] \(with `NA` rows and columns for aliased coefficients, as `vcov(lm_fit)`\), #link("https://rdrr.io/r/base/sigma.html")[`sigma()`], #link("https://rdrr.io/r/base/df.residual.html")[`df.residual()`], #link("https://rdrr.io/r/base/nobs.html")[`nobs()`], #link("https://rdrr.io/r/base/deviance.html")[`deviance()`] \(the weighted residual sum of squares\), #link("https://rdrr.io/r/base/logLik.html")[`logLik()`] \(hence #link("https://rdrr.io/r/base/AIC.html")[`AIC()`] and #link("https://rdrr.io/r/base/BIC.html")[`BIC()`]\), #link("https://rdrr.io/r/base/formula.html")[`formula()`], #link("https://rdrr.io/r/base/terms.html")[`terms()`], #link("https://rdrr.io/r/base/model.frame.html")[`model.frame()`], #link("https://rdrr.io/r/base/residuals.html")[`residuals()`], #link("https://rdrr.io/r/base/fitted.html")[`fitted()`], #link("https://rdrr.io/r/base/predict.html")[`predict()`], `print()` and #link("https://rdrr.io/r/base/summary.html")[`summary()`] \(a coefficient table, R-squared and F statistic, as `summary.lm()`\)\.

#strong[Differences with `lm()`\.] Multivariate responses and `method = "model.frame"` are not supported\. The result is not an `lm` object: it has no `qr`, `effects`, `assign` or `residuals` components, so functions that expect those \(such as `plot.lm()` or `anova.lm()`\) do not apply\. In `predict()`, aliased coefficients count as zero, and predictions on new data that are not estimable from a rank-deficient fit are `NA`, with a warning \(`rankdeficient = "NAwarn"` of `predict.lm()`\)\. Prediction intervals use the variance `sigma^2 / w` at the training rows of a weighted fit and `sigma^2` for new observations, with the warnings of `predict.lm()`; its `weights`, `pred.var` and `scale` arguments are not supported\. With `na.action = na.omit`, offsets of `newdata` are taken for the rows kept \(`predict.lm()` would recycle them\)\.

#strong[Remote data\.] `lm_xtx()` also fits remote data \(a DuckDB table or query, or Parquet files\), with the same arguments and methods; see the #emph[Linear models] section of #link("remote_fits.typ")[`remote_fits`] for how, and for the differences from in-memory fits\.

== Fixed effects

A formula `y ~ x | fe1 + fe2` absorbs fixed effects, as #link("https://rdrr.io/pkg/fixest/man/feols.html")[`fixest::feols()`]: each FE is a variable \(a factor, character, integer or numeric column\) or an interaction `a^b`\. The fit never builds the FE dummies or demeans the rows: it solves the normal equations from sums over rows \(the design's products, sums within FE levels and crosstabs of the FE\), eliminating the largest FE exactly and solving for the others densely \(up to 1000 levels\) or by conjugate gradients: on the reduced system formed, with an incomplete Cholesky preconditioner, when it is sparse \(low mobility, such as workers who rarely change firms\), and without forming it, with a Jacobi preconditioner, when it would be near dense \(high mobility\)\. The same code fits data in memory and remote data \(DuckDB tables or queries, Parquet\), with identical results up to rounding\.

- Rows: rows with a missing value in the model, FE or cluster variables are dropped, as are rows with zero weight and, with `fixef.rm = "singleton"` \(default\), singletons \(rows alone in a level of some FE\), iteratively, as `feols()` does\. `fit$nsingletons` counts them\.
- Coefficients: the intercept is absorbed; a column collinear with the FE or earlier columns is `NA` \(fixest's rule: its residual variance is below `tol` times its variance; `tol` defaults to `1e-10`, as `collin.tol`\)\. The residual variance comes from sums over rows; for a column within rounding of that threshold \(about `2e-12` times its sum of squares\) it is computed again from the data, one pass over the rows, so that the decision does not depend on rounding\. One step of iterative refinement from the residuals is taken, and numeric variables far from zero are centered when that is an exact transform of the design, so that coefficients and standard errors stay accurate\.
- Standard errors: `vcov = "iid"` \(default\) or `vcov = ~cl` \(clustered by one variable\), with fixest's default small sample corrections \(`ssc(K.adj = TRUE, K.fixef = "nonnested", G.adj = TRUE)`\); clustered tests in `summary()` use `t(G - 1)`\. The cluster scores are corrected for FE not nested in the clusters\. `fit$fe_dof` is the number of FE parameters counted \(`sum(levels) - (number of FE - 1)`\); `fit$fe_dof_exact` the exact rank of the FE \(fixest's `K.exact = TRUE`\), exact for one or two FE and for three or more when the reduced system was factored \(otherwise `NA`\)\.
- Methods: those of a `lm_xtx` fit \(class `c("xtx_lm_fe", "xtx_lm")`\); `summary()` reports R-squared and within R-squared\. `fixef()` \(the generic of fixest or nlme, when loaded\) returns the FE values, also in `fit$fixef`: one named vector per FE, normalized with one level per connected set of levels at zero in every FE but the largest \(differences within a connected set are those of fixest\)\. `predict(fit, newdata)` is `X b` plus the FE values of the rows' levels \(`NA` for levels not in the fit\)\. With remote data, `fitted()` and `predict()` return tables and `residuals()` a remote design, as for remote fits without FE\. These methods read the source rows when called; keep the training source available and unchanged for in-sample results\.
- Offsets \(the `offset` argument and `offset()` terms\) are subtracted from the response, as in #link("https://rdrr.io/pkg/fixest/man/feols.html")[`fixest::feols()`]; `fitted()` and `predict()` include them\.
- Not supported with FE: `subset`, `compress = TRUE`, `na.action` other than `na.omit` and `na.exclude` \(in memory and remotely\), `.` in the formula, `se.fit` and intervals in `predict()`, clustering by more than one variable\. For remote data, the response and `weights` are column names, `DBI::SQL()` or R expressions, as without FE \(expressions with no SQL translation are evaluated in R, in batches\)\.
- Remote data: the largest FE stays in the database\. The rows are selected there \(singletons too\), the sums within its levels and its crosstab with the other FE are temporary tables, and R solves the system over the levels of the other FE: memory grows with the coefficients, the levels of the other FE and the links between them, not with the rows nor with the levels of the largest FE\. The values of an FE are read at the end when it has at most `getOption("xtx.fe_levels_max", 5e6)` levels; beyond, they stay in the database: `fit$fixef` holds the name of their table \(columns `key`, the level, and `val`\), and `fitted()`, `residuals()` and `predict()` on remote data of the same connection look them up there\. When the other FE are so densely linked by the levels of the largest that the system over them cannot be formed \(a join of more than `getOption("xtx.fe_join_max", 2e9)` rows\), the fit is an error\. The same holds for #link("glm_xtx.typ")[`glm_xtx()`], whose FE values stay in the database during the iterations\.

== Examples

```r
## fixed effects (examples below: without)
set.seed(2)
panel <- data.frame(id = rep(1:100, each = 10), t = rep(1:10, 100), x = rnorm(1000))
panel$y <- panel$x + panel$id / 50 + rnorm(1000)
fe <- lm_xtx(y ~ x | id + t, panel, vcov = ~id)
summary(fe)
head(fe$fixef$id)

set.seed(1)
d <- data.frame(y = rnorm(1000), x = rnorm(1000),
                f = factor(sample(letters[1:5], 1000, TRUE)),
                g = factor(sample(1:3, 1000, TRUE)))
fit <- lm_xtx(y ~ x * f, d)
all.equal(coef(fit), coef(lm(y ~ x * f, d)))
all.equal(vcov(fit), vcov(lm(y ~ x * f, d)))
summary(fit)

## factors only: fitted on 15 patterns instead of 1000 rows
fit <- lm_xtx(y ~ f * g, d)
fit$patterns
all.equal(sigma(fit), sigma(lm(y ~ f * g, d)))

predict(fit, newdata = d[1:3, ])

## a design matrix and a response, as lm.fit()
X <- xtx::design_matrix(~ x * f, d)
fit <- lm_xtx(X, y = d$y)
all.equal(coef(fit), coef(lm(y ~ x * f, d)))
```
