#import "/.calepin/calepin.typ" as calepin
#show: calepin.document

#set document(title: "Generalized linear models without the model matrix")
#metadata((
  title: "glm_xtx",
  pagetitle: "Generalized linear models without the model matrix — glm_xtx",
  aliases: ("glm_xtx", "predict.xtx_glm"),
  "source-files": ("R/glm.R", "R/methods.R"),
)) <website-metadata>

#calepin.setup(echo: true, eval: true, results: "verbatim", warning: false, message: false)

#calepin.chunk(echo: false, results: "hide")[
```r
library(xtx)
library(xtx.fit)
```
]

= Generalized linear models without the model matrix

== Description

`glm_xtx()` fits a generalized linear model by iteratively reweighted least squares \(IRLS\), like #link("https://rdrr.io/pkg/stats/man/glm.html")[`stats::glm()`], without forming the dense model matrix\. It takes a formula and a data frame, a #link("https://rdrr.io/pkg/xtx/man/design_matrix.html")[`xtx::design_matrix()`] and a response `y` \(like #link("https://rdrr.io/pkg/stats/man/glm.fit.html")[`stats::glm.fit()`]\), or remote data \(see #emph[Remote data] below\)\.

== Usage

```r
#| eval: false
glm_xtx(
  formula,
  family = stats::gaussian,
  data,
  weights,
  na.action,
  start = NULL,
  offset,
  control = list(),
  contrasts = NULL,
  compress = "auto",
  model = TRUE,
  y = NULL,
  vcov = "iid",
  fixef.rm = "perfect_fit",
  ...
)

# S3 method for class 'xtx_glm'
predict(
  object,
  newdata = NULL,
  type = c("link", "response"),
  na.action = stats::na.pass,
  ...
)
```

== Arguments

#table(
  columns: 2,
  table.header([Argument], [Description]),
  [`formula`], [A model formula, as in #link("https://rdrr.io/pkg/stats/man/glm.html")[`stats::glm()`], or a #link("https://rdrr.io/pkg/xtx/man/design_matrix.html")[`xtx::design_matrix()`] \(then give the response as `y`\)\.],
  [`family`], [A family object, function or name, as in #link("https://rdrr.io/pkg/stats/man/glm.html")[`stats::glm()`]\.],
  [`data, weights, na.action, start, offset, contrasts`], [As in #link("https://rdrr.io/pkg/stats/man/glm.html")[`stats::glm()`] \(`na.action`: `na.omit`, `na.exclude` or `na.fail`; with fixed effects `na.omit` or `na.exclude`\)\. With a design matrix, `weights` and `offset` are numeric vectors with one value per row, and `data`, `na.action` and `contrasts` do not apply\. For remote data, `weights` and `offset` are column names, `DBI::SQL()` expressions or R expressions over the columns \(evaluated in batches when they have no SQL translation\), as in #link("lm_xtx.typ")[`lm_xtx()`]\.],
  [`control`], [A list of arguments of #link("https://rdrr.io/pkg/stats/man/glm.control.html")[`stats::glm.control()`]\.],
  [`compress`], [`"auto"`, `TRUE` or `FALSE`: fit on the design patterns \(see Details\)\.],
  [`model`], [Keep the model frame and its design matrix, for #link("https://rdrr.io/r/base/residuals.html")[`residuals()`], #link("https://rdrr.io/r/base/fitted.html")[`fitted()`] and #link("https://rdrr.io/r/base/predict.html")[`predict()`] without `newdata`\.],
  [`y`], [With a design matrix as first argument, the response: a numeric vector, a factor, or a two-column matrix of successes and failures \(as in #link("https://rdrr.io/r/base/glm.fit.html")[`glm.fit()`]\)\.],
  [`vcov`], [With fixed effects: `"iid"` or a one-sided formula naming a cluster variable \(`~firm`\)\. Without fixed effects only `"iid"`\.],
  [`fixef.rm`], [With fixed effects: `"perfect_fit"` \(default\), `"singleton"`, `"infinite_coef"` or `"none"`; see #emph[Fixed effects]\.],
  [`...`], [For remote data: `table`, `sql`, `engine` and `chunk_size`; see #link("https://rdrr.io/pkg/xtx/man/design_remote.html")[`xtx::design_remote`]\. Otherwise not used\.],
  [`object`], [A fit of `glm_xtx()`\.],
  [`newdata`], [A data frame of new observations, or remote data \(a DuckDB connection with `table` or `sql`, or Parquet paths\)\. Omitted, predictions are those of the fit's rows \(which need `model = TRUE`\)\. From a rank-deficient fit, predictions on new data that are not estimable are `NA`, with a warning\.],
  [`type`], [`"link"` \(the linear predictor\) or `"response"` \(the fitted means\)\.],
)

== Value

An object of class `c("xtx_glm", "xtx_fit")` \(with fixed effects, `"xtx_glm_fe"` first\): a list with, among others, `coefficients`, `deviance`, `null.deviance`, `aic`, `iter`, `converged`, `boundary`, `df.residual`, `df.null`, `nobs`, `rank`, `cov.unscaled`, `family`, `compressed`, `patterns`, `rows`, `spec`, `call` and `terms`\.

== Details

#strong[Method\.] The iterations are those of #link("https://rdrr.io/r/base/glm.fit.html")[`glm.fit()`]: starting values from `family$initialize` or `start`, a weighted least squares problem at each iteration, solved as #link("lm_xtx.typ")[`lm_xtx()`] solves it \(the normal equations of the design with numeric variables centered, once for the fit; aliased columns by its rule, with tolerance `1e-9`\), step halving when the deviance is not finite or the fitted values are invalid, the convergence rule `|dev - dev_old| / (|dev| + 0.1) < epsilon`, and the same warnings \(non-convergence, boundary, fitted probabilities numerically 0 or 1\)\. The fit reaches the optimum of `glm()`, with the same degrees of freedom and the same aliased coefficients for columns that are exactly collinear \(see #emph[Aliased columns] in #link("lm_xtx.typ")[`lm_xtx()`]\)\. Both stop within `epsilon` of it, not at the same point: at the default `epsilon` coefficients agree to about `1e-4` standard errors, and to rounding as `epsilon` gets smaller\.

#strong[Compression\.] Rows with the same design pattern and the same offset have the same fitted value\. For every family that `glm_xtx()` supports, the deviance has the form `w * (a(y) + y * h(mu) + k(mu))`, so at fixed coefficients a pattern contributes to each iteration only through the sum of its prior weights and its weighted mean response, and to the deviance through a constant that is computed once\. The fit then runs on the patterns, from their starting values, and reaches the optimum of the fit on every row \(it runs on the rows when the starting values of the patterns are not valid and those of the rows are\)\. The AIC, the dispersion and the null deviance are computed from the rows\. `compress` works as in #link("lm_xtx.typ")[`lm_xtx()`]\.

#strong[Families\.] The families of #link("https://rdrr.io/pkg/stats/man/family.html")[`stats::family`]: `binomial`, `poisson`, `gaussian`, `Gamma`, `inverse.gaussian`, `quasibinomial`, `quasipoisson` and `quasi` \(with variance `"constant"`, `"mu(1-mu)"`, `"mu"`, `"mu^2"` or `"mu^3"`\), with the links `identity`, `log`, `logit`, `probit`, `cloglog`, `cauchit`, `inverse`, `sqrt` and `1/mu^2`\. These are the families and links that also have a SQL translation, so that every fit works in memory and on remote data alike; other family objects are an error\.

#strong[Remote data\.] `data` can be a DuckDB connection \(with `table` or `sql`\) or Parquet paths \(see #link("https://rdrr.io/pkg/xtx/man/design_remote.html")[`xtx::design_remote`]\); `weights` and `offset` are then column names, `DBI::SQL()` or R expressions over the columns\. The fit is the same code as in memory: its columns \(response, weights, linear predictor, working weights\) are SQL expressions, and its products with the design run in the database\. See the section #emph[GLMs] of #link("remote_fits.typ")[`remote_fits`]\.

#strong[Fixed effects\.] A formula `y ~ x | fe1 + fe2` absorbs fixed effects \(FE\), as #link("https://rdrr.io/pkg/fixest/man/feglm.html")[`fixest::feglm()`] and #link("https://rdrr.io/pkg/fixest/man/fepois.html")[`fixest::fepois()`]: each FE is a variable or an interaction `a^b`\. The iterations are the same \(starting values, step halving, deviance convergence rule and warnings\), with the slopes and the FE values as parameters; the weighted least squares step of each iteration is the estimator of linear models with FE of #link("lm_xtx.typ")[`lm_xtx()`] \(sums over rows, the largest FE eliminated exactly, the others solved densely or by conjugate gradients\), so the dummies are never formed\. The same code fits data in memory and remote data; remotely the largest FE stays in the database, with the values of the FE during the iterations \(see #emph[Fixed effects] in #link("lm_xtx.typ")[`lm_xtx()`]\)\.

- Rows: rows with a missing value in the model, FE or cluster variables, or a zero weight, are dropped; then, iteratively, rows whose FE level the FE fit perfectly, as fixest's `fixef.rm`: `"perfect_fit"` \(default\) drops singletons and, when the inverse link maps onto \(0, Inf\) or \(0, 1\) \(fixest's rule: log, logit, probit, cloglog links\), levels whose response is all 0 or all 1 \(their FE would be infinite\); `"singleton"`, `"infinite_coef"` and `"none"` drop only the first, only the second, or nothing\. `fit$nsingletons` and `fit$nperfect` count the rows dropped\.
- Coefficients: the intercept is absorbed; a column collinear with the FE or earlier columns is `NA` \(fixest's rule, `collin.tol = 1e-9`\)\.
- Standard errors: `vcov = "iid"` \(default: the inverse information of the last weighted least squares step, times the dispersion\) or `vcov = ~cl` \(clustered by one variable\), with fixest's default small sample corrections, `(n - 1) / (n - K)` with `K` the number of slopes and FE parameters, and `G / (G - 1)` for clusters \(FE nested in the clusters not counted in `K`\)\.
- Methods: those of a `glm_xtx` fit \(class `c("xtx_glm_fe", "xtx_glm")`\); `logLik()`, `AIC()` and `BIC()` count the FE parameters; `fixef()` \(the generic of fixest or nlme, when loaded\) returns the FE values, also in `fit$fixef`, normalized as for #link("lm_xtx.typ")[`lm_xtx()`]; `predict()` gives the link or the response \(`NA` for levels not in the fit\); remotely, per-row results are tables\. The in-sample methods read the training source when called\. `fit$fe_iter` holds the iterations of conjugate gradients of each step \(0 for a direct solve\)\.
- Not supported with FE: `start`, `compress = TRUE` \(compression is turned off\), `na.action` other than `na.omit` and `na.exclude`, `.` in the formula, clustering by more than one variable, a design matrix as first argument\.

#strong[What is kept\.] As in #link("lm_xtx.typ")[`lm_xtx()`]: the coefficients, the unscaled covariance of the last weighted least squares problem \(as in #link("https://rdrr.io/r/base/summary.glm.html")[`summary.glm()`]\), the deviances, AIC, degrees of freedom, the number of iterations, and, with `model = TRUE`, the model frame and design, from which #link("https://rdrr.io/r/base/residuals.html")[`residuals()`], #link("https://rdrr.io/r/base/fitted.html")[`fitted()`] and #link("https://rdrr.io/r/base/predict.html")[`predict()`] are recomputed\.

#strong[Methods\.] #link("https://rdrr.io/r/base/coef.html")[`coef()`], #link("https://rdrr.io/r/base/vcov.html")[`vcov()`] \(`NA` for aliased coefficients, scaled by the dispersion\), #link("https://rdrr.io/r/base/summary.html")[`summary()`] \(coefficient table, dispersion, null and residual deviance, AIC, iterations\), #link("https://rdrr.io/r/base/predict.html")[`predict()`] \(`type = "link"` or `"response"`\), #link("https://rdrr.io/r/base/residuals.html")[`residuals()`] \(types `"deviance"`, `"pearson"`, `"working"` and `"response"`\), #link("https://rdrr.io/r/base/fitted.html")[`fitted()`], #link("https://rdrr.io/r/base/deviance.html")[`deviance()`], #link("https://rdrr.io/r/base/df.residual.html")[`df.residual()`], #link("https://rdrr.io/r/base/nobs.html")[`nobs()`], #link("https://rdrr.io/r/base/logLik.html")[`logLik()`] \(hence #link("https://rdrr.io/r/base/AIC.html")[`AIC()`] and #link("https://rdrr.io/r/base/BIC.html")[`BIC()`]\), #link("https://rdrr.io/r/base/family.html")[`family()`], #link("https://rdrr.io/r/base/formula.html")[`formula()`], #link("https://rdrr.io/r/base/terms.html")[`terms()`], #link("https://rdrr.io/r/base/model.frame.html")[`model.frame()`] and `print()`\.

#strong[Differences with `glm()`\.] The result is not a `glm` object: it has no `qr`, `effects`, `residuals`, `fitted.values` or `weights` components, so functions that expect those \(such as `anova.glm()` or `plot.lm()`\) do not apply\. Every option works in memory and on remote data alike; so `subset`, `etastart`, `mustart`, `method`, `x`, `singular.ok = FALSE`, `na.action` other than `na.omit`, `na.exclude` and `na.fail`, and `se.fit` and `type = "terms"` in `predict()` are not supported\. The null deviance of a fit with an offset and an intercept is refitted, as `glm()` does \(also for a design matrix, where `glm.fit()` would not\)\. Extra arguments do not go to #link("https://rdrr.io/r/base/glm.control.html")[`glm.control()`]: give `control = list(maxit = 50)`\.

== Examples

```r
set.seed(1)
n <- 2000
d <- data.frame(x = rnorm(n), f = factor(sample(letters[1:4], n, TRUE)),
                g = factor(sample(1:3, n, TRUE)))
d$y <- rbinom(n, 1, plogis(d$x + as.integer(d$f) / 4))
fit <- glm_xtx(y ~ x * f, binomial, d)
all.equal(coef(fit), coef(glm(y ~ x * f, binomial, d)))
all.equal(vcov(fit), vcov(glm(y ~ x * f, binomial, d)))
summary(fit)

## factors only: fitted on 12 patterns instead of 2000 rows
d$k <- rpois(n, exp(0.5 + as.integer(d$g) / 3))
fit <- glm_xtx(k ~ f * g, poisson, d)
fit$patterns
all.equal(deviance(fit), deviance(glm(k ~ f * g, poisson, d)))

predict(fit, newdata = d[1:3, ], type = "response")

## fixed effects
panel <- data.frame(id = rep(1:100, each = 10), t = rep(1:10, 100), x = rnorm(1000))
panel$y <- rpois(1000, exp(0.5 * panel$x + panel$id / 100))
fe <- glm_xtx(y ~ x | id + t, poisson, panel, vcov = ~id)
summary(fe)
```
