#import "/.calepin/calepin.typ" as calepin
#show: calepin.document

#set document(title: "Fitting models on remote data")
#metadata((
  title: "remote_fits",
  pagetitle: "Fitting models on remote data — remote_fits",
  aliases: ("remote_fits",),
  "source-files": ("R/remote_fit.R",),
)) <website-metadata>

#calepin.setup(echo: true, eval: true, results: "verbatim", warning: false, message: false)

#calepin.chunk(echo: false, results: "hide")[
```r
library(xtx)
library(xtx.fit)
```
]

= Fitting models on remote data

== Description

#link("lm_xtx.typ")[`lm_xtx()`] and #link("glm_xtx.typ")[`glm_xtx()`] fit remote data \(a DuckDB table or query, or Parquet files; see #link("https://rdrr.io/pkg/xtx/man/design_remote.html")[`xtx::design_remote`]\) with the same arguments and methods as data in memory\. The data stay in the database: the fits run on the remote products of 'xtx'\. A fit takes the remote data and a formula; a remote design made with #link("https://rdrr.io/pkg/xtx/man/design_matrix.html")[`xtx::design_matrix()`] is not fitted directly\.

== Linear models

#link("lm_xtx.typ")[`lm_xtx()`] works on remote data as on a data frame: `lm_xtx(y ~ x + f, con, table = "flights", weights = w, subset = year > 2000)`\. `weights`, `offset` and `subset` are column names \(strings or bare names\), `DBI::SQL()` expressions, or R expressions over the columns, as are the response and `offset()` terms of the formula\. R expressions with a SQL translation \(arithmetic as for the SQL engine; for `subset`, also comparisons, `&`, `|`, `!`, `%in%` with constants and `is.na()`\) run in the database\. Any other R expression, such as `weights = myfun(x)`, is evaluated in R, batch by batch, and its values are written with the columns of the remote data to a temporary table, on which the fit runs \(one more pass over the data, with a message naming the expression; `engine = "sql"` refuses\)\. Compression stays exact, since the computed values are columns of that table\. `na.action` can be #link("https://rdrr.io/r/base/na.omit.html")[`na.omit()`] \(default\), #link("https://rdrr.io/r/base/na.exclude.html")[`na.exclude()`], #link("https://rdrr.io/r/base/na.fail.html")[`na.fail()`] or #link("https://rdrr.io/r/base/na.pass.html")[`na.pass()`], with the results or errors of #link("https://rdrr.io/r/base/lm.html")[`lm()`] \(with fixed effects, `na.omit` or `na.exclude`; other na\.action functions are refused, in memory too\)\. The rows used are those with no missing model variable, response, weight or offset, and factor levels are those of these rows, as in memory\.

The fit is the same algorithm as in memory, on remote primitives: with few design patterns it is compressed \(see #emph[Compression] in #link("https://rdrr.io/pkg/xtx/man/design_remote.html")[`xtx::design_remote`]\); otherwise it solves the normal equations of the remote design \(`X'WX` in the database, suspect columns judged from products with the data\), with numeric variables shifted when that is an exact transform of the design \(so that `x + 1e7` keeps its accuracy\)\. Coefficients, covariance, residual standard error, degrees of freedom, `logLik()`, `AIC()`, `BIC()` and `summary()` are those of #link("https://rdrr.io/r/base/lm.html")[`lm()`] on the same rows, up to rounding\.

Results that have one value per row are tables in the database \(names of temporary views, of class `xtx_remote_table`, with all the columns of the remote data\): `fitted(fit)` \(columns `.fitted` and `.resid`; with `na.exclude`, the rows with missing values too, with `NULL`\), and `predict(fit)` or `predict(fit, newdata = con, table = ...)` \(column `.fitted`; with `se.fit = TRUE`, `.se.fit`, and the attributes `df` and `residual.scale`; with intervals, `.lwr` and `.upr`; from a rank-deficient fit, `NULL` where a prediction on new data is not estimable, with a warning\), computed by the same code as in memory, in SQL\. `residuals(fit)` returns the remote design of the fit with two more columns, `.fitted` and `.resid`, which can scale its rows: for heteroskedasticity-consistent and clustered standard errors,

```
R <- residuals(fit)
meat_hc1 <- crossprod(R, DBI::SQL('".resid" * ".resid"') * R)
scores <- colSums(".resid" * R, group = "cluster")
meat_cl <- crossprod(scores)
```

\(with `bread <- vcov(fit) / sigma(fit)^2`, the HC1 covariance is `n / (n - k) * bread %*% meat_hc1 %*% bread`\)\.

#strong[Differences from in-memory fits\.]

- `model.frame()` is not available for a remote fit, and `model` is ignored: the data stay in the database \(`fit$data$design` is the remote design of the fit\)\.
- Results with one value per row are tables in the database, not vectors: `fitted()` \(`.fitted`, `.resid`\), `predict()` \(`.fitted`; `se.fit` as a column `.se.fit` with attributes `df` and `residual.scale`, instead of a list\), and `residuals()`, a remote design\.
- These tables have no row order: rows are matched to the data through its columns, which every table keeps\.
- `compress = "auto"` decides from other evidence than in memory \(see #emph[Compression] in #link("https://rdrr.io/pkg/xtx/man/design_remote.html")[`xtx::design_remote`]\); the results are the same either way\.

== GLMs

#link("glm_xtx.typ")[`glm_xtx()`] accepts remote data like #link("lm_xtx.typ")[`lm_xtx()`]: `glm_xtx(y ~ x + f, binomial, con, table = "d", weights = w, offset = log(e))`, with `weights` and `offset` column names, `DBI::SQL()` or R expressions over the columns, and `offset()` terms in the formula\. The fit is the same code as in memory, on the remote methods of the products: with few design patterns \(see #emph[Compression] in #link("https://rdrr.io/pkg/xtx/man/design_remote.html")[`xtx::design_remote`]\), the iterations run on the patterns; otherwise each iteration computes its products in the database\. `predict(fit, newdata = con, table = ..., type = "response")`, `fitted()` and `residuals(fit, type = )` return a temporary view with a `.fitted` \(and `.resid`\) column next to the columns of the remote data\. Fixed effects \(`y ~ x | fe`, see #link("glm_xtx.typ")[`glm_xtx()`]\) work remotely as in memory: each iteration's sums over rows are queries, the working weights and responses a temporary table\.

== Examples

```r
if (requireNamespace("duckdb", quietly = TRUE)) {
    con <- DBI::dbConnect(duckdb::duckdb())
    set.seed(1)
    d <- data.frame(x = rnorm(500), f = gl(5, 100))
    d$y <- d$x + as.integer(d$f) + rnorm(500)
    d$k <- rpois(500, exp(0.3 * d$x))
    DBI::dbWriteTable(con, "d", d)
    fit <- lm_xtx(y ~ x + f, con, table = "d")
    all.equal(coef(fit), coef(lm(y ~ x + f, d)))
    summary(fit)
    pois <- glm_xtx(k ~ x + f, poisson, con, table = "d")
    all.equal(coef(pois), coef(glm(k ~ x + f, poisson, d)))
    DBI::dbDisconnect(con)
}
```
