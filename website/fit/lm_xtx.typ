#import "/.calepin/calepin.typ" as calepin
#show: calepin.document

#set document(title: [`lm_xtx()`])
#metadata((title: "lm_xtx()")) <website-metadata>

#import "/benchmarks/bench.typ": bench, bench-table, bench-code

#calepin.setup(echo: true, eval: true, results: "verbatim", warning: false, message: false)

#calepin.chunk(echo: false, results: "hide")[
```r
options(width = 200)
```
]

#title()

`lm_xtx()` fits a linear model like `lm()`, without the dense model matrix. The same model with both, on a dataset that comes with R:

```r
library(xtx.fit)

lm(mpg ~ wt + hp + factor(cyl), mtcars)
fit <- lm_xtx(mpg ~ wt + hp + factor(cyl), mtcars)
fit
```

The coefficients, covariance matrix, residual standard error and degrees of freedom agree with those of `lm()`, aliased columns included. `summary()`, `predict()`, `residuals()`, `fitted()`, `vcov()` and `logLik()` work as for `lm()`.

```r
summary(fit)
```

When many rows share a design pattern (factors and discrete numeric variables), the fit runs on the patterns, with their counts and the sums of the response, and gives the same results as the fit on every row. Fixed effects have their own page: #link("fixed_effects.html")[Fixed effects].

= Benchmark

`lm()` and `lm_xtx()` on 200,000 rows and 119 columns, on one thread, on top of the #link("../benchmarks.html#setup")[benchmark setup].

#bench("lm")

= Least squares

Least squares coefficients only need `X'X` and `X'y`, so `qr()` on a design matrix gives them without forming the dense matrix. With weights, use `qr(sqrt(w) * D)` and `sqrt(w) * y`.

```r
library(xtx)
D <- design_matrix(mpg ~ wt + hp + factor(cyl), mtcars)
all.equal(coef(lm(mpg ~ wt + hp + factor(cyl), mtcars)), qr.coef(qr(D), mtcars$mpg))
```

#calepin.chunk(echo: false, results: "typst")[
```r
## (the number of rows, from the results: the run sets it, OLS_BENCH_N)
ols_res <- read.csv(file.path(normalizePath(file.path(Sys.getenv("XTX_ROOT", "../.."), "bench", "ols")), "results.csv"))
ols_n <- unique(ols_res$n)[1]
cat("The benchmark below compares `qr()` on a design matrix with eight other implementations of least squares, on ",
    format(ols_n, big.mark = ",", scientific = FALSE),
    " rows and designs from ", min(ols_res$p), " to ", max(ols_res$p), " columns. It runs outside this page, because a fair comparison needs a fresh R process for each method and design: `bench/ols/run.sh` runs it and writes `bench/ols/results.csv`, which the tables below read. The protocol:\n", sep = "")
```
]

- Two comparisons. _From formula_: every method starts from the formula and the data, as a user would. _Prebuilt_: the model matrix (or design matrix) is built before timing, and only the fit is timed.
- One thread everywhere, including the BLAS. The model frame is built exactly as `model.matrix()` builds it (option `xtx.skip_na_copy = FALSE`).
- Memory is the peak resident memory of the first call in a fresh process, so it counts every allocation, including those made in C++.
- Every result is compared with the coefficients of `lm()`, including which columns are aliased (`NA`). A method that gives a different answer is marked "no".

#calepin.chunk(echo: false, results: "hide")[
```r
bench_dir <- normalizePath(file.path(Sys.getenv("XTX_ROOT", "../.."), "bench", "ols"))
res <- read.csv(file.path(bench_dir, "results.csv"))
res$p <- as.integer(res$p)
ols_table <- function(mode) {
    x <- res[res$mode == mode, ]
    x$cell <- ifelse(is.na(x$ms), "error",
                     sprintf("%s ms, %s MB%s", format(x$ms, big.mark = ",", trim = TRUE),
                             x$MB, ifelse(x$matches_lm %in% TRUE, "", " (no)")))
    tab <- reshape(x[, c("method", "p", "cell")], idvar = "method", timevar = "p", direction = "wide")
    names(tab) <- c("Method", paste(sort(unique(x$p)), "columns"))
    tab <- tab[order(tab$Method != "xtx qr()", tab$Method), ]
    tinytable::tt(tab)
}
```
]

*From formula*

#calepin.chunk(echo: false, results: "typst")[
```r
ols_table("formula")
```
]

*Prebuilt*

#calepin.chunk(echo: false, results: "typst")[
```r
ols_table("prebuilt")
```
]

#calepin.chunk(echo: false, results: "verbatim")[
```r
cat(readLines(file.path(bench_dir, "machine.txt")), sep = "\n")
```
]

From 20 columns up, `qr()` on a design matrix is the fastest method that matches `lm()` in both comparisons, on every design with a factor. From formula, it is 1.1 to 6 times faster than the next such method, and it uses the least memory, because it never builds the dense matrix. Prebuilt, it is 1.2 to 4 times faster, though some methods need even less memory, since they read the prebuilt matrix in place. The gains are largest when most columns are factor dummies or their interactions, whose zeros the design matrix skips; with many numeric variables and one factor, `qr()` is only slightly ahead (tied prebuilt at 100 columns). On the design with 49 columns, 6 of them aliased, the only other methods that match `lm()` are `speedglm` (3 to 4 times slower) and `.lm.fit()`, the routine that `lm()` uses (6 times slower from formula, 12 times slower prebuilt); the faster alternatives return a different solution or stop with an error.

With numeric variables only, there are no zeros to skip. At 20 columns, `RcppEigen`'s Cholesky decomposition is faster from formula and tied prebuilt; at 50 and 100 columns, `qr()` is 1 to 1.3 times faster than the next method.

With 5 or 6 columns, the dense matrix is small and the gains disappear: `fixest`, which does not build a model frame, is the fastest from formula (by 2 to 5 ms) and usually uses the least memory, and `RcppEigen` is the fastest prebuilt (by 1 to 2 ms).

= On a database

The same call fits data that stay in a DuckDB table or in Parquet files (see #link("../out-of-core.html")[out-of-core data]).

```r
library(xtx)
library(xtx.fit)
con <- DBI::dbConnect(duckdb::duckdb())
set.seed(1)
n <- 10000
d <- data.frame(
    id = seq_len(n),
    x = rnorm(n),
    z = sample(1:5, n, TRUE),
    f = sample(c("a", "b", "c"), n, TRUE),
    w = runif(n, 0.5, 2),
    cl = sample(1:50, n, TRUE))
d$y <- d$x + d$z + rnorm(n) * (1 + abs(d$x))
d$k <- rpois(n, exp(0.2 * d$z))
DBI::dbWriteTable(con, "d", d)
```

`lm_xtx()` fits a linear model on remote data with the same algorithm as in memory: with few design patterns, on the patterns with their counts and sums of the response; otherwise, from `qr()` of the remote design, with `X'X` computed in the database. Either way, the coefficients, aliased columns, covariance and residual standard error are those of `lm()` on the same rows, and `weights` names a column.

```r
fit <- lm_xtx(y ~ x + f, con, table = "d", weights = "w")
all.equal(coef(fit), coef(lm(y ~ x + f, d, weights = w)))
all.equal(vcov(fit), vcov(lm(y ~ x + f, d, weights = w)))

fit_z <- lm_xtx(y ~ factor(z) * f, con, table = "d")
fit_z
```

`weights`, `offset` and `subset` are columns, `DBI::SQL()` expressions or R expressions with a SQL translation, `offset()` terms work, and `na.action` can be `na.omit`, `na.exclude`, `na.fail` or `na.pass`, as for `lm()`. Results with one value per row are tables in the database: `fitted(fit)` (with `.fitted` and `.resid`), and `predict(fit, newdata = con, table = )`, with `se.fit`, confidence and prediction intervals, and the non-estimable cases of a rank-deficient fit, as `predict.lm()`.

```r
fit_s <- lm_xtx(y ~ x + f, con, table = "d", subset = z > 2, offset = log(w))
p <- predict(fit_s, newdata = con, table = "d", interval = "confidence")
head(as.data.frame(p))
```

Expressions with no SQL translation, such as `weights = myfun(x)` or a response `my_transform(y)`, are evaluated in R, batch by batch, into a temporary table on which the fit runs; with `engine = "auto"`, a message names the expression.

#calepin.chunk(echo: false, results: "hide")[
```r
DBI::dbDisconnect(con)
```
]
