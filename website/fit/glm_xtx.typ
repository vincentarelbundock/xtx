#import "/.calepin/calepin.typ" as calepin
#show: calepin.document

#set document(title: [`glm_xtx()`])
#metadata((title: "glm_xtx()")) <website-metadata>

#import "/benchmarks/bench.typ": bench, bench-table, bench-code

#calepin.setup(echo: true, eval: true, results: "verbatim", warning: false, message: false)

#calepin.chunk(echo: false, results: "hide")[
```r
options(width = 200)
```
]

#title()

`glm_xtx()` fits a generalized linear model like `glm()`, without the dense model matrix. The same model with both, on a dataset that comes with R:

```r
library(xtx.fit)

glm(breaks ~ wool + tension, poisson, warpbreaks)
fit <- glm_xtx(breaks ~ wool + tension, poisson, warpbreaks)
fit
```

The iterations are those of `glm.fit()`, step for step (starting values, step halving, convergence rule, warnings), so the coefficients, covariance matrix, deviances, AIC and number of iterations agree with those of `glm()`.

```r
summary(fit)
```

Rows that share a design pattern and an offset share their fitted value, so the fit runs on the patterns: here, 6 patterns for 54 rows.

The families are those of `stats` (binomial, Poisson, Gaussian, Gamma, inverse Gaussian and the quasi families) with their standard links. Every option of `glm_xtx()` also works on remote data (see #link("../out-of-core.html")[out-of-core data]), which is why `subset`, `etastart`, `mustart` and standard errors of predictions are not offered. Fixed effects have their own page: #link("fixed_effects.html")[Fixed effects].

= Benchmark

`glm()` and `glm_xtx()`, a logistic regression on 200,000 rows and 119 columns, on one thread, on top of the #link("../benchmarks.html#setup")[benchmark setup].

#bench("glm")

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

`glm_xtx()` takes remote data as `lm_xtx()` does, with `weights` and `offset` naming columns. With few design patterns, one `GROUP BY` returns the sums of each pattern (including the working weights of the first iteration and the parts of the deviance and AIC that do not depend on the coefficients), and the iterations run in memory. Otherwise, each iteration computes the linear predictor, the working weights and responses, `X'WX` and `X'Wz` in the database (the family and link are translated to SQL). Either way, the results are those of `glm()` on the same rows.

```r
d2 <- transform(d, k = rpois(n, exp(0.2 * z)), yb = as.integer(y > 5))
DBI::dbWriteTable(con, "d2", d2)
fit <- glm_xtx(k ~ x + f, poisson, con, table = "d2")
all.equal(coef(fit), coef(glm(k ~ x + f, poisson, d2)))
fit$method

fit_z <- glm_xtx(yb ~ factor(z) * f, binomial, con, table = "d2")
all.equal(coef(fit_z), coef(glm(yb ~ factor(z) * f, binomial, d2)))
fit_z$method
```

`predict(fit, newdata = con, table = , type = "response")` scores in the database, and `fitted(fit)` and `residuals(fit, type = )` return tables with `.fitted` (and `.resid`) columns next to the columns of the data.

#calepin.chunk(echo: false, results: "hide")[
```r
DBI::dbDisconnect(con)
```
]
