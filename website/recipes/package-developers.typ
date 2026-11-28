#import "/.calepin/calepin.typ" as calepin
#show: calepin.document

#set document(title: [Building estimators])
#metadata((title: "Building estimators")) <website-metadata>

#calepin.setup(echo: true, eval: true, results: "verbatim", warning: false, message: false)

#calepin.chunk(echo: false, results: "hide")[
```r
options(width = 200)
```
]

#title()

These recipes build estimators from the #link("../operations.html#columns")[columns and groupings]. Each one is written once and runs on a data frame and on a database table. The remote examples need the DBI and duckdb packages.

= Data

The same data, as a data frame and as a DuckDB table: an outcome `y` with a few gross outliers, a numeric and a categorical regressor, and 40 clusters.

```r
library(xtx)
set.seed(48)
n <- 5000
d <- data.frame(
    id = seq_len(n),
    x = rnorm(n),
    f = sample(c("a", "b", "c"), n, TRUE),
    cl = sample(1:40, n, TRUE))
d$y <- 1 + d$x + (d$f == "b") + rnorm(40)[d$cl] + rnorm(n)
out <- sample(n, 50)
d$y[out] <- d$y[out] + 20

con <- DBI::dbConnect(duckdb::duckdb())
DBI::dbWriteTable(con, "d", d)

D_mem <- design_matrix(y ~ x + f, d)
D_db  <- design_matrix(y ~ x + f, con, table = "d")
```

= Least squares with clustered standard errors

The coefficients come from $X'X$ and $X'y$; the residuals are a column; the scores of each cluster are the column sums of the design scaled by the residuals, $sum_(i in g) e_i x_i$, one row per cluster. The sandwich $V = G/(G-1) (X'X)^(-1) S'S (X'X)^(-1)$ is then a small computation in R.

```r
ols_cluster <- function(D, cl) {
    y <- design_column(D, "y")
    XtX <- crossprod(D)
    b <- drop(solve(XtX, crossprod(D, y)))
    e <- y - drop(D %*% b)                  # a column
    S <- colSums(e * D, group = cl)         # G x p cluster scores
    G <- nrow(S)
    bread <- solve(XtX)
    V <- G / (G - 1) * bread %*% crossprod(S) %*% bread
    fit <- cbind(estimate = b, std.error = sqrt(diag(V)))
    attr(fit, "rss") <- sum(e^2)
    fit
}
```

The only difference between the two calls is how the clusters are named: in memory a vector, remotely a column of the table (a design in memory keeps only the variables of its formula).

```r
fit_mem <- ols_cluster(D_mem, d$cl)
fit_db  <- ols_cluster(D_db, "cl")
fit_db
all.equal(fit_mem, fit_db)
```

The same numbers in base R, from the dense model matrix:

```r
X <- model.matrix(y ~ x + f, d)
b <- solve(crossprod(X), crossprod(X, d$y))
e <- drop(d$y - X %*% b)
S <- rowsum(e * X, d$cl)
bread <- solve(crossprod(X))
V <- 40 / 39 * bread %*% crossprod(S) %*% bread
all.equal(unname(fit_db[, "std.error"]), unname(sqrt(diag(V))))
```

= A one-step robust estimator

One step of iteratively reweighted least squares for Huber's M-estimator, from the least squares fit: the residuals are scaled, the rows with large residuals are downweighted with `col_pmin()`, and a weighted least squares step gives the new coefficients. `col_compute()` stores the residuals in a temporary table, since they are used several times (in memory it returns the vector). `col_sums()` reads several sums in one query.

```r
huber_step <- function(D, b, k = 1.345) {
    y <- design_column(D, "y")
    e <- col_compute(y - drop(D %*% b))
    s <- sqrt(sum(e^2) / (nrow(D) - ncol(D)))
    u <- abs(e) / s
    w <- col_pmin(1, k / u)                 # Huber weights
    b1 <- drop(solve(crossprod(D, w * D), crossprod(D, w * y)))
    list(coef = b1,
         scale = s,
         range = range(e),
         sums = col_sums(downweighted = u > k, weight = w,
                         big = col_ifelse(abs(e) > 10, 1, 0)))
}

b0 <- fit_mem[, "estimate"]
step_mem <- huber_step(D_mem, b0)
step_db  <- huber_step(D_db, b0)
step_db
all.equal(step_mem, step_db)
```

The least squares slope is pulled by the outliers; one Huber step brings the coefficients back near those of the data without them:

```r
rbind(ols = b0, huber = step_db$coef,
      clean = coef(lm(y ~ x + f, d[-out, ])))
```

= The within transformation

Fixed-effect estimators work from groupings, described in #link("../operations.html#groupings")[Operations].

The within (one-way fixed effect) slope of `y` on `x` demeans both columns by cluster, $tilde(y) = y - "group_lookup"("group_sums"(y, g) \/ n_g, g)$, and regresses one on the other. The key is given as in the clustered example: in memory by its values, remotely by the name of a column. `g$G` is the number of groups and `g$levels` the key values of each group, in code order.

```r
within_slope <- function(D, g) {
    n_g <- group_sums(1, g)[, 1]              # rows per cluster
    demean <- function(v) v - group_lookup(group_sums(v, g)[, 1] / n_g, g)
    yt <- demean(design_column(D, "y"))
    xt <- demean(design_column(D, "x"))
    col <- col_sums(xy = xt * yt, xx = xt^2)  # one query remotely
    col[["xy"]] / col[["xx"]]
}

g_mem <- design_groups(D_mem, d["cl"])
g_db  <- design_groups(D_db, "cl")
g_db$G
head(g_db$levels)
c(memory = within_slope(D_mem, g_mem), duckdb = within_slope(D_db, g_db),
  lm = coef(lm(y ~ x + factor(cl), d))[["x"]])
all.equal(within_slope(D_mem, g_mem), within_slope(D_db, g_db))
```

= Lazy columns

Printing a remote column shows its dimensions and its SQL, without reading a row. Arithmetic builds a new expression, and the values are computed only when a summary, a product or `as.vector()` needs them.

```r
y <- design_column(D_db, "y")
r <- y - drop(D_db %*% b0)
r
z <- col_ifelse(abs(r) > 10, NA, r)
c(sum = sum(z), sum_na_rm = sum(z, na.rm = TRUE), mean = mean(z, na.rm = TRUE))
head(as.data.frame(round(r, 2), keep = "id", nm = "resid"))
```

The connection stays open while columns are in use. Close it at the end.

```r
DBI::dbDisconnect(con)
```
