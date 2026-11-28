#import "/.calepin/calepin.typ" as calepin
#show: calepin.document

#set document(title: [Algebra])
#metadata((title: "Algebra")) <website-metadata>

#import "/benchmarks/bench.typ": bench, bench-table, bench-code

#calepin.setup(echo: true, eval: true, results: "verbatim", warning: false, message: false)

#calepin.chunk(echo: false, results: "hide")[
```r
options(width = 200)
```
]

#title()

These examples go with the tables on the #link("../operations.html")[Operations] page. Each one checks a small result against base R, with `identical()` when the result is exact and `all.equal()` when sums run in a different order. The benchmarks run on the #link("../benchmarks.html#setup")[benchmark setup], in a fresh R process each.

```r
library(xtx)
library(Matrix)
d <- data.frame(
    y = c(0.3, -1.2, 0.8, 1.5, -0.4, 0.9, -2.1, 0.2),
    x = c(-1.5, 0.2, 1.1, -0.3, 2.0, 0.7, -0.9, 1.4),
    f = factor(c("a", "b", "c", "a", "b", "c", "a", "b")),
    cl = c(1, 1, 2, 2, 3, 3, 4, 4))
M <- model.matrix(y ~ x * f, d)
D <- design_matrix(y ~ x * f, d)
```

= Building a design

A design matrix holds the prepared variables and the layout of the columns, but not the dense matrix. `as.matrix()` builds the latter.

```r
D
a <- model.matrix(y ~ x * f, d)
b <- as.matrix(D)
identical(a, b, num.eq = FALSE)
```

= Threads

Writing the dense output accounts for much of the work in both implementations. Additional threads can speed things up (set with `options(xtx.threads = )`: an integer, `0` for all available threads, or a fraction of them; the default is half). The memory column includes both the dense output and temporary allocations made while building it.

#bench("threads")

= Reusing a design on new data

To make predictions, new rows must be turned into columns exactly as the training rows were. Here, `poly(x, 2)` must use the training polynomial, and `f` must keep one column per training level even though the new data only contain `"b"`. Pass the training design and the new data to `design_matrix()`, then use `as.matrix()` to display the columns.

```r
Dp <- design_matrix(y ~ poly(x, 2) + f, d)
nd <- data.frame(x = c(0, 1), f = "b")
as.matrix(design_matrix(Dp, nd))
```

A fitted model works too, and gives the same result:

```r
fit <- lm(y ~ poly(x, 2) + f, d)
identical(as.matrix(design_matrix(Dp, nd)), as.matrix(design_matrix(fit, newdata = nd)), num.eq = FALSE)
```

When you have many batches, prepare once with `spec <- design_matrix(fit, newdata = NULL)`, then call `design_matrix(spec, batch)` for each batch. The benchmark below does this for a 100-row batch.

#bench("reuse")

= Arithmetic

Multiply or divide a design matrix by a number to scale it, or by a vector to scale each row. You can also add or subtract two designs with matching rows and columns. The result is another design matrix.

```r
Dx <- design_matrix(~ x, data.frame(x = c(2, 3)))
as.matrix(c(1, 2) * Dx)
```

= Matrix products

`D %*% B` computes `X B` from the factor codes and numeric variables. `B` can have several columns.

```r
beta <- c(1, 0.5, -1, 2, 0.1, 0.3)

a <- M %*% beta
b <- D %*% beta
all.equal(a, b, check.attributes = FALSE)
```

`crossprod(D, W)` computes `X'W`, such as `X'y` in least squares or weighted sums of the columns of `X`. `crossprod(D)` computes `X'X`, and `crossprod(D, w * D)` the weighted `X'WX`; both return ordinary matrices, which go straight to `chol()`, `chol2inv()`, `eigen()` or `solve()`.

```r
w <- seq(-1, 1, length.out = 8)

a <- crossprod(M, w)
b <- crossprod(D, w)
all.equal(a, b, check.attributes = FALSE)

all.equal(crossprod(M), crossprod(D))
all.equal(crossprod(M, w * M), crossprod(D, w * D))
```

For `X'X`, the design matrix adds up the products of each row's few non-zero entries, instead of all `p^2 / 2` products of the dense row.

= QR decomposition

`qr(D)` decomposes a design matrix without forming `X` or `Q`, and the usual functions accept the result: `qr.coef()`, `qr.fitted()`, `qr.resid()`, `qr.qty()`, `qr.qy()`, `qr.Q()`, `qr.R()` and `solve()`. `R` is the Cholesky factor of `crossprod(D)`, with columns pivoted as in `qr()` and `lm()`, and `Q` is `X R^-1`, never formed. Coefficients are refined from the residuals `y - X b`, which recovers much of the accuracy that the normal equations lose on badly conditioned designs. Columns whose residual norm is close to the tolerance are checked against the data, so aliased columns are those of `qr()`.

Results agree with `qr(as.matrix(D))` up to floating point rounding, not bit for bit. Signs of the rows of `R` (and of the columns of `Q`) can differ, and `qr.Q()` and `qr.qty()` return only the first `rank` columns or rows, since the complete `Q` is never formed. See #link("../reference/design_qr.html")[`design_qr`]. The #link("../fit/lm_xtx.html#least-squares-with-qr")[`lm_xtx()` page] benchmarks it against other least squares implementations.

```r
a <- qr(M)
b <- qr(D)
b
all.equal(qr.coef(a, d$y), qr.coef(b, d$y))
all.equal(qr.resid(a, d$y), qr.resid(b, d$y), check.attributes = FALSE)
all.equal(abs(qr.R(a)), abs(qr.R(b)), check.attributes = FALSE)
```

= Sums and means

The row and column summaries of base R work on design matrices. Pass one group label per row to `colSums()` or `colMeans()` to get one row per group, as with `base::rowsum()`. Group sums are a building block of cluster-robust standard errors and subgroup averages.

```r
sites <- data.frame(site = c("B", "A", "B", "B", "A"), x = c(2, 10, 4, 6, 20))
Ds <- design_matrix(~ x, sites)

colMeans(Ds)
colSums(Ds, group = sites$site)
colMeans(Ds, group = sites$site)

all.equal(colSums(Ds, group = sites$site), rowsum(model.matrix(~ x, sites), sites$site))
```

To weight the rows, multiply the design by a vector first, as in `colSums(w * D, group = g)`.

= Quadratic diagonal

A common task is to extract the diagonal of `X V X'`. For example, in a linear regression context, with `V` the variance-covariance matrix of the coefficients and `X` the design matrix, the diagonal holds the variance of each fitted value, and its square root is the standard error. With `V = solve(crossprod(X))` instead, the diagonal holds the leverage of each observation. `quadratic_diagonal()` computes this diagonal without forming `X` or the full `n` by `n` product.

```r
V <- diag(6) / 10

a <- rowSums((M %*% V) * M)
b <- quadratic_diagonal(D, V)
all.equal(a, b, check.attributes = FALSE)
```

= Sparse matrices

`as(D, "CsparseMatrix")` builds a sparse matrix without forming the dense one. It is identical to converting the dense model matrix.

```r
dense_to_sparse <- function(M) as(as(as(M, "dMatrix"), "generalMatrix"), "CsparseMatrix")
a <- dense_to_sparse(M)
b <- as(D, "CsparseMatrix")
identical(a, b)
```

`Matrix::sparse.model.matrix()` and `MatrixModels::model.Matrix()` are shown for reference: they have their own rules for contrasts and column names. This benchmark uses a factor interaction and checks that all four methods produce the same dimensions and entries.

#bench("sparse")

= Repeated rows

When rows have identical predictor values, `xtx` can compute the result once per distinct row and reuse it for every matching row. Row-wise results stay identical bit for bit. The option `xtx.dedup` controls this: `"auto"` (default), `TRUE` (always) or `FALSE` (never).

Repeated rows are found once per design matrix and remembered, so the first call pays for the search. The first two timings include design construction. The cached timing reuses both the design and its groups. Here, 5,000 patterns repeat over 200,000 rows.

#bench("dedup")

= DuckDB and Parquet

Every operation above also runs when the data stay in a DuckDB table or in Parquet files (see #link("../out-of-core.html")[out-of-core data]): the design takes a connection and a table, a `sql` query, or the path of the files, instead of a data frame. These examples need the DBI and duckdb packages.

```r
library(xtx)
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

`design_matrix(formula, con, table = )` returns a remote design. Factor levels are pinned in R, with R's sorting rules rather than the database collation, and rows with a missing value in a model variable are dropped, as `na.omit()` drops them.

```r
D <- design_matrix(y ~ x * f + log(z), con, table = "d")
D
```

The same works with `sql = "SELECT ..."`, or with Parquet files, which a private DuckDB connection reads:

```r
pq <- tempfile(fileext = ".parquet")
DBI::dbExecute(con, sprintf("COPY d TO '%s' (FORMAT PARQUET)", pq))
design_matrix(y ~ x * f + log(z), pq)
```

The SQL engine handles factors, character and logical columns, numeric columns, interactions, and arithmetic with `log()`, `exp()`, `sqrt()` and `abs()`. Any other term, such as `poly()` or `splines::ns()`, is evaluated in R once, batch by batch, into columns of a temporary table, and the design then reads those columns in SQL like any others. `summary()` shows which terms SQL computes, and `summary(D, sql = TRUE)` the queries.

```r
D2 <- design_matrix(y ~ splines::ns(x, 3) + f, con, table = "d")
summary(D2)
```

To fit in R and score in the database, pass a design or a fitted model instead of the formula: its levels, contrasts and transformations are used, and rows with missing values are kept, with `NULL` predictions, as in memory.

= Algebra in the database

`crossprod()`, `colSums()`, `colMeans()`, `rowSums()`, `quadratic_diagonal()`, `qr()` and its methods, `as.matrix()` and sparse coercion (for small data), row selection with `D[DBI::SQL("x > 0"), ]` and row scaling work as in memory. Rows are scaled by a column name or a SQL expression (`DBI::SQL()`) instead of a vector, so `crossprod(D, "w" * D)` is `X'WX`:

```r
X <- model.matrix(y ~ x * f + log(z), d)
all.equal(crossprod(D), crossprod(X), check.attributes = FALSE)
all.equal(crossprod(D, "w" * D), crossprod(X, d$w * X), check.attributes = FALSE)
all.equal(crossprod("w" * D, "y"), crossprod(X, d$w * d$y), check.attributes = FALSE)
all.equal(colSums(D, group = "cl"), rowsum(X, d$cl), check.attributes = FALSE)
```

When the rows have few distinct design patterns (factors and discrete numeric variables), these products and `lm_xtx()` run one `GROUP BY` over the input columns of the design, and compute the result from the design of the patterns, in memory. The option `xtx.remote_compress` controls this.

= Columns in the database

`D %*% beta`, `rowSums(D)` and `quadratic_diagonal(D, V)` return a remote column, a `design_column`: a SQL expression over the rows of the design, which nothing evaluates until its values are asked for. `design_column(D, "y")` makes one from a column of the data (or from a `DBI::SQL()` expression), and arithmetic on columns (`+`, `*`, `exp()`, `log()`, comparisons, `is.na()`) builds new ones. `as.vector()` and `as.matrix()` read the values into R, `as.data.frame()` reads them as a data frame, with columns of the data given by `keep`, and `col_compute()` stores a column in a temporary table, so that later queries read it instead of recomputing it; with `name = "preds"`, in a permanent table of the connection, for use outside the package (join it to other data on a key column such as `id`). Printing a column shows its SQL and its connection, and reads no row.

```r
beta <- solve(crossprod(D), crossprod(D, "y"))
p <- D %*% beta
p
r <- design_column(D, "y") - p
head(as.data.frame(r, keep = "id", nm = "resid"))
r2 <- col_compute(r^2)
sum(r2)
## a permanent table, for use outside the package
rt <- col_compute(r, name = "resid", nm = "resid")
head(DBI::dbGetQuery(con, "SELECT id, resid FROM resid ORDER BY id"))
```

A column is valid while its connection is open: the views and temporary tables it reads belong to the session of the connection. To keep values, store them in a permanent table with `col_compute(r, name = )`, or read them with `as.data.frame()` and write them where they belong.

#calepin.chunk(echo: false, results: "hide")[
```r
DBI::dbDisconnect(con)
```
]
