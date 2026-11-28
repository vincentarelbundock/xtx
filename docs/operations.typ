#import "/.calepin/calepin.typ" as calepin
#show: calepin.document

#set document(title: [Operations])
#metadata((title: "Operations")) <website-metadata>


#calepin.setup(echo: true, eval: true, results: "verbatim", warning: false, message: false)

#calepin.chunk(echo: false, results: "hide")[
```r
options(width = 200)
```
]

#title()

This page lists what you can do with a design matrix, with columns and with groupings. The #link("recipes/algebra.html")[Algebra recipes] check each result against base R, and #link("benchmarks.html")[Benchmarks] times the operations.

= Building designs

A design matrix holds the prepared variables and the layout of the columns, but not the dense matrix. `as.matrix()` builds the latter, which defers the cost of allocating and filling it.

#table(
  columns: (auto, 1fr, auto),
  [*Operation*], [*What it computes*], [*Usage*],
  [`design_matrix()`], [A design from a formula and data, without the dense matrix. Same columns as `model.matrix()`.], [`design_matrix(y ~ x * f, d)`],
  [`as.matrix()`], [The dense matrix, identical to `model.matrix()`.], [`as.matrix(D)`],
  [Threads], [Writing the dense output can use several threads. The option takes an integer, `0` for all available threads, or a fraction of them; the default is half.], [`options(xtx.threads = 4)`],
  [New data], [The design of new rows, encoded as the training rows were: same `poly()` basis, same factor levels.], [`design_matrix(D, nd)`],
  [Fitted model], [The same, from a fitted model.], [`design_matrix(fit, newdata = nd)`],
  [Specification], [Prepare once, then reuse for each batch of new data.], [`spec <- design_matrix(fit, newdata = NULL)`],
  [Sparse matrix], [A `CsparseMatrix` without forming the dense matrix, identical to converting the dense one.], [`as(D, "CsparseMatrix")`],
  [Repeated rows], [Computes once per distinct row and reuses the result. Row-wise results stay identical bit for bit. The option is `"auto"` (default), `TRUE` or `FALSE`.], [`options(xtx.dedup = "auto")`],
)

Examples and benchmarks of threads, #link("recipes/algebra.html#new-data")[reuse on new data], #link("recipes/algebra.html#sparse-matrices")[sparse matrices] and #link("recipes/algebra.html#repeated-rows")[repeated rows] are in the #link("recipes/algebra.html")[Algebra recipes].

= Algebra

None of these forms the dense matrix.

#table(
  columns: (auto, 1fr, auto),
  [*Operation*], [*What it computes*], [*Usage*],
  [Arithmetic], [Scales a design by a number, or each row by a vector, and adds or subtracts two designs with matching rows and columns. The result is a design matrix.], [`c(1, 2) * D`],
  [`%*%`], [`X B` from the factor codes and numeric variables. `B` can have several columns.], [`D %*% beta`],
  [`crossprod()`], [`X'X`, as an ordinary matrix that goes straight to `chol()`, `chol2inv()`, `eigen()` or `solve()`.], [`crossprod(D)`],
  [`crossprod()` with a vector], [`X'W`, such as `X'y` in least squares or weighted sums of the columns.], [`crossprod(D, w)`],
  [`crossprod()` weighted], [The weighted `X'WX`.], [`crossprod(D, w * D)`],
  [`qr()`], [A QR decomposition without forming `X` or `Q`. `R` is the Cholesky factor of `crossprod(D)`, pivoted as in `qr()` and `lm()`. Coefficients are refined from the residuals. Agrees with `qr(as.matrix(D))` up to rounding; signs of `R` and `Q` can differ, and `qr.Q()` and `qr.qty()` return only the first `rank` columns or rows. See #link("reference/design_qr.html")[`design_qr`].], [`qr(D)`],
  [`qr` methods], [The usual functions accept the result: `qr.coef()`, `qr.fitted()`, `qr.resid()`, `qr.qty()`, `qr.qy()`, `qr.Q()`, `qr.R()` and `solve()`.], [`qr.coef(qr(D), y)`],
  [`colSums()`, `colMeans()`], [Column sums and means. With one group label per row, one row per group, as with `base::rowsum()`. Group sums are a building block of cluster-robust standard errors and subgroup averages. Multiply by a vector first to weight the rows.], [`colSums(D, group = g)`],
  [`rowSums()`, `rowMeans()`], [Row sums and means.], [`rowSums(D)`],
  [`quadratic_diagonal()`], [The diagonal of `X V X'` without forming `X` or the `n` by `n` product. With `V` the variance-covariance matrix of the coefficients it holds the variance of each fitted value; with `V = solve(crossprod(X))`, the leverages.], [`quadratic_diagonal(D, V)`],
)

The #link("recipes/algebra.html")[Algebra recipes] check each of these against base R. The #link("fit/lm_xtx.html#least-squares")[`lm_xtx()` page] benchmarks `qr()` against other least squares implementations.

= Columns

A column holds one value per row. In memory it is a numeric vector. For a remote design (a DuckDB table or Parquet files) it is a `design_column`: a lazy SQL expression over the rows, which nothing evaluates until a sum, a product or `as.vector()` asks for its values. Only small results (numbers, $p times p$ matrices, one row per cluster) come back to R. An estimator written with these functions and the products above runs unchanged on a data frame and on a database table: see #link("recipes/package-developers.html")[Estimators] for worked examples.

#table(
  columns: 2,
  [*Task*], [*Functions*],
  [Columns], [`design_column(D, "y")`, arithmetic, comparisons, `exp()`, `log()`, `round()`, `lgamma()`, `is.na()`, ...],
  [Combining columns], [`col_ifelse()`, `col_pmin()`, `col_pmax()`],
  [Summaries], [`sum()`, `mean()`, `min()`, `max()`, `range()`, `any()`, `all()`, `col_sums()` (several sums, one query)],
  [Values], [`as.vector()`, `as.matrix()`, `as.data.frame(keep = )`, `col_compute()` (a temporary table, or a permanent one with `name =`)],
  [Compression], [`design_patterns()` (a design and a response compressed to sufficient statistics)],
)

== Lazy columns

A remote column is a SQL expression. These points hold for it, and the in-memory vector gives the same results.

#table(
  columns: 2,
  [*Topic*], [*Behaviour*],
  [Printing], [Shows the dimensions ($n times k$, as `dim()`) and the SQL, without reading a row.],
  [Arithmetic], [Builds a new expression. The values are computed only when a summary, a product or `as.vector()` needs them.],
  [Missing values], [Follow R: SQL `NULL` is `NA` and `NaN` stays `NaN`. `is.na()` is true for both. `sum()`, `mean()`, `min()` and the other summaries give `NA` where base R does, and `na.rm = TRUE` drops both.],
  [Combining columns], [Columns of one design combine with each other, also when they come from separate products, read or not. A product read on its own is made again, as one more view, where it is combined.],
  [What does not combine], [Columns of two designs, and two columns stored separately with `col_compute()`.],
  [Validity], [A column is valid while its connection is open.],
)

The #link("recipes/package-developers.html#lazy-columns")[lazy columns] section of the recipes shows printing, missing values and `as.data.frame(keep = )` on a remote design.

= Groupings

Fixed-effect estimators work from groupings: each row gets the number $1..G$ of its group, from the values of one or more key variables.

#table(
  columns: 3,
  [*Function*], [*Returns*], [*Note*],
  [`design_groups(D, key)`], [A grouping of the rows of `D` by `key`], [In memory the key is given by its values, remotely by the name of a column.],
  [`group_sums(col, g)`], [A $G times k$ matrix of the sums of a column (or several) within each group], [One query remotely. `group_sums(1, g)` counts the rows of each group.],
  [`group_crosstab(g1, g2, w)`], [A sparse $G_1 times G_2$ matrix], [The (weighted) number of rows in each pair of groups.],
  [`group_lookup(values, g)`], [The column whose value in a row is the value of its group, `values[codes]`], [With `group_sums()`, it demeans a column by group: $tilde(y) = y - "group_lookup"("group_sums"(y, g) \/ n_g, g)$.],
  [`group_nested(g1, g2)`], [`TRUE` or `FALSE`], [Whether each group of `g1` lies within one group of `g2`.],
)

A grouping is a list with these elements:

#table(
  columns: 2,
  [*Element*], [*Contents*],
  [`codes`], [The group of each row: an integer vector in memory, a `design_column` remotely.],
  [`G`], [The number of groups.],
  [`levels`], [A data frame of the key values of each group, in code order, sorted as the levels of fixed effects.],
  [`n`], [The number of rows.],
  [`design`], [The design the grouping was made from.],
)

Results with one row per group are in code order: row `l` is `g$levels[l, ]`. Every row must have a group: a missing key (`NA`, `NaN`, SQL `NULL`) is an error on both backends, so select the rows with a key first (remotely `D[DBI::SQL("cl IS NOT NULL"), ]`). The #link("recipes/package-developers.html#the-within-transformation")[within transformation] is a worked example.
