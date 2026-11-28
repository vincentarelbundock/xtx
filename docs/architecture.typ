#import "/.calepin/calepin.typ" as calepin
#show: calepin.document

#set document(title: [Architecture])
#metadata((title: "Architecture")) <website-metadata>

#calepin.setup(echo: true, eval: false)

#title()

This page describes how `xtx` is built, for readers who want to know where a result comes from or where to look in the source. File names refer to the #link("https://github.com/vincentarelbundock/xtx")[repository].

= Overview

A model matrix is highly redundant. A factor with 1,000 levels fills 999 columns, but each row holds one integer code. The package keeps that compact description (factor codes, numeric variables, contrast matrices and the layout of the columns) and computes with it.

Three choices follow:

+ *A design is a description of the matrix.* The dense matrix is built only by `as.matrix()`, and then it is identical to the result of `model.matrix()`.
+ *Products are sums over rows.* $X'X$, $X'y$, column sums and grouped sums are computed from codes and numbers, without the dense matrix.
+ *Each operation has two methods.* One runs in C on data in memory. The other writes SQL for data in DuckDB tables or Parquet files. Code written with these operations runs on both.

= Layers

#table(
  columns: (auto, 1fr, auto),
  [*Layer*], [*Role*], [*Source*],
  [Frame and plan], [R's own model frame, and what new data must reuse: transformations, levels, contrasts], [`R/frame.R`, `R/prep.R`, `R/plan.R`, `R/parts.R`],
  [Design], [Checks, contrasts and column names as in R. The dense fill.], [`R/design_matrix.R`, `src/modelmatrix.c`],
  [Kernels], [Products with a design in memory, from codes and numbers], [`R/algebra.R`, `R/structured.R`, `src/structured.c`, `src/kernels.c`, `src/groups.c`],
  [Primitives], [The operations on designs and columns, with a method in memory and a remote method], [`R/primitives.R`, `R/primitives_lookup.R`, `R/primitives_fused.R`],
  [Columns and groupings], [The public column type and the grouping functions], [`R/column.R`, `R/columns_public.R`, `R/groups_public.R`, `R/group_codes.R`, `src/fe.c`],
  [Remote], [Designs over DuckDB tables, queries and Parquet files. The SQL engine.], [`R/remote.R`, `R/remote_sql.R`, `R/remote_compress.R`, `R/remote_resolve.R`],
  [Least squares], [`qr()` of a design, and exact centering], [`R/qr.R`, `R/center.R`],
)

= Designs

*Preparation.* `model.matrix()` and `model.frame()` do a good deal of work in R before they reach C: non-standard evaluation, character columns turned into factors, contrasts, `xlev`, `na.action`. The package does not copy this code. At load time it takes the installed `stats:::model.matrix.default()` and makes it return what its call to C would have received (`R/prep.R`). `model.frame.default()` is patched in the same way (`R/frame.R`). If the installed code is not the code a patch was written for, the function of stats is used instead.

*Setup and fill.* R's C routine `modelmatrix()` is ported and split in two (`src/modelmatrix.c`). `xtx_setup()` performs every check, evaluates the contrasts, emits the warnings and computes the column names, in the order of R. `xtx_fill()` allocates the matrix and fills it. A design is the result of the setup plus the variables.

The fill performs the same floating point operations as R, cell by cell. Only the order in which rows are visited changes. This is why `as.matrix(D)` is identical to `model.matrix()`, including `NA`, `NaN`, `Inf` and the sign of zero.

*Plans.* A plan keeps what a design fixes for later data: the terms, the transformations that depend on the data (`poly()`, `ns()`, `scale()`) with their training parameters, the factor levels and the contrasts (`R/plan.R`). `design_matrix(plan, newdata)` applies it to new rows. Remote designs carry a plan too. `design_parts()` is a generic that extracts these parts from a fitted model (`R/parts.R`); see the #link("recipes/package-developers.html")[recipes for package developers].

= Products

Take term $t$ of the design. On row $i$, let $F_t (i)$ be the combined level of its factors, $N_t (i)$ the product of its numeric variables, and $C_t$ the Kronecker product of the contrast matrices of its factors, with one row per combined level. The columns of the term are

$ X_t [i, dot] = N_t (i) thin C_t [F_t (i), dot]. $

A product with the design then separates in two parts: a sum over rows, which involves only codes and numbers, and a small dense product with the $C_t$. For two terms $s$ and $t$,

$ (X' W X)_(s t) = C_s ' M_(s t) C_t, quad M_(s t) [a, b] = sum_(i : F_s (i) = a, thin F_t (i) = b) w_i N_s (i) N_t (i). $

$M_(s t)$ is a weighted crosstab. In memory it is one pass over the rows in C. In a database it is a `GROUP BY` over the factors of the two terms. $X' W y$ and column sums have the same form with one term.

In memory there are two sets of kernels:

#table(
  columns: (auto, 1fr, auto),
  [*Kernels*], [*Used when*], [*Source*],
  [Structured], [Every input is finite and no product can overflow. They use the factorization above, and agree with the dense computation up to rounding.], [`src/structured.c`],
  [Row], [Every other case, Jacobians that combine the rows of two designs, and results with one value per row such as $X beta$. They list the columns of each row that can be non-zero, with values computed in the order of the dense fill, so that `0 * Inf` is `NaN` where the dense matrix has it.], [`src/kernels.c`],
)

*Repeated rows.* Rows with the same factor codes and the same numeric values give the same row of the design. When a design has far fewer such patterns than rows, the kernels compute once per pattern (`src/groups.c`, `R/algebra.R`, option `xtx.dedup`). `design_patterns()` returns the patterns with the sufficient statistics of a response (`R/exports.R`).

*Threads.* The dense fill and the kernels can run on several threads (`options(xtx.threads = )`). The default is half of the available threads.

The functions are listed in #link("operations.html#algebra")[Operations], and timed in #link("benchmarks.html")[Benchmarks].

= Remote data

A remote design is made from a DuckDB table, a SQL query or Parquet files. The rows stay in the database. See #link("out-of-core.html")[Out-of-core].

== Primitives

`R/primitives.R` defines a small set of operations on a *design* and on *columns*. A column holds one value per row of a design.

#table(
  columns: (auto, 1fr, 1fr),
  [], [*In memory*], [*Remote*],
  [Design], [Class `xtx_matrix`], [Class `xtx_remote`: a source, a plan, and how to compute each variable in SQL],
  [Column], [Numeric vector or matrix], [Class `xtx_column`, made by `design_column()`: SQL expressions over the rows of the design],
  [Method], [C kernels], [SQL sent to DuckDB],
)

The operations are:

- products with the design: $X'X$, $X'c$, $X B$ and $op("rowSums")((X V) circle.small X)$;
- sums and grouped sums of columns, and crosstabs of two groupings;
- lookups of one value per group, and subsets of rows;
- arithmetic on columns: `+`, `*`, comparisons, `exp()`, `log()` and other functions.

Missing values follow R: SQL `NULL` is `NA`, and a sum returns `NA` where `sum()` does.

== Columns

The remote column is lazy. Arithmetic on it builds SQL, and nothing runs until values or a sum are asked for. `D %*% b` returns a column on a remote design and a vector in memory, so the same expression works on both.

An operation that reads rows takes all its columns at once and sends one statement. For example, `sum()` of several columns is one query. Columns combine when they come from the same design.

`as.vector()`, `as.matrix()` and `as.data.frame()` read the values into R. `col_compute()` stores a column in a table of the database. Base generics such as `sum()`, `mean()`, `range()` and `round()` have methods for columns, and what base R cannot dispatch on is exported as `col_ifelse()`, `col_pmin()`, `col_pmax()` and `col_sums()` (`R/column.R`, `R/columns_public.R`). See #link("operations.html#columns")[Columns].

== Groupings

A grouping assigns an integer code to each row (`R/groups_public.R`). In memory the codes come from one pass in C. Remotely the distinct values of the key are read, sorted in R and written to a temporary table, and a join turns them into codes.

#table(
  columns: (auto, 1fr),
  [*Function*], [*Result*],
  [`design_groups(D, key)`], [The codes, the number of groups and the levels],
  [`group_sums(col, g)`], [Sums of columns within groups: one `GROUP BY` remotely, one pass in C in memory],
  [`group_crosstab(g1, g2, w)`], [A sparse matrix of weighted counts for each pair of groups],
  [`group_lookup(values, g)`], [The column `values[codes]`],
  [`group_nested(g1, g2)`], [Whether each group of `g1` lies within one group of `g2`],
)

See #link("operations.html#groupings")[Groupings].

== SQL engine

For a remote design, the sums of a term, or of a pair of terms for $X' W X$, are `GROUP BY` aggregates over the factor columns as stored (`R/remote_sql.R`). All of them come from one statement (`GROUP BY GROUPING SETS`). R then turns levels into codes and applies the contrasts $C_t$ to the sums that come back.

What comes back is bounded by the cells of the terms involved. It does not grow with the number of rows.

A variable that SQL cannot compute, such as a call to an R function with no translation, is evaluated in R once, batch by batch, into a column of a temporary table. It is then a column like any other.

== Compression

When the data have few patterns, the scan also groups by the numeric variables of the design. It returns each cell with its number of rows and the sums of the weights and responses, and the products are computed in memory from the cells (`R/remote_compress.R`).

An estimate of the number of patterns, from a sample and from `approx_count_distinct()`, decides whether to try. When there are too many cells, the ordinary scan runs instead.

= Least squares

`qr()` of a design decomposes $X P = Q R$ without the dense $X$ or the $n times p$ matrix $Q$ (`R/qr.R`).

- $R$ is the Cholesky factor of $X'X$, from `crossprod()`.
- Columns are pivoted with the rule of `lm()`. Columns close to the tolerance are checked against the data, because $X'X$ alone is not accurate enough to judge them.
- Coefficients solve the normal equations. They are then refined from the residuals $y - X b$. Each step costs one $X b$ and one $X' r$.
- Numeric variables far from zero are centered when that is an exact linear transform of the design (`R/center.R`). This keeps $X'X$ accurate. Remote designs are centered whenever possible. Designs in memory are centered only when $X'X$ is badly conditioned.

The cost is a few passes over the rows, and memory of order $p^2$. The result works with `qr.coef()`, `qr.fitted()`, `qr.resid()` and `solve()`.

= Costs

- Only the dense matrix is exact. Products and `qr()` agree with the dense computation up to rounding, not bit for bit.
- Least squares goes through $X'X$. On badly conditioned designs the coefficients agree with a Householder QR to fewer digits. In memory, `qr()` falls back to the dense decomposition when $X'X$ cannot represent the scale of the columns.
- Memory is of order $p^2$ for $X'X$. A design with very many columns is not the intended use.
- Remote data must be in DuckDB, or in Parquet files that DuckDB reads.

The companion package `xtx.fit` fits models with the exported functions of `xtx` only: see #link("fit/lm_xtx.html")[`lm_xtx()`] and #link("fit/fixed_effects.html")[fixed effects].
