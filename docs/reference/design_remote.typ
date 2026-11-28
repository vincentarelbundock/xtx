#import "/.calepin/calepin.typ" as calepin
#show: calepin.document

#set document(title: "Out-of-core designs on DuckDB")
#metadata((
  title: "coerce,xtx_remote,CsparseMatrix-method",
  pagetitle: "Out-of-core designs on DuckDB — design_remote",
  aliases: ("coerce,xtx_remote,CsparseMatrix-method", "coerce,xtx_remote,dgCMatrix-method", "colMeans,xtx_remote-method", "colSums,xtx_remote-method", "design_remote", "rowMeans,xtx_remote-method", "rowSums,xtx_remote-method", "summary.xtx_remote"),
  "source-files": ("R/remote.R",),
)) <website-metadata>

#calepin.setup(echo: true, eval: true, results: "verbatim", warning: false, message: false)

#calepin.chunk(echo: false, results: "hide")[
```r
library(xtx)
```
]

= Out-of-core designs on DuckDB

== Description

With a DBI connection to DuckDB, or the path of Parquet files, as `data`, #link("design_matrix.typ")[`design_matrix()`] returns a #emph[remote design]: the rows stay in the database, and products with the design run there, so that data larger than memory can be used\.

== Usage

```r
#| eval: false
# S4 method for signature 'xtx_remote'
colSums(x, na.rm = FALSE, dims = 1, ..., group = NULL, reorder = TRUE)

# S4 method for signature 'xtx_remote'
colMeans(x, na.rm = FALSE, dims = 1, ..., group = NULL, reorder = TRUE)

# S4 method for signature 'xtx_remote'
rowSums(x, na.rm = FALSE, dims = 1, ...)

# S4 method for signature 'xtx_remote'
rowMeans(x, na.rm = FALSE, dims = 1, ...)

# S3 method for class 'xtx_remote'
summary(object, sql = FALSE, ...)
```

== Arguments

#table(
  columns: 2,
  table.header([Argument], [Description]),
  [`x, object`], [A remote design\.],
  [`na.rm, dims`], [Only the defaults are supported\.],
  [`...`], [Not used\.],
  [`group`], [For `colSums()` and `colMeans()`, the name of a column of the remote data whose values define the groups, or a SQL expression over its columns \(`DBI::SQL()`\)\.],
  [`reorder`], [Only `TRUE` \(groups sorted\) is supported\.],
  [`sql`], [For `summary()`: also print the generated SQL \(the query that computes the model variables, and the aggregates of `crossprod()`\)\.],
)

== Details

```
design_matrix(formula, con, table = "flights")
design_matrix(formula, con, sql = "SELECT * FROM flights WHERE year > 2000")
design_matrix(formula, "path/to/*.parquet")
design_matrix(spec, con, table = "new_flights")   # fit in R, score in the database
```

== Sources

`con` must be a #link("https://rdrr.io/pkg/DBI/man/dbConnect.html")[`DBI`] connection to DuckDB \(packages DBI and duckdb\)\. Give either `table`, the name of a table or view \(a character string or a #link("https://rdrr.io/pkg/DBI/man/Id.html")[`DBI::Id()`]\), or `sql`, a `SELECT` query\. A character vector of Parquet paths \(globs allowed\) opens a private in-memory DuckDB connection that reads them with `read_parquet()`; it is closed when the design is garbage collected\.

== Training and new data

With a formula, the design is learned from the remote data, as `design_matrix(formula, data)` learns it from a data frame:

- Levels: those of an `ENUM` column \(a factor in R\) are its declared levels; those of a character column are its distinct values in the rows used, sorted as #link("https://rdrr.io/r/base/factor.html")[`factor()`] sorts them in R \(not by SQL collation\)\. Logical columns are factors with levels `FALSE` and `TRUE`\.
- Rows: rows with a missing value \(SQL `NULL`, or `NaN`\) in any model variable, response included, are dropped, as #link("https://rdrr.io/pkg/stats/man/na.omit.html")[`stats::na.omit()`] drops them\.
- Transformations whose parameters depend on the data, such as #link("https://rdrr.io/pkg/stats/man/poly.html")[`stats::poly()`], #link("https://rdrr.io/pkg/splines/man/ns.html")[`splines::ns()`] or #link("https://rdrr.io/r/base/scale.html")[`scale()`], are fitted as #link("https://rdrr.io/pkg/stats/man/model.frame.html")[`stats::model.frame()`] fits them, on the whole of their input columns: those columns \(only those\) are read into R once\. To avoid this, fit the specification in R, on a sample or on other data, and pass it instead of the formula\.
- Variables computed by an R expression that yields a factor, such as `factor(cyl)` or `interaction(a, b)`, take the levels that the expression gives on the distinct values of its input columns\.

With a specification \(a design, a plan, a fitted model or another remote design\), the levels, contrasts and transformation parameters come from it, and the remote data is only scored\. Levels absent from the specification raise the error that #link("https://rdrr.io/pkg/stats/man/model.frame.html")[`stats::model.frame()`] raises, and columns of another type than in training \(a numeric variable stored as text, say\) the error of #link("https://rdrr.io/pkg/stats/man/predict.lm.html")[`stats::predict.lm()`]\. As in memory, the rows with missing values are kept by default \(`na.action = na.pass`\): their predictions are `NULL`, and products involving their columns are `NA`\. With `na.action = na.omit`, they are dropped\.

== Engine

Products with a remote design are computed in DuckDB\. SQL computes factors, character and logical columns, numeric columns, interactions, `factor(x)`, `as.factor(x)` and `as.character(x)` of a text, factor or integer column \(with the levels R gives them: in numeric order for `factor()` of an integer column, sorted as text for `as.character()`\), and numeric expressions built from `+`, `-`, `*`, `/`, `^`, #link("https://rdrr.io/r/base/I.html")[`I()`], #link("https://rdrr.io/r/base/log.html")[`log()`] \(with base 2, 10 or a number\), #link("https://rdrr.io/r/base/log2.html")[`log2()`], #link("https://rdrr.io/r/base/log10.html")[`log10()`], #link("https://rdrr.io/r/base/exp.html")[`exp()`], #link("https://rdrr.io/r/base/sqrt.html")[`sqrt()`] and #link("https://rdrr.io/r/base/abs.html")[`abs()`]\. Any other variable, such as #link("https://rdrr.io/pkg/stats/man/poly.html")[`stats::poly()`], #link("https://rdrr.io/pkg/splines/man/ns.html")[`splines::ns()`], `factor()` of a numeric expression or with `levels`, or a function of yours, is evaluated in R once, in batches of `chunk_size` rows, and written with the columns of the remote data to a temporary table \(one pass over the data, and a copy of it in the database\), whose columns the design then reads like any others; a message names the first such term\. The expression is evaluated on the rows where its input columns are present, and a row where it gives `NA` is a row with a missing value\. With `engine = "sql"`, such a formula is refused instead\. `summary(D)` shows, for each term, whether SQL computes it, and `summary(D, sql = TRUE)` the generated SQL\.

A term whose factors have more than about eight million combined level-column cells is refused; for a factor with many levels, use fixed effects \(`y ~ x | f` in `xtx.fit::lm_xtx()` and `xtx.fit::glm_xtx()`\)\.

== Operations

The operations of a design in memory \(see #link("design_algebra.typ")[`design_algebra`] and #link("design_arithmetic.typ")[`design_arithmetic`]\), computed in the database:

- `crossprod(D)`: `X'X`; `crossprod(D, "y")`: `X'y` for one or more columns of the remote data \(or a `DBI::SQL()` expression\)\.
- Scaled rows: `"w" * D` \(a column name\), `DBI::SQL("1 / v") * D` \(an expression over the columns\), `2 * D`, `D / "w"` and `-D` scale each row, as `w * D` with a numeric vector in memory\. So `crossprod(D, "w" * D)` is `X'WX`, `crossprod("w" * D, "y")` is `X'Wy` and `colSums("w" * D, group = "g")` the weighted sums within groups\. Only one design is involved: `crossprod(x, y)` needs `y` to be `x` with other scales\.
- `colSums(D)`, `colMeans(D)`, and with `group = "g"` \(a column, or a `DBI::SQL()` expression\), the sums or means within groups \(sorted as #link("https://rdrr.io/pkg/base/man/rowsum.html")[`base::rowsum()`] sorts them\)\.
- `D %*% beta`, `rowSums(D)`, `rowMeans(D)` and `quadratic_diagonal(D, V)`: a remote column \(#link("design_column.typ")[`design_column`]\), a lazy SQL expression over the rows of the design, where a design in memory returns a vector\. Arithmetic on it builds another column; `as.vector(p)` reads its values into R, `as.data.frame(p, keep = "id")` reads them with columns of the remote data \(to join them back\), and #link("design_column.typ")[`col_compute()`] stores it in a table \(temporary, or permanent with a `name`\)\. `design_column(D, "y")` makes a column from a column of the remote data\. Columns of one design combine, separate products included \(`D %*% b1 + D %*% b2`\): a product is SQL that each query reading it computes over the rows of the design, with the other columns of that query; reading a column creates nothing in the database\.
- `as.matrix(D)` collects the rows and builds the dense matrix, only when there are at most `getOption("xtx.remote_max_rows", 1e5)` rows or with `force = TRUE`\.
- `design_matrix(D, newdata)` applies the design to a data frame, and `design_matrix(D, con, table = ...)` to other remote data\.
- `summary(D)`: the engine and, for each term, whether SQL computes it\.

== Compression

When the rows have few distinct #emph[design patterns] \(the same factor levels and the same values of the numeric variables of the design, as with factors and discrete numeric variables\), `crossprod()`, `colSums()`, `colMeans()`, `xtx.fit::lm_xtx()` and `xtx.fit::glm_xtx()` can compress: the one scan of the rows also groups by the numeric variables, which returns each pattern with its number of rows and the sums of the weights and responses, and the result is computed from the design of the patterns, in memory\. The option `xtx.remote_compress` controls it\. `TRUE` always compresses \(unless every row is distinct\); `FALSE` never does\. With `"auto"` \(default\), the choice depends on the operation:

- `colSums()` and `colMeans()` are never compressed \(grouping by the numeric variables costs more than it saves\);
- `crossprod()` is compressed only for designs with many numeric columns \(currently at least 8 distinct products of numeric variables; the threshold is provisional\), where compression saves many per-row products, and when the patterns are at most 1% of the rows\. Other designs skip the decision entirely;
- `xtx.fit::glm_xtx()` is compressed when the patterns are at most 20% of the rows \(each IRLS iteration then reads the patterns instead of the rows\);
- `xtx.fit::lm_xtx()` without weights is compressed when the patterns are at most 0\.1% of the rows;
- weighted `xtx.fit::lm_xtx()` is compressed only for designs with many numeric columns \(as `crossprod()`\), when the patterns are at most 0\.1% of the rows \(the sum of squares within the patterns takes a second pass over the rows\)\.

These thresholds are provisional\. The option `xtx.remote_compress_ratio` \(unset by default\), when set, caps every one of these ratios: under `"auto"`, the patterns must then also be at most that fraction of the rows\.

In every case the number of patterns must also be at most `getOption("xtx.remote_compress_max", 1e6)`\. To decide, a design with a continuous numeric column \(a `DOUBLE` column with many distinct values in a sample of the rows\) is not compressed, at almost no cost; otherwise the number of patterns is estimated over all the rows \(`approx_count_distinct()`\)\. Designs that keep rows with missing values \(`na.action = na.pass`\) are not compressed\.

Remote designs keep only the #emph[Algebra] promise of the package: results agree with the in-memory computation up to floating point rounding\. With infinite values in the data, the placement of `NaN` and `Inf` in sums can differ from the dense computation, and `NaN` may appear where the dense result has `NA`\.

== Fits

Linear and generalized linear models on remote data, with or without fixed effects, are in the package 'xtx\.fit': `xtx.fit::lm_xtx()` and `xtx.fit::glm_xtx()`; its help page `remote_fits` describes how they run on the remote products and how they differ from fits in memory\.

== Arguments of remote data

With remote data, #link("design_matrix.typ")[`design_matrix()`], #link("design_frame.typ")[`design_frame()`], `xtx.fit::lm_xtx()`, `xtx.fit::glm_xtx()` and `predict()` take, through `...`: `table` or `sql` \(the remote rows: a table name, or a `SELECT` query; not with Parquet paths\), `engine` \(`"auto"`, or `"sql"` to refuse variables that must be evaluated in R\) and `chunk_size` \(rows per batch when variables are evaluated in R, default `getOption("xtx.chunk_size", 250000)`\)\.

== Examples

```r
if (requireNamespace("duckdb", quietly = TRUE)) {
    con <- DBI::dbConnect(duckdb::duckdb())
    d <- data.frame(y = rnorm(100), x = rnorm(100), f = gl(4, 25))
    DBI::dbWriteTable(con, "d", d)
    D <- design_matrix(y ~ x * f, con, table = "d")
    D
    beta <- solve(crossprod(D), crossprod(D, "y"))
    p <- D %*% beta
    p
    head(as.vector(p))
    r <- design_column(D, "y") - p
    head(as.data.frame(r, nm = "resid"))
    summary(D)
    DBI::dbDisconnect(con)
}
```
