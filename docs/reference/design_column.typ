#import "/.calepin/calepin.typ" as calepin
#show: calepin.document

#set document(title: "Columns of a design")
#metadata((
  title: "[.xtx_column",
  pagetitle: "Columns of a design — design_column",
  aliases: ("[.xtx_column", "as.data.frame.xtx_column", "as.matrix.xtx_column", "as.vector.xtx_column", "col_compute", "design_column", "format.xtx_column", "print.xtx_column"),
  "source-files": ("R/column.R",),
)) <website-metadata>

#calepin.setup(echo: true, eval: true, results: "verbatim", warning: false, message: false)

#calepin.chunk(echo: false, results: "hide")[
```r
library(xtx)
```
]

= Columns of a design

== Description

A #emph[column] holds one value per row of a design\. In memory it is a plain numeric vector\. For a remote design \(see #link("design_remote.typ")[`design_remote`]\), it is a `design_column`: a lazy SQL expression over the rows of the design, which nothing evaluates until its values are asked for\. Products with a remote design return one, arithmetic on it builds a new one, and the base generics #link("https://rdrr.io/r/base/as.vector.html")[`as.vector()`], #link("https://rdrr.io/r/base/as.matrix.html")[`as.matrix()`] and #link("https://rdrr.io/r/base/as.data.frame.html")[`as.data.frame()`] evaluate it, as they would a vector in memory\. So the same code runs on a data frame and on a DuckDB table:

== Usage

```r
#| eval: false
design_column(x, column)

# S3 method for class 'xtx_column'
as.vector(x, mode = "any")

# S3 method for class 'xtx_column'
as.matrix(x, ...)

col_compute(x, name = NULL, overwrite = FALSE, nm = NULL, table = TRUE)

# S3 method for class 'xtx_column'
as.data.frame(
  x,
  row.names = NULL,
  optional = FALSE,
  ...,
  keep = NULL,
  nm = NULL
)

# S3 method for class 'xtx_column'
x[i, j, ..., drop = TRUE]

# S3 method for class 'xtx_column'
format(x, ...)

# S3 method for class 'xtx_column'
print(x, ...)
```

== Arguments

#table(
  columns: 2,
  table.header([Argument], [Description]),
  [`x`], [For `design_column()`, a design: a `design_matrix` or a remote design\. For the other functions, a column\.],
  [`column`], [The name of a column \(a string\), or for a remote design a SQL expression over its columns \(`DBI::SQL()`\)\.],
  [`mode`], [As in #link("https://rdrr.io/r/base/as.vector.html")[`as.vector()`]\.],
  [`...`], [Not used\.],
  [`name`], [For `col_compute()`: `NULL` \(default\) for a temporary table, or the name of a permanent table to create on the connection \(a string, quoted as an identifier\)\.],
  [`overwrite`], [For `col_compute()` with a `name`: replace a table of that name if it exists \(otherwise that is an error\)\.],
  [`nm`], [Names of the value columns \(one per column of `x`\): of the data frame for `as.data.frame()`, of the table for `col_compute()` with a `name`\.],
  [`table`], [For `col_compute()`: `TRUE` \(default\) to store the values in a table, `FALSE` for a projection \(see #emph[Getting values]; not with a `name`\)\.],
  [`row.names, optional`], [Not used\.],
  [`keep`], [Names of columns of the remote data to read with the values \(for instance a key, to match the rows with other data\)\.],
  [`i`], [Not supported for a remote column \(remote rows have no order\)\.],
  [`j`], [Indices \(numeric or logical\) of columns of `x`\.],
  [`drop`], [Not used\.],
)

== Value

`design_column()`: a double vector in memory, a `design_column` remotely\. `as.vector()`: a double vector\. `as.matrix()`: a double matrix\. `as.data.frame()`: a data frame\. `col_compute()`: a `design_column` that reads the table \(in memory, `x`\)\.

== Details

```
eta <- D %*% beta                 # vector in memory, design_column remotely
mu  <- exp(eta)
r   <- design_column(D, "y") - mu
as.vector(r)                      # the values, in R
```

== Making columns

- `design_column(D, "y")`: the column `y` of the remote data, as a double \(numeric and logical columns; a logical column gives 0 and 1\)\. `design_column(D, DBI::SQL("y / w"))`: a SQL expression over the columns of the remote data, kept as SQL computes it \(a logical expression stays logical, and can select rows: `D[col, ]`\)\. The column is over the rows of `D` \(its filter, for a design with `na.action = na.omit` or a subset\), not its row scale\.
- In memory, `design_column(D, "y")` returns the variable `y` of the model frame of the design \(the response and the variables of the formula, by their names in the model frame, such as `"log(z)"`\) as a double vector\. A design in memory does not keep the other columns of its data, so any other name is an error: use the vector itself\. A `DBI::SQL()` expression is an error in memory \(there is no database to evaluate it\)\.
- `D %*% b`, `rowSums(D)`, `rowMeans(D)` and `quadratic_diagonal(D, V)` of a remote design return a column \(`b` with `k` columns gives a column with `k` columns\)\.
- Arithmetic: `+`, `-`, `*`, `/`, `^` with numbers and other columns of the same design, comparisons, `&`, `|`, `!`, #link("https://rdrr.io/r/base/is.na.html")[`is.na()`], #link("https://rdrr.io/r/base/is.nan.html")[`is.nan()`], #link("https://rdrr.io/r/base/is.finite.html")[`is.finite()`], #link("https://rdrr.io/r/base/is.infinite.html")[`is.infinite()`] and the functions of the `Math` group \(`exp()`, `log()` with a base \(a number or a column\), `sqrt()`, `abs()`, \.\.\.\), with R's missing values \(SQL `NULL` is `NA`, and `NaN` stays `NaN`\)\. `col_ifelse()`, `col_pmin()` and `col_pmax()` combine columns, and `sum()`, `mean()`, `min()`, `range()`, `any()`, \.\.\. and `col_sums()` summarize them, in one query: see #link("column_arithmetic.typ")[`column_arithmetic`]\.
- `e * D` \(or `D / w`\), with `e` a column of a remote design `D`: the design with rows scaled by the column, as for a vector in memory \(for instance `colSums(e * D, group = "cl")`, or `crossprod(D, w * D)`\)\.

== Getting values

- `as.vector(x)`: the values in R, a double vector, for a column with one column \(a column with several columns is an error: use `as.matrix()`\)\. `as.matrix(x)`: an `n x k` double matrix, without dimnames \(as the products in memory\), for one column or several\. Remotely both read every row: only for results that fit in memory\. In memory they are base R's functions on the vector or matrix\.
- `as.data.frame(x, keep = "id")`: a data frame with one column per column of `x` \(named `.fitted` for `D %*% b`, `.rowsum`, `.rowmean`, `.quad`, the column name for `design_column(D, "y")`, else `.value`; or `nm`\), after the columns `keep` of the remote data, read in the same query\. Remote rows have no order: to match the values with other data, keep a key column\. Two calls to `as.vector()` need not return the rows in the same order\.
- `col_compute(x)`: the column computed once into a temporary table, and a column that reads it \(later queries scan the table instead of recomputing the expression and the products and lookups under it\)\. `col_compute(x, name = "preds")`: the same, into a permanent table `preds` of the connection \(`CREATE TABLE ... AS`; an error if the table exists, unless `overwrite = TRUE`, which replaces it\), for use outside the package\. The table holds the rows of the design \(after its filter: `na.action`, subsets\), each with the columns of the data, the values \(named `nm`: by default as `as.data.frame()` names them\) and internal columns of the package \(names starting with `__xtx_`: the products and lookups that the values were computed from, and values stored earlier\), which can be ignored\. The column returned reads that table and stays aligned with the rows of its design, as with a temporary table: the rows keep, next to the values, the columns that other columns of the design read, so a product made earlier or later is computed over the table when it is combined with the stored values\. The table holds the values as they were when it was written: later changes to the source do not reach it\. If the statement that writes the table fails, no table is made and a table of that name is left as it was \(inside a transaction of yours, DuckDB then aborts that transaction, as for any failed statement\)\. To match the values with other data outside the package, join on a key column of the data \(such as `id`\)\. In memory `col_compute()` returns `x` \(there is no table\)\. `col_compute(x, table = FALSE)`: the column computed once in a projection instead of a table: later expressions refer to its values by name instead of repeating its expressions \(an iterative fit keeps its SQL small this way\), and each query that reads them computes them once per row\. Nothing is created in the database and no statement is run\. Columns computed separately this way combine with each other and with the other columns of their design: `col_compute(a, table = FALSE) + col_compute(b, table = FALSE)`\. \(Names may change before the first release\.\)
- `x[, j]`: the columns `j` of a column with several columns \(as for the matrix in memory\), for instance of `D %*% cbind(b1, b2)`\. Rows cannot be selected: remote rows have no order\.
- `dim(x)` is `c(n, k)`: the rows of the design and the number of columns of `x`, so `nrow()`, `ncol()`, `NROW()` and `NCOL()` work, and `length(x)` is `n * k`, as for a matrix in memory \(`n` is the row count kept with the design: no query\)\. `print()` and `format()` show `n x k`, the SQL and the connection, without reading any row\.

== Lifetime

A remote column is valid while its connection is open\. The temporary tables it reads \(the values of lookups, and the tables of `col_compute()`\) live in the session of that connection, and disappear when it is closed; for a Parquet source, the column keeps its private connection open\. To keep values beyond the session, store them in a permanent table with `col_compute(x, name = )` \(on the connection of the design: for a Parquet source, a private in-memory connection, so read them with `as.data.frame()` instead\), or read them with `as.data.frame()` and write them where they belong\. Columns of one design combine with each other, including the results of separate products \(`D %*% b1 + D %*% b2`, or `(D %*% b) / sqrt(quadratic_diagonal(D, V))`\): a product is a SQL expression over the rows of the design, which each query that reads it computes with the other columns of that query, over the same rows\. Reading a column \(with `as.vector()`, `as.data.frame()`, a sum, \.\.\.\) creates nothing in the database and changes no column, so columns combine in the same way whatever was read before\. A column stored with `col_compute()` combines with the other columns of its design, which are then computed over the stored table\. Two columns each stored with `col_compute()`, separately, are the exception: the database has no key to align the rows of two tables, so combine the columns before storing them \(store one column made from both, or one product `D %*% cbind(b1, b2)` and then `[, 1]` and `[, 2]`\), or read them with `as.data.frame(keep = "id")` and match them in R\. Columns of different designs do not combine \(a subset `D[i, ]` is another design, even when it keeps every row, as is a design made again from the same data\)\. A column made by `design_column()` combines with any column of its design\.

== Examples

```r
d <- data.frame(y = rnorm(100), x = rnorm(100), f = gl(4, 25))
D <- design_matrix(y ~ x * f, d)
b <- solve(crossprod(D), crossprod(D, d$y))
r <- design_column(D, "y") - D %*% b
head(as.vector(r))

if (requireNamespace("duckdb", quietly = TRUE)) {
    con <- DBI::dbConnect(duckdb::duckdb())
    d$id <- seq_len(nrow(d))
    DBI::dbWriteTable(con, "d", d)
    Dr <- design_matrix(y ~ x * f, con, table = "d")
    r <- design_column(Dr, "y") - Dr %*% b
    r
    head(as.data.frame(r, keep = "id", nm = "resid"))
    r2 <- col_compute(r^2)
    sum(r2)
    rt <- col_compute(r, name = "resid", nm = "resid")
    head(DBI::dbGetQuery(con, "SELECT id, resid FROM resid ORDER BY id"))
    DBI::dbDisconnect(con)
}
```
