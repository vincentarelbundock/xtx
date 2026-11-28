#import "/.calepin/calepin.typ" as calepin
#show: calepin.document

#set document(title: "Summaries and arithmetic of columns")
#metadata((
  title: "Summary.xtx_column",
  pagetitle: "Summaries and arithmetic of columns — column_arithmetic",
  aliases: ("Summary.xtx_column", "cbind.xtx_column", "col_const", "col_ifelse", "col_pmax", "col_pmin", "col_sums", "column_arithmetic", "mean.xtx_column", "rowSums,xtx_column-method"),
  "source-files": ("R/columns_public.R",),
)) <website-metadata>

#calepin.setup(echo: true, eval: true, results: "verbatim", warning: false, message: false)

#calepin.chunk(echo: false, results: "hide")[
```r
library(xtx)
```
]

= Summaries and arithmetic of columns

== Description

Base R functions work on the columns of a design \(see #link("design_column.typ")[`design_column`]\): in memory a column is a numeric vector, and remotely a `design_column` has methods for the same functions, which compute in the database and return the same values as base R on the vector\. Together with `col_ifelse()`, `col_pmin()`, `col_pmax()` and `col_sums()`, which cover what base R cannot dispatch on, the same code runs on a data frame and on a database table\.

== Usage

```r
#| eval: false
# S3 method for class 'xtx_column'
Summary(..., na.rm = FALSE)

# S3 method for class 'xtx_column'
mean(x, trim = 0, na.rm = FALSE, ...)

col_ifelse(test, yes, no, missing = NULL)

col_pmin(a, b, checked = TRUE)

col_pmax(a, b, checked = TRUE)

col_sums(..., na.rm = FALSE)

col_const(like, value)

# S3 method for class 'xtx_column'
cbind(..., deparse.level = 1)

# S4 method for signature 'xtx_column'
rowSums(x, na.rm = FALSE, dims = 1, ...)
```

== Arguments

#table(
  columns: 2,
  table.header([Argument], [Description]),
  [`...`], [For `col_sums()`, columns or numbers, optionally named\. For `cbind()`, columns or numbers\. For the methods, as in base R\.],
  [`na.rm`], [For `col_sums()`, drop `NA` and `NaN` values \(as `sum()`\): `TRUE` or `FALSE`, or one value per argument \(names may change before the first release\)\. For `rowSums()`, only the default\.],
  [`x`], [A remote column\.],
  [`trim`], [Only `0` is supported for a remote column\.],
  [`test`], [A logical column \(in memory a logical vector\)\.],
  [`yes, no, a, b`], [Columns \(vectors in memory, `design_column`s remotely\) or numbers\.],
  [`missing`], [For `col_ifelse()`: `NULL` \(default: a missing `test` gives a missing value, as #link("https://rdrr.io/r/base/ifelse.html")[`ifelse()`]\), or one number, the value where `test` is missing\.],
  [`checked`], [For `col_pmin()` and `col_pmax()`: `TRUE` \(default\) follows #link("https://rdrr.io/r/base/pmin.html")[`pmin()`] and #link("https://rdrr.io/r/base/pmax.html")[`pmax()`] on missing values \(`NA` or `NaN` if either operand is\)\. `FALSE` is for operands known to have no missing values: remotely plain `least()` and `greatest()`, whose SQL is shorter; in memory the same values\.],
  [`like`], [A column, or a design: the rows of the constant column\.],
  [`value`], [One number\.],
  [`deparse.level`], [Not used\.],
  [`dims`], [For `rowSums()`, only the default\.],
)

== Value

`col_ifelse()`, `col_pmin()`, `col_pmax()`: a column \(a double vector in memory, a `design_column` remotely\)\. `col_sums()`: a double vector\. The `Summary` methods and `mean()`: a number, two for `range()`, a logical value for `any()` and `all()`\.

== Base generics on remote columns

- Arithmetic and comparisons \(`+`, `-`, `*`, `/`, `^`, `==`, `<`, \.\.\., `&`, `|`, `!`\), #link("https://rdrr.io/r/base/is.na.html")[`is.na()`], #link("https://rdrr.io/r/base/is.nan.html")[`is.nan()`], #link("https://rdrr.io/r/base/is.finite.html")[`is.finite()`], #link("https://rdrr.io/r/base/is.infinite.html")[`is.infinite()`] and the `Math` functions `exp()`, `log()` \(`log(x, base)`: `base` a number or a column of the design, R's `log(x, 1)` and negative bases included\), `sqrt()`, `abs()`, `log1p()`, `expm1()`, `log2()`, `log10()`, `cos()`, `sin()`, `tan()`, `floor()`, `ceiling()`, `sign()`, `trunc()`, `round(x, digits)` \(R's algorithm: `round(2.5)` is 2, `round(0.15, 1)` is 0\.1\), `signif(x, digits)`, `gamma()` and `lgamma()` build a new column, lazily \(no query\)\.
- `sum()`, `prod()`, `min()`, `max()`, `range()` \(with `finite`\), `any()`, `all()` \(the `Summary` group\) and `mean()` return numbers \(or `TRUE`/`FALSE`/`NA`\), over all the values of all their arguments, as base R \(a column with several columns counts all its values; numbers among the arguments count too\)\. Each call is one query, whatever the number of arguments\. `na.rm = TRUE` drops `NA` and `NaN`\.

== Missing values

SQL `NULL` is `NA` and IEEE `NaN` is `NaN`, with R's rules: arithmetic propagates both, comparisons with either give `NA`, `is.na()` is `TRUE` for both and `is.nan()` for `NaN` only\. `sum()`, `prod()` and `mean()` give `NA` when a value is `NA`, else `NaN` when one is `NaN` \(where R's result for a mix of `NA` and `NaN` depends on the platform, the remote one is `NA`\)\. `min()`, `max()` and `range()` give `NA` when a value is `NA`, else `NaN` \(as base R: an `NA` trumps a `NaN`\), and with no value left, `Inf` or `-Inf` with base R's warning\. `any()` and `all()` are three-valued, a `NaN` counting as `NA`\. `col_ifelse()` gives `NA` where the test is `NA`; `col_pmin()` and `col_pmax()` give `NA` or `NaN` where either value is, as `pmin()` and `pmax()`\.

Remote values are doubles: a logical remote column \(a comparison\) sums to a double, where base R gives an integer, and `as.vector()` returns it as 0 and 1\. Sums over many rows can differ from R's in the last bits \(the database adds in another order, and R in extended precision\); `gamma()` and `lgamma()` use the C library, R its own code \(equal to rounding\)\.

== Helpers

- `col_ifelse(test, yes, no)`: `yes` where `test` is `TRUE`, `no` where it is `FALSE`, `NA` where it is `NA`; `yes` and `no` are columns or numbers\. In memory as #link("https://rdrr.io/r/base/ifelse.html")[`ifelse()`], but always a double vector, without the attributes of `test`\.
- `col_pmin(a, b)`, `col_pmax(a, b)`: the parallel minimum and maximum of two columns, or of a column and a number, as #link("https://rdrr.io/r/base/pmin.html")[`pmin()`] and #link("https://rdrr.io/r/base/pmax.html")[`pmax()`]\.
- `col_sums(...)`: the sums of several columns \(each over all its values, as `sum()`\), in one query remotely: a named double vector, one value per argument, named after the argument \(`col_sums(n = 1 + 0 * r, rss = r^2)`\), or after the variable for a bare name\. `na.rm` can differ by argument: `col_sums(dev, col_ifelse(cond, 1, 0), na.rm = c(FALSE, TRUE))` counts the rows where `cond` is `TRUE` \(`NA` counting as `FALSE`\) in the query that sums `dev`\.
- `col_const(like, value)`: a constant column with the rows of `like` \(a column or a design\): `rep(value, n)` in memory \(names may change before the first release\)\.
- `cbind(...)`: columns side by side, one column with several columns \(numbers recycle\), as `cbind()` of vectors in memory; `x[, j]` selects them again\.
- `rowSums(x)`: the sums across the columns of a column with several columns, a column\.
- `col_compute(x, name = NULL)` \(see #link("design_column.typ")[`design_column`]\): the column stored once in a table, a temporary one or, with a `name`, a permanent table of the connection, and a column that reads it, for a column that later queries read several times\. In memory it returns `x`\.

Arguments can mix remote columns and numbers\. Remote columns must come from the same design \(columns of two designs are an error\); in memory the vectors recycle as in base R\. The `Summary` functions dispatch on their first argument: put a column first \(`sum(r, 1)`, not `sum(1, r)`\)\.

== Examples

```r
d <- data.frame(y = c(1, 4, 6, 2, 8, 3), x = c(0, 1, 1, 2, 3, 3))
D <- design_matrix(y ~ x, d)
y <- design_column(D, "y")
z <- col_ifelse(y > 5, NA, y)
sum(z)
sum(z, na.rm = TRUE)
range(z, na.rm = TRUE)
mean(round(y / 4, 1))
col_sums(n = !is.na(z), total = col_ifelse(is.na(z), 0, z))
col_pmax(y, 3)

if (requireNamespace("duckdb", quietly = TRUE)) {
    con <- DBI::dbConnect(duckdb::duckdb())
    DBI::dbWriteTable(con, "d", d)
    Dr <- design_matrix(y ~ x, con, table = "d")
    yr <- design_column(Dr, "y")
    zr <- col_ifelse(yr > 5, NA, yr)
    sum(zr)
    sum(zr, na.rm = TRUE)
    range(zr, na.rm = TRUE)
    mean(round(yr / 4, 1))
    col_sums(n = !is.na(zr), total = col_ifelse(is.na(zr), 0, zr))
    sort(as.vector(col_pmax(yr, 3)))
    DBI::dbDisconnect(con)
}
```
