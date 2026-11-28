#import "/.calepin/calepin.typ" as calepin
#show: calepin.document

#set document(title: "Compress a design and a response to sufficient statistics")
#metadata((
  title: "design_patterns",
  pagetitle: "Compress a design and a response to sufficient statistics — design_patterns",
  aliases: ("design_patterns",),
  "source-files": ("R/exports.R",),
)) <website-metadata>

#calepin.setup(echo: true, eval: true, results: "verbatim", warning: false, message: false)

#calepin.chunk(echo: false, results: "hide")[
```r
library(xtx)
```
]

= Compress a design and a response to sufficient statistics

== Description

`design_patterns()` groups the rows of a design that have the same row of the model matrix \(and, with `keys`, the same values of the key columns\) into patterns, and summarizes a response `z` within each pattern: the sum of the weights, the weighted mean, and the weighted sum of squares around that mean\. A weighted least squares fit of the pattern means on the pattern design, with the summed weights as weights, has the same coefficients as the fit on every row; its residual sum of squares plus `within` is the residual sum of squares of the fit on every row\. This is the compression `xtx.fit::lm_xtx()` uses when a design has few distinct rows \(factors and discrete numeric variables\)\.

== Usage

```r
#| eval: false
design_patterns(
  x,
  z,
  w = NULL,
  keys = NULL,
  compress = "auto",
  fit = c("lm", "glm"),
  check = TRUE
)
```

== Arguments

#table(
  columns: 2,
  table.header([Argument], [Description]),
  [`x`], [A design: a `design_matrix` built from a data frame, or a remote design\. Scaled designs \(`s * D`\) are not accepted\.],
  [`z`], [The response\. In memory, a numeric vector with one value per row of `x` \(such as a #link("design_column.typ")[`design_column()`] of `x`, or a product with it\)\. For a remote design, the name of a numeric column of the remote data, a `DBI::SQL()` expression over its columns, or a #link("design_column.typ")[`design_column()`] of the design \(such as `design_column(x, "y") - x %*% b`\)\. Missing or infinite values are an error\.],
  [`w`], [Optional weights, in the same form as `z`: finite and non-negative\. Rows with zero weight do not enter the statistics\.],
  [`keys`], [Optional extra grouping columns: patterns are then also split by their values\. In memory, a numeric vector, matrix or data frame of numeric columns with one row per row of `x`\. For a remote design, the names of numeric columns of the remote data, a `DBI::SQL()` expression, or a #link("design_column.typ")[`design_column()`] of the design\.],
  [`compress`], [`"auto"` \(default\), `TRUE` or `FALSE`\. With `"auto"`, the patterns are returned only when compression pays \(few patterns relative to rows\): in memory, following the option `xtx.dedup` \(see #link("design_algebra.typ")[`design_algebra`]\) and at most half as many patterns as rows; remotely, following the option `xtx.remote_compress` \(see #link("design_remote.typ")[`design_remote`]\)\. The two backends can therefore decide differently with `"auto"`\. With `TRUE`, the patterns are returned unless every row is a pattern of its own\. With `FALSE`, `NULL` is returned\.],
  [`fit`], [`"lm"` \(default\) or `"glm"`: the fit the patterns are for, which sets the rule of `compress = "auto"` for a remote design \(a GLM iterates on the patterns, so compression pays sooner\)\. In memory both give the same result\.],
  [`check`], [Check the inputs \(default\)\. `FALSE` is for callers that pass columns of `x` they have already checked, as `xtx.fit::lm_xtx()` and `xtx.fit::glm_xtx()` do: `z`, `w` and `keys` are then a numeric vector \(or matrix, for `keys`\) with one value per row in memory, #link("design_column.typ")[`design_column()`]s of `x` remotely, and they are used as given \(remotely, without the query that checks them\)\. Names may change before the first release\.],
)

== Value

`NULL` when compression does not apply: with `compress = FALSE`; when every row is its own pattern; with `"auto"`, when compression does not pay; with fewer than two rows; in memory, when the per-pattern sums overflow; remotely, for a design with `na.action = na.pass`\. Otherwise a list with elements, one row or value per pattern \(patterns whose weights sum to zero are dropped\) in no particular order:

#terms(
  terms.item([`x`], [The design of the patterns: an in-memory `design_matrix` with the columns of `x`, one row per pattern \(`as.matrix(x)` gives the matrix\)\.]),
  terms.item([`z`], [The weighted mean of `z` in each pattern\.]),
  terms.item([`w`], [The sum of the weights in each pattern \(the number of rows, without `w`\)\. Always a numeric vector\.]),
  terms.item([`within`], [One number: the weighted sum of squares of `z` around the pattern means, summed over patterns\.]),
  terms.item([`rows`], [The number of rows of `x` summarized\.]),
  terms.item([`keys`], [A numeric matrix with one row per pattern and one column per key column \(named after the key columns\), or `NULL` without `keys`\.]),
)

== Details

In memory, the patterns are found from the factor codes and numeric variables of the design, without forming the dense matrix\. For a remote design \(see #link("design_remote.typ")[`design_remote`]\), the database groups the rows by the factors and numeric variables of the design in one scan, and only the patterns are returned to R\. Both give the same patterns and the same statistics, up to rounding and the order of the patterns\.

== Examples

```r
d <- data.frame(y = rnorm(100), f = gl(4, 25), x = rep(1:5, 20))
D <- design_matrix(y ~ f + x, d)
p <- design_patterns(D, d$y, compress = TRUE)
nrow(p$x)
fit <- lm.wfit(as.matrix(p$x), p$z, p$w)
all.equal(unname(coef(fit)), unname(coef(lm(y ~ f + x, d))))
# residual sum of squares of the fit on every row
sum(p$w * fit$residuals^2) + p$within

if (requireNamespace("duckdb", quietly = TRUE)) {
    con <- DBI::dbConnect(duckdb::duckdb())
    DBI::dbWriteTable(con, "d", d)
    R <- design_matrix(y ~ f + x, con, table = "d")
    pr <- design_patterns(R, "y", compress = TRUE)
    nrow(pr$x)
    DBI::dbDisconnect(con)
}
```
