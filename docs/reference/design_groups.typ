#import "/.calepin/calepin.typ" as calepin
#show: calepin.document

#set document(title: "Groupings of the rows of a design")
#metadata((
  title: "design_groups",
  pagetitle: "Groupings of the rows of a design — design_groups",
  aliases: ("design_groups", "format.xtx_groups", "group_crosstab", "group_lookup", "group_nested", "group_sums", "print.xtx_groups"),
  "source-files": ("R/groups_public.R",),
)) <website-metadata>

#calepin.setup(echo: true, eval: true, results: "verbatim", warning: false, message: false)

#calepin.chunk(echo: false, results: "hide")[
```r
library(xtx)
```
]

= Groupings of the rows of a design

== Description

A #emph[grouping] gives each row of a design the number of its group, from the values of one or more key variables: integer codes `1..G`, one per level of the key\. Grouped sums, crosstabs, lookups of a value per group and nesting tests then work from the codes, in memory and on a remote design \(see #link("design_remote.typ")[`design_remote`]\) with the same code and the same results\. These are the operations fixed-effect estimators are built from: for instance, the within transformation of a column is `y - group_lookup(group_sums(y, g) / n_g, g)`\.

== Usage

```r
#| eval: false
design_groups(x, key, levels = NULL)

group_sums(col, g, total = NULL)

group_crosstab(g1, g2, w = NULL)

group_lookup(values, g)

group_nested(g1, g2)

# S3 method for class 'xtx_groups'
format(x, ...)

# S3 method for class 'xtx_groups'
print(x, ...)
```

== Arguments

#table(
  columns: 2,
  table.header([Argument], [Description]),
  [`x`], [For `design_groups()`, a design: a `design_matrix` or a remote design \(a row scale, `s * D`, is ignored\)\. For `format()` and `print()`, a grouping\.],
  [`key`], [The key: one or more variables whose combinations of values are the groups\. A character vector is a vector of variable #emph[names]: of the remote data for a remote design \(columns of any type\); in memory, of the variables of the model frame of the design \(the response and the variables of the formula\), the only data a design in memory keeps\. In memory, a key can also be given by its values: a vector \(numeric, logical, factor\) with one value per row of `x`, or a data frame \(such as `d["firm"]`\) or a list of such vectors \(strings as values must be in a data frame or a list, since a character vector is a vector of names\)\. For a remote design, also a `DBI::SQL()` expression over the columns of the remote data, or a #link("design_column.typ")[`design_column()`] of `x` \(each of its columns a key variable\)\. A list combines several of these\.],
  [`levels`], [`NULL`, or the levels of a grouping made before \(its `levels`, a data frame with one column per key variable\): the rows of `x` are then coded against those levels, as new data against the levels of a fit, instead of the levels that occur in `x`\. Rows whose key is missing or not among `levels` get the code `NA` \(no error\): such a grouping works with `group_lookup()` \(`NA` for those rows\) and the predictions of `xtx.fit::lm_xtx()` with fixed effects, not with `group_sums()`, `group_crosstab()` or `group_nested()`\. Names may change before the first release\.],
  [`col`], [A column: in memory a numeric \(or logical\) vector with one value per row, or a matrix or data frame of such columns; for a remote design a #link("design_column.typ")[`design_column()`] \(with one or more columns\)\. A character vector of names, or a `DBI::SQL()` expression, is taken as `design_column(g$design, .)` \(in memory, variables of the model frame\)\. A list of columns is bound into one with several columns\. A single number is a constant column: `group_sums(1, g)` counts the rows of each group\.],
  [`g, g1, g2`], [Groupings made by `design_groups()`, of the same design, or bare codes \(see #link("design_key.typ")[`key_codes()`]\), which are used as given, without the checks\. For `group_sums()`, `g` can also be a list of groupings: the sums within each of them, from one pass over the rows\.],
  [`total`], [For `group_sums()` with a list of groupings: `NULL`, or a column whose sum over all rows is computed in the same pass\.],
  [`w`], [Optional weights, a column in the same form as `col` \(one column\)\. Without `w`, the counts of rows\.],
  [`values`], [The value of each group: a numeric vector of length `G`, or a matrix \(or data frame\) with `G` rows for several columns, in code order\. Missing values give missing values\.],
  [`...`], [Not used\.],
)

== Value

`design_groups()`: a grouping \(see #emph[The grouping object]\)\. `group_sums()` with a list of groupings: a list with `groups` \(one matrix of sums per grouping\) and `total`\. `group_sums()`: a `G x k` numeric matrix, `k` the number of columns of `col`\. `group_crosstab()`: a `G1 x G2` #link("https://rdrr.io/pkg/Matrix/man/dgCMatrix-class.html")[`Matrix::dgCMatrix`] of the sums of `w` \(counts without `w`\) in each pair of groups; pairs with no row are not stored\. `group_lookup()`: a column \(in memory a double vector, or a matrix for a matrix of `values`; remotely a #link("design_column.typ")[`design_column()`]\)\. `group_nested()`: `TRUE` or `FALSE`\.

== Details

- `design_groups(x, key)` makes the grouping of the rows of design `x` by `key`\.
- `group_sums(col, g)`: the sums of a column \(or several\) within each group, a `G x k` matrix \(one query remotely\)\.
- `group_crosstab(g1, g2, w)`: the \(weighted\) number of rows in each pair of groups, a sparse `G1 x G2` #link("https://rdrr.io/pkg/Matrix/man/dgCMatrix-class.html")[`Matrix::dgCMatrix`]\.
- `group_lookup(values, g)`: the column whose value in a row is the value of its group, `values[codes]`\.
- `group_nested(g1, g2)`: is `g1` nested in `g2` \(does each group of `g1` lie within one group of `g2`\)?

== The grouping object

A list of class `xtx_groups` with elements

#terms(
  terms.item([`codes`], [The group of each row, integers `1..G`: in memory an integer vector with one value per row of the design; for a remote design, a #link("design_column.typ")[`design_column()`] \(lazy: a join of the rows with a table of the levels, made when its values are first needed\)\. `codes` is a column like any other: `as.vector(g$codes)` reads it\.]),
  terms.item([`G`], [The number of groups \(an integer\)\.]),
  terms.item([`levels`], [A data frame with `G` rows, one column per key variable: the values of the key in each group, in code order \(row `l` is group `l`\)\.]),
  terms.item([`n`], [The number of rows of the design \(remotely, `NULL` when not yet known\)\.]),
  terms.item([`design`], [The design\.]),
)

Treat it as read-only: the functions check that the codes are within `1..G` before they use them\.

== Levels and their order

The groups are the combinations of key values that occur in the rows of the design, sorted as the levels of fixed effects are sorted \(as `xtx.fit::lm_xtx()` numbers them\): by the first key variable, then the second, and so on; numbers in increasing order, strings in the C locale \(as `sort(method = "radix")`: upper case before lower case\), `FALSE` before `TRUE`\. Logical keys give integer levels `0` and `1`\. Factors in memory are sorted in the order of their levels; a remote data source has no factors \(strings, and DuckDB `ENUM` columns, are sorted as strings\): for the same groups on both backends, use character keys, or factors whose levels are sorted\.

== Missing keys

Every row of the design must have a group: a key with a missing value \(`NA` or `NaN` in memory, SQL `NULL` or `NaN` remotely\) in a row of the design is an error, on both backends\. Fits with fixed effects drop such rows; a grouping has no rows of its own, so select the rows first: for a remote design, `D[DBI::SQL("k IS NOT NULL"), ]`; in memory, make the design from the rows with a key\. Missing values in the #emph[columns] are kept: as #link("https://rdrr.io/r/base/rowsum.html")[`rowsum()`] does, `group_sums()` gives `NA` for a group with a missing value, and `group_crosstab()` an `NA` count for a pair with a missing weight\.

== Rows, designs and backends

In memory, columns are plain vectors \(or matrices\) with one value per row of the design, in the order of its rows\. For a remote design, columns are #link("design_column.typ")[`design_column()`]s of that design \(or, where a column is expected, the name of a column of the remote data or a `DBI::SQL()` expression\): remote rows have no order, so a vector in R cannot be matched with them, and results with one value per row \(`group_lookup()`\) are columns\. A grouping combines only with columns and groupings of the design it was made from \(a row selection `D[i, ]` is another design; in memory, where columns are plain vectors, only their length is checked\); an in-memory grouping does not combine with remote columns, nor the converse\. Results with one row per group \(`group_sums()`, `group_crosstab()`\) are in code order, in R, without dimension names: row `l` is `g$levels[l, ]`\. Both backends give the same values, up to rounding of the sums\.

A remote grouping is valid while its connection is open \(see #link("design_column.typ")[`design_column()`], #emph[Lifetime]\): `design_groups()` reads the levels \(one `SELECT DISTINCT`\) and writes them to a temporary table, and each use joins that table to the rows, in the query that reads them\.

== Examples

```r
set.seed(1)
d <- data.frame(y = rnorm(200), x = rnorm(200),
                firm = sample(c("a", "b", "c", "d"), 200, TRUE))
d$region <- ifelse(d$firm %in% c("a", "b"), "north", "south")
D <- design_matrix(y ~ x, d)

# the same code for both backends: the within (one-way fixed effect)
# slope of y on x
within_slope <- function(D, g) {
    n_g <- group_sums(1, g)[, 1]      # rows per group
    demean <- function(v) v - group_lookup(group_sums(v, g)[, 1] / n_g, g)
    yt <- demean(design_column(D, "y"))
    xt <- demean(design_column(D, "x"))
    sum(xt * yt) / sum(xt * xt)
}

g <- design_groups(D, d["firm"])
g
group_sums(design_column(D, "y"), g)
within_slope(D, g)
coef(lm(y ~ x + firm, d))[["x"]]
group_nested(g, design_groups(D, d["region"]))

if (requireNamespace("duckdb", quietly = TRUE)) {
    con <- DBI::dbConnect(duckdb::duckdb())
    DBI::dbWriteTable(con, "d", d)
    R <- design_matrix(y ~ x, con, table = "d")
    gr <- design_groups(R, "firm")
    gr$levels
    group_sums("y", gr)
    within_slope(R, gr)
    group_nested(gr, design_groups(R, "region"))
    group_crosstab(gr, design_groups(R, c("region", "firm")))
    DBI::dbDisconnect(con)
}
```
