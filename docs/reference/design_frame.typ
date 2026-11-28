#import "/.calepin/calepin.typ" as calepin
#show: calepin.document

#set document(title: "Build a model frame without redundant copies, or the remote frame of a model")
#metadata((
  title: "design_frame",
  pagetitle: "Build a model frame without redundant copies, or the remote frame of a model — design_frame",
  aliases: ("design_frame",),
  "source-files": ("R/frame.R",),
)) <website-metadata>

#calepin.setup(echo: true, eval: true, results: "verbatim", warning: false, message: false)

#calepin.chunk(echo: false, results: "hide")[
```r
library(xtx)
```
]

= Build a model frame without redundant copies, or the remote frame of a model

== Description

`design_frame()` takes the arguments of #link("https://rdrr.io/pkg/stats/man/model.frame.html")[`stats::model.frame()`] and returns the same model frame, faster on large data\. #link("design_matrix.typ")[`design_matrix()`] uses it to build its model frames\. When `data` is remote data \(a DuckDB connection, or Parquet paths\), it returns the #emph[remote frame] of the model instead: the remote design of the rows the model uses, with its response, weights and offset as columns over those rows\.

== Usage

```r
#| eval: false
design_frame(...)
```

== Arguments

#table(
  columns: 2,
  table.header([Argument], [Description]),
  [`...`], [Arguments of #link("https://rdrr.io/pkg/stats/man/model.frame.html")[`stats::model.frame()`]: `formula`, `data`, `subset`, `na.action`, `drop.unused.levels`, `xlev`, and further variables such as `weights` or `offset`\. With remote data, see #emph[Remote data]\.],
)

== Value

In memory, a model frame, as returned by #link("https://rdrr.io/pkg/stats/man/model.frame.html")[`stats::model.frame()`]\.

With remote data, a list of class `xtx_remote_frame`:

- `x`: the remote design of the rows used\.
- `response`: a list with `kind` and `ncol`, and by kind: `"numeric"` \(a numeric column or an expression\) and `"logical"` \(a logical column, as 0/1\): `y`, a #link("design_column.typ")[`design_column()`] \(`ncol` is more than 1 when an expression evaluated in R gave a matrix; `y` is then its first column\); `"factor"` \(a character or factor column\): `y`, its values as a text column, which compares with a string \(`y == "a"`\), and `levels`, the levels on the rows used; `"matrix"` \(`cbind()` of `ncol` expressions\): `cols`, a list of its first two columns \(`cbind()` of one expression is that expression, but for a factor column, whose codes are then the one column of `cols`, as in memory\); `"other"` \(a column of another type, such as a date\) and `"none"` \(a one-sided formula\): nothing\.
- `weights`: a column, or `NULL`\.
- `offset`: the sum of the `offset()` terms of the formula and of the `offset` argument, a column, or `NULL`\.
- `ncol`: the numbers of columns of the `weights` and of the `offset`, more than 1 when an expression evaluated in R gave a matrix \(which a caller refuses, as #link("https://rdrr.io/pkg/stats/man/lm.html")[`stats::lm()`] does\)\.
- `types`: the kind of each column of the data \(`"double"`, `"integer"`, `"logical"`, `"character"`, `"factor"` or `"other"`\)\.
- `terms`, `contrasts`, `xlevels`: the terms of the formula, and the contrasts and factor levels of the design\.
- `na.action`: `"omit"`, `"exclude"`, `"fail"` or `"pass"`\.

== Details

#strong[In memory\.] The result is identical to that of `model.frame()` in the strict sense described in #link("design_matrix.typ")[`design_matrix()`]: same values, same attributes in the same order, same row names; warnings and errors are the same too\. `design_frame()` skips two steps of `model.frame()` when they would return their input unchanged:

- #link("https://rdrr.io/pkg/stats/man/na.omit.html")[`stats::na.omit()`] and #link("https://rdrr.io/pkg/stats/man/na.exclude.html")[`stats::na.exclude()`] copy every column of the frame and check its row names for duplicates, even when no row has a missing value\. `design_frame()` skips the copy when the frame has nothing to omit, and makes the check with a faster one\.
- With `xlev`, `model.frame()` recodes every factor twice, even when it already has the requested levels\. `design_frame()` leaves such factors as they are\.

Anything else runs the code of `model.frame()`: other `na.action` functions, columns with a class or with attributes that subsetting would drop or reorder, factors whose levels differ from `xlev`, and objects with their own `model.frame()` method \(the method is chosen by #link("https://rdrr.io/pkg/stats/man/model.frame.html")[`stats::model.frame()`] itself\)\.

`design_frame()` is built when the package is loaded, from the `model.frame.default()` of the running version of R\. If that function, the data frame methods of `na.omit()` and `na.exclude()` or `[.data.frame` are not the ones this package was written for, or if a self-test at load finds a difference, the corresponding shortcut is turned off; without shortcuts, `design_frame()` is `model.frame()`\. Option `xtx.skip_na_copy = FALSE` turns both shortcuts off\.

#strong[Remote data\.] `data` is looked at first \(it is evaluated before the other arguments, after `formula`\)\. When it is a DuckDB connection or a character vector of Parquet paths, the arguments are, all by name except `formula` and `data`:

- `formula`: a model formula\. The design is that of its right-hand side \(see #link("design_remote.typ")[`design_remote`]\)\.
- `weights`, `offset`, `subset`: taken unevaluated, as `model.frame()` takes them: the name of a column, a `DBI::SQL()` expression, or an R expression over the columns\. Expressions with a SQL translation run in the database; the others are evaluated in R, batch by batch, into a temporary table \(with a message; `engine = "sql"` refuses\)\. R expressions are evaluated in the environment of the formula\.
- `na.action`: #link("https://rdrr.io/pkg/stats/man/na.omit.html")[`stats::na.omit()`] \(the default\) or #link("https://rdrr.io/pkg/stats/man/na.exclude.html")[`stats::na.exclude()`] \(the same rows\), #link("https://rdrr.io/pkg/stats/man/na.fail.html")[`stats::na.fail()`] or #link("https://rdrr.io/pkg/stats/man/na.pass.html")[`stats::na.pass()`] \(the errors of #link("https://rdrr.io/pkg/stats/man/lm.html")[`stats::lm()`] when a row has a missing value\), or their names\.
- `contrasts`: as the `contrasts.arg` of #link("https://rdrr.io/pkg/stats/man/model.matrix.html")[`stats::model.matrix()`]\.
- `table`, `sql`, `engine`, `chunk_size`: as #link("design_matrix.typ")[`design_matrix()`] takes them for remote data \(see #link("design_remote.typ")[`design_remote`]\)\.
- other arguments: bare names of columns that are not in the model but must not be missing on the rows used \(the variables of fixed effects or clusters, say\), as `model.frame()` keeps the extra variables it is given\.

The rows used are those of the subset where no variable of the model, response, weight, offset or extra column is missing, and the levels of factors are those of these rows, as `model.frame(drop.unused.levels = TRUE)` gives them in memory\. The response is returned whatever its type, as `model.frame()` does; a column given as `weights` or `offset` must be numeric \(the errors of `lm()`: `'weights' must be a numeric vector`, `'offset' must be numeric`\)\.

== Examples

```r
d <- data.frame(y = rnorm(10), x = rnorm(10), f = gl(2, 5))
identical(design_frame(y ~ x * f, d), model.frame(y ~ x * f, d))

if (requireNamespace("duckdb", quietly = TRUE)) {
    con <- DBI::dbConnect(duckdb::duckdb())
    d$w <- runif(10)
    d$id <- rep(1:5, 2)
    DBI::dbWriteTable(con, "d", d)
    fr <- design_frame(y ~ x * f, con, table = "d", weights = 2 * w, subset = x > -1, id)
    fr
    dim(fr$x)
    sum(fr$weights)
    b <- solve(crossprod(fr$x), crossprod(fr$x, fr$response$y))
    DBI::dbDisconnect(con)
}
```
