#import "/.calepin/calepin.typ" as calepin
#show: calepin.document

#set document(title: "Select rows of a design matrix")
#metadata((
  title: "[.xtx_matrix",
  pagetitle: "Select rows of a design matrix — [.xtx_matrix",
  aliases: ("[.xtx_matrix",),
  "source-files": ("R/subset.R",),
)) <website-metadata>

#calepin.setup(echo: true, eval: true, results: "verbatim", warning: false, message: false)

#calepin.chunk(echo: false, results: "hide")[
```r
library(xtx)
```
]

= Select rows of a design matrix

== Description

`x[i, ]` keeps the rows `i` of a design matrix, without forming the dense matrix: the factor codes and numeric variables of the selected rows are copied, and every product and summary works on the result as on any other design\. For a design built by applying a plan to `newdata`, `design_matrix(plan, newdata)[i, ]` is identical to `design_matrix(plan, newdata[i, ])`, row names included, in the strict sense described in #link("design_matrix.typ")[`design_matrix()`]\. Scaled rows and combinations of two designs keep their scales for the selected rows\.

== Usage

```r
#| eval: false
# S3 method for class 'xtx_matrix'
x[i, j, ..., drop = FALSE]
```

== Arguments

#table(
  columns: 2,
  table.header([Argument], [Description]),
  [`x`], [A `design_matrix`\.],
  [`i`], [Rows to keep\.],
  [`j`], [Not supported: columns cannot be selected\.],
  [`...`], [Not used\.],
  [`drop`], [Ignored: the result is always a design matrix\.],
)

== Value

A `design_matrix` with `length(i)` rows \(or the number selected\)\.

== Details

Only rows can be selected\. `i` can be positive integers \(repeats allowed\), negative integers or a logical vector, as for a matrix; missing values are not allowed\.

== Examples

```r
d <- data.frame(x = rnorm(6), f = gl(3, 2))
D <- design_matrix(~ x * f, d)
as.matrix(D[c(2, 5), ])
```
