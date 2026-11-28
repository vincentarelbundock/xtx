#import "/.calepin/calepin.typ" as calepin
#show: calepin.document

#set document(title: "Sparse model matrices")
#metadata((
  title: "coerce,xtx_matrix,CsparseMatrix-method",
  pagetitle: "Sparse model matrices — sparse",
  aliases: ("coerce,xtx_matrix,CsparseMatrix-method", "coerce,xtx_matrix,dgCMatrix-method", "sparse"),
  "source-files": ("R/sparse.R",),
)) <website-metadata>

#calepin.setup(echo: true, eval: true, results: "verbatim", warning: false, message: false)

#calepin.chunk(echo: false, results: "hide")[
```r
library(xtx)
```
]

= Sparse model matrices

== Description

`as(x, "CsparseMatrix")` \(or `as(x, "dgCMatrix")`\) builds a sparse column-compressed matrix \(a `Matrix::dgCMatrix`\) from a design matrix, without forming the dense matrix\. The non-zero values are those of #link("https://rdrr.io/r/base/as.matrix.html")[`as.matrix()`], bit for bit; exact zeros are dropped and `NA` and `NaN` are kept\. The result is identical to converting `as.matrix(x)` to a general `dgCMatrix` \(through `"dMatrix"`, `"generalMatrix"` and `"CsparseMatrix"`\); unlike a plain `as(<matrix>, "CsparseMatrix")`, it never returns a triangular or symmetric class for square matrices\. Dropping zeros loses the sign of `-0`, which the dense matrix keeps\.

== Arguments

#table(
  columns: 2,
  table.header([Argument], [Description]),
  [`from`], [A `design_matrix`\.],
)

== Value

A `dgCMatrix` with the row and column names of `as.matrix(from)`\.

== Details

This is not the same matrix as #link("https://rdrr.io/pkg/Matrix/man/sparse.model.matrix.html")[`Matrix::sparse.model.matrix()`], which has its own rules for contrasts and column names; it follows #link("https://rdrr.io/pkg/stats/man/model.matrix.html")[`stats::model.matrix()`]\.

== Examples

```r
d <- data.frame(x = rnorm(20), f = gl(4, 5), g = gl(2, 10))
D <- design_matrix(~ x + f * g, d)
as(D, "CsparseMatrix")
```
