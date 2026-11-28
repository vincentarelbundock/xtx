#import "/.calepin/calepin.typ" as calepin
#show: calepin.document

#set document(title: "Algebra on design matrices without forming the dense matrix")
#metadata((
  title: "%*%.xtx_matrix",
  pagetitle: "Algebra on design matrices without forming the dense matrix — design_algebra",
  aliases: ("%*%.xtx_matrix", "colMeans,xtx_matrix-method", "colSums,xtx_matrix-method", "crossprod.xtx_matrix", "design_algebra", "quadratic_diagonal", "rowMeans,xtx_matrix-method", "rowSums,xtx_matrix-method", "rowsum.xtx_matrix"),
  "source-files": ("R/algebra.R",),
)) <website-metadata>

#calepin.setup(echo: true, eval: true, results: "verbatim", warning: false, message: false)

#calepin.chunk(echo: false, results: "hide")[
```r
library(xtx)
```
]

= Algebra on design matrices without forming the dense matrix

== Description

Products with a `design_matrix` are computed directly from the factor codes and numeric variables, with the usual base R operators:

== Usage

```r
#| eval: false
# S3 method for class 'xtx_matrix'
rowsum(x, group, reorder = TRUE, na.rm = FALSE, ...)

# S4 method for signature 'xtx_matrix'
colSums(x, na.rm = FALSE, dims = 1, ..., group = NULL, reorder = TRUE)

# S4 method for signature 'xtx_matrix'
colMeans(x, na.rm = FALSE, dims = 1, ..., group = NULL, reorder = TRUE)

# S4 method for signature 'xtx_matrix'
rowSums(x, na.rm = FALSE, dims = 1, ...)

# S4 method for signature 'xtx_matrix'
rowMeans(x, na.rm = FALSE, dims = 1, ...)

quadratic_diagonal(x, V)

# S3 method for class 'xtx_matrix'
x %*% y

# S3 method for class 'xtx_matrix'
crossprod(x, y = NULL, ...)
```

== Arguments

#table(
  columns: 2,
  table.header([Argument], [Description]),
  [`x`], [A `design_matrix`\.],
  [`group`], [A vector or factor with one value per row, as in #link("https://rdrr.io/pkg/base/man/rowsum.html")[`base::rowsum()`]\.],
  [`reorder`], [As in #link("https://rdrr.io/pkg/base/man/rowsum.html")[`base::rowsum()`]: if `TRUE`, the groups are sorted; otherwise they are in order of first appearance\.],
  [`na.rm, dims`], [Only the defaults are supported\.],
  [`...`], [Not used\.],
  [`V`], [A square numeric matrix with `ncol(x)` rows\.],
  [`y`], [For `%*%`, a numeric vector of length `ncol(x)` or a matrix with `ncol(x)` rows \(the coefficients\)\. For `crossprod()`, a numeric vector with one value per row or a matrix with `nrow(x)` rows \(the weights\), another `design_matrix` with the same rows, or `NULL` for `t(X) %*% X`\.],
)

== Value

`rowsum()`: a matrix with one row per group and the columns of `x`\. `colSums()` and `colMeans()`: a vector with one value per column, or with `group`, a matrix with one row per observed group \(sorted by default\)\. `rowSums()`, `rowMeans()` and `quadratic_diagonal()`: a vector with one value per row \(for a remote design, a column: see #link("design_column.typ")[`design_column`]\)\.

== Details

- `x %*% beta`: `X %*% beta`, an `nrow(x)` by `ncol(beta)` matrix\.
- `crossprod(x, w)`: `t(X) %*% w`, an `ncol(x)` by `ncol(w)` matrix\.
- `crossprod(x)` and `crossprod(x, y)` with `y` another `design_matrix` with the same rows: `t(X) %*% X` and `t(X) %*% Y`, ordinary base R matrices\. Weighted crossproducts are `crossprod(x, w * x)` or `crossprod(sqrt(w) * x)`\. The result can go straight to #link("https://rdrr.io/r/base/chol.html")[`chol()`], #link("https://rdrr.io/r/base/chol2inv.html")[`chol2inv()`], #link("https://rdrr.io/r/base/eigen.html")[`eigen()`] or #link("https://rdrr.io/r/base/solve.html")[`solve()`]; see also #link("design_qr.typ")[`qr()`] for least squares\.
- `colSums(x, group = g)` and `colMeans(x, group = g)` summarize each design column within groups of rows\. Group order, row names and missing group handling follow #link("https://rdrr.io/pkg/base/man/rowsum.html")[`base::rowsum()`]\. For weighted sums, use `colSums(w * x, group = g)`\. For design matrices, `rowsum(x, g)` is deprecated in favor of `colSums(x, group = g)`\.
- `colSums(x)`, `colMeans(x)`, `rowSums(x)` and `rowMeans(x)`: as for the dense matrix \(without `na.rm` or `dims`\)\.
- `quadratic_diagonal(x, V)`: `rowSums((X %*% V) * X)`, the diagonal of `X V X'`\. For a Jacobian and a coefficient covariance matrix, these are the variances of unit-level estimates\.

The results agree with the same computation on the dense matrix up to floating point rounding \(the order of the sums differs\), with `NA` and `NaN` in the same places\. Only `design %*% beta`, `crossprod(design, w)` and `crossprod(design, design)` are supported, with a design first\.

`beta` and `V` must be finite: the kernels skip structural zeros, so `0 * Inf` terms that the dense product would turn into `NaN` are not computed\.

Designs that the kernels cannot handle \(non-formula objects, terms with more than 32 variables\) fall back to dense algebra\.

#strong[Identical rows\.] Many data sets have far fewer distinct covariate patterns than rows\. The kernels can group rows whose inputs \(factor codes, numeric values, scales and group\) have the same bits, compute once per pattern, and expand\. Row-wise results \(`%*%`, `quadratic_diagonal()`\) are then identical bit for bit to computing every row; `crossprod()` and grouped column summaries first sum the weights within each pattern, which changes only the rounding\. The groups are found once per object and cached\. The option `xtx.dedup` controls this: `"auto"` \(default\) groups when the work per row is large enough to pay for hashing it and at most half the rows are distinct, `TRUE` always groups and `FALSE` never does\.

== Examples

```r
d <- data.frame(x = rnorm(100), f = gl(4, 25), g = gl(2, 1, 100))
D <- design_matrix(~ x * f, d)
X <- model.matrix(~ x * f, d)
b <- rnorm(ncol(D))
all.equal(D %*% b, X %*% b, check.attributes = FALSE)
crossprod(D, rep(1, 100))
colMeans(D, group = d$g)
colSums(runif(100) * D, group = d$g)
colMeans(D)
quadratic_diagonal(D, diag(ncol(D)))
```
