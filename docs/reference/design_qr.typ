#import "/.calepin/calepin.typ" as calepin
#show: calepin.document

#set document(title: "QR decomposition of a design matrix without forming the dense matrix")
#metadata((
  title: "design_qr",
  pagetitle: "QR decomposition of a design matrix without forming the dense matrix — design_qr",
  aliases: ("design_qr", "qr.Q,xtx_qr-method", "qr.R,xtx_qr-method", "qr.coef,xtx_qr,ANY-method", "qr.fitted,xtx_qr,ANY-method", "qr.qty,xtx_qr,ANY-method", "qr.qy,xtx_qr,ANY-method", "qr.resid,xtx_qr,ANY-method", "qr.xtx_matrix", "solve.xtx_qr"),
  "source-files": ("R/qr.R",),
)) <website-metadata>

#calepin.setup(echo: true, eval: true, results: "verbatim", warning: false, message: false)

#calepin.chunk(echo: false, results: "hide")[
```r
library(xtx)
```
]

= QR decomposition of a design matrix without forming the dense matrix

== Description

`qr(x)` decomposes a `design_matrix` as `X P = Q R` without forming the dense matrix `X` or the `n` by `p` matrix `Q`\. It returns an object of class `"xtx_qr"`, which the usual base R functions accept: `qr.coef()`, `qr.fitted()`, `qr.resid()`, `qr.qty()`, `qr.qy()`, `qr.Q()`, `qr.R()` and `solve()`\. This gives least squares coefficients, fitted values and residuals at the cost of a few passes over the rows\.

== Usage

```r
#| eval: false
# S3 method for class 'xtx_matrix'
qr(x, tol = 1e-07, ...)

# S4 method for signature 'xtx_qr,ANY'
qr.coef(qr, y)

# S4 method for signature 'xtx_qr,ANY'
qr.fitted(qr, y, k = qr$rank)

# S4 method for signature 'xtx_qr,ANY'
qr.resid(qr, y)

# S4 method for signature 'xtx_qr,ANY'
qr.qty(qr, y)

# S4 method for signature 'xtx_qr,ANY'
qr.qy(qr, y)

# S4 method for signature 'xtx_qr'
qr.Q(qr, complete = FALSE, Dvec)

# S4 method for signature 'xtx_qr'
qr.R(qr, complete = FALSE)

# S3 method for class 'xtx_qr'
solve(a, b, ...)
```

== Arguments

#table(
  columns: 2,
  table.header([Argument], [Description]),
  [`x`], [A `design_matrix`\.],
  [`tol`], [The tolerance for detecting aliased columns, as in #link("https://rdrr.io/r/base/qr.html")[`qr()`]\.],
  [`...`], [Not used, except `y` for a remote design: a column of the design \(a #link("design_column.typ")[`design_column()`]\) whose products `X'y` are computed in the scan that computes `X'X`, and kept for `qr.coef(qr, y)` on that column\.],
  [`qr`], [A `"xtx_qr"` object, from `qr()` on a `design_matrix`\.],
  [`y`], [A numeric vector or matrix: the response for `qr.coef()`, `qr.fitted()`, `qr.resid()` and `qr.qty()` \(`nrow(x)` rows\), or the values to multiply by `Q` for `qr.qy()` \(`rank` rows\)\.],
  [`k`], [The rank for `qr.fitted()`; only `qr$rank` is supported\.],
  [`complete, Dvec`], [Not supported; only the defaults are accepted\.],
  [`a, b`], [For `solve()`: a `"xtx_qr"` object and the right-hand side, as in #link("https://rdrr.io/r/base/solve.qr.html")[`solve.qr()`]\.],
)

== Value

An object of class `"xtx_qr"`: a list with `design` \(`x`\), `R` \(the `p` by `p` upper triangular factor, columns in pivoted order\), `rank`, `pivot` and `tol`\. For a remote design, also the sums that its scan computed: `weight_sum` \(the sum of the weights of the rows: their number for an unscaled design\) and, with `y`, `response_ss` \(the sum of squares of `y`\)\. `qr.coef()` on a remote design returns the coefficients with an attribute `"moments"`: a list with `b` \(the coefficients of the kept columns\), `h` \(their `X'y`\) and `G` \(their `X'X`\), from which sums of squares follow without another pass over the rows \(`y'y - sum(b * h)` is the residual sum of squares\); these are those of the centered design when the factor is\.

== Details

#strong[Method\.] `R` is the Cholesky factor of `X'X`, computed by `crossprod(x)` without forming `X`\. Columns are pivoted as in `qr(as.matrix(x))` \(LINPACK's limited pivoting, used by #link("https://rdrr.io/r/base/lm.html")[`lm()`]\): they are taken in order, and a column whose residual norm is below `tol` times its norm is moved to the end and counts as aliased\. From `X'X`, a residual norm is accurate only to about `sqrt(.Machine$double.eps)` times the column norm, which cannot separate an aliased column from one just above `tol`; the residual norms of such columns are computed from the data \(one product with `x` per block of columns\)\. `Q` is `X[, pivot] %*% solve(R)`, never formed\. Coefficients solve the normal equations and are then refined from the residuals `y - X b`, each step costing one `x %*% b` and one `crossprod(x, r)`\. Refinement recovers much of the accuracy that the normal equations lose on badly conditioned designs\. When `X'X` is badly conditioned and numeric variables are far from zero \(`x + 1e7`\), the values of `R` come from the `X'X` of the design with those variables centered, when that is an exact linear transform `X T` of the design \(`R` is then `R~ T^-1`\); the columns kept and their order still come from `X'X`\.

#strong[Agreement with `qr()`\.] Results agree with those of `qr(as.matrix(x))` up to floating point rounding, not bit for bit\. On badly conditioned designs the coefficients agree to fewer digits than with a Householder QR\. A column whose residual norm is within rounding of `tol` can be aliased by one decomposition and not the other\. Rows of `R` can differ in sign \(the diagonal here is positive\), and so can the columns of `Q` and the values of `qr.qty()`\. For an aliased design, the block of `R` for the aliased columns is zero\.

#strong[What is not formed\.] The complete `Q` \(`n` by `n`\) is never formed, so:

- `qr.Q(q)` returns the `n` by `rank` matrix of the first `rank` columns of `Q`, which spans the columns of `X`; `complete = TRUE` is not supported\.
- `qr.qty(q, y)` returns the first `rank` rows of `t(Q) %*% y`\.
- `qr.qy(q, y)` takes `y` with `rank` rows and returns `Q[, 1:rank] %*% y`\.
- `qr.fitted()` only supports `k = q$rank`\.

Non-finite entries in `x` are an error, as in #link("https://rdrr.io/r/base/qr.html")[`qr()`]\. `y` must be finite\. When `X'X` cannot represent the columns' scale \(entries whose squares overflow or underflow, around `1e154` and `1e-154`\), `qr()` falls back to the dense QR of `qr(as.matrix(x))`, kept in the object, and the methods give its results\.

#strong[Calling the methods\.] `qr.coef()`, `qr.fitted()`, `qr.resid()`, `qr.qty()`, `qr.qy()`, `qr.Q()` and `qr.R()` are S4 generics created by this package \(as the Matrix package does for its decompositions\): they dispatch when xtx is attached, or when called as `xtx::qr.coef()`\. Calls to `base::qr.coef()` \(for example from code in another package that imports base's function\) do not dispatch, and stop because the object is not a base QR decomposition\. `qr()` and `solve()` are S3 generics and always dispatch\.

== Examples

```r
set.seed(1)
d <- data.frame(y = rnorm(100), x = rnorm(100), f = gl(4, 25))
D <- design_matrix(y ~ x * f, d)
q <- qr(D)
all.equal(qr.coef(q, d$y), coef(lm(y ~ x * f, d)))
all.equal(qr.fitted(q, d$y), unname(fitted(lm(y ~ x * f, d))))
```
