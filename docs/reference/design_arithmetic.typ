#import "/.calepin/calepin.typ" as calepin
#show: calepin.document

#set document(title: "Arithmetic on design matrices")
#metadata((
  title: "Ops.xtx_matrix",
  pagetitle: "Arithmetic on design matrices — design_arithmetic",
  aliases: ("Ops.xtx_matrix", "design_arithmetic"),
  "source-files": ("R/algebra.R",),
)) <website-metadata>

#calepin.setup(echo: true, eval: true, results: "verbatim", warning: false, message: false)

#calepin.chunk(echo: false, results: "hide")[
```r
library(xtx)
```
]

= Arithmetic on design matrices

== Description

Scaling rows and adding or subtracting design matrices gives another `design_matrix`, still without forming the dense matrix\. Row `i` of the result is `s1[i] * x1_i + s2[i] * x2_i`, where `x1_i` and `x2_i` are rows of at most two designs with the same rows and columns\.

== Usage

```r
#| eval: false
# S3 method for class 'xtx_matrix'
Ops(e1, e2)
```

== Arguments

#table(
  columns: 2,
  table.header([Argument], [Description]),
  [`e1, e2`], [A `design_matrix`, or a number or numeric vector for `*` and `/`\.],
)

== Value

A `design_matrix`\.

== Details

- `s * D` and `D * s`: scale the rows of `D` by `s`, a number or a numeric vector with one value per row\.
- `D / s`: the same with `1 / s`\.
- `-D`: the same with `-1`\.
- `D1 + D2` and `D1 - D2`: add or subtract two design matrices with the same rows and columns\.

This covers the Jacobians used for delta-method standard errors: the Jacobian of predictions from a model with a linear predictor is `dmu * D`, where `dmu` is the derivative of the inverse link, and the Jacobian of a comparison between two counterfactual data sets is `dmu1 * D1 - dmu0 * D0`\. It also covers weighted designs, such as `w * D` in `crossprod(D, w * D)`\.

Scales are combined before any product is taken: `as.matrix(s2 * (s1 * D))` is `(s2 * s1) * as.matrix(D)`\. A result that would combine more than two designs is an error\.

== Examples

```r
d <- data.frame(x = rnorm(50), f = gl(5, 10))
p <- design_matrix(~ x * f, d)
X0 <- design_matrix(p, d)
X1 <- design_matrix(p, transform(d, x = x + 1))
J <- X1 - X0
V <- diag(10) / 100
quadratic_diagonal(J, V)
```
