#import "/.calepin/calepin.typ" as calepin
#show: calepin.document

#set document(title: "Build a design without forming the model matrix")
#metadata((
  title: "design_matrix",
  pagetitle: "Build a design without forming the model matrix — design_matrix",
  aliases: ("design_matrix",),
  "source-files": ("R/design_matrix.R",),
)) <website-metadata>

#calepin.setup(echo: true, eval: true, results: "verbatim", warning: false, message: false)

#calepin.chunk(echo: false, results: "hide")[
```r
library(xtx)
```
]

= Build a design without forming the model matrix

== Description

`design_matrix()` takes the same arguments as #link("https://rdrr.io/pkg/stats/man/model.matrix.html")[`stats::model.matrix()`] and returns a `design_matrix` object: the prepared variables plus everything needed to build the matrix \(column names, contrasts, term layout\)\. The dense matrix is built by #link("https://rdrr.io/r/base/as.matrix.html")[`as.matrix()`], and algebra such as `design %*% beta` runs directly on the variables\.

== Usage

```r
#| eval: false
design_matrix(
  object,
  data = environment(object),
  contrasts.arg = NULL,
  xlev = NULL,
  na.action = stats::na.pass,
  ...,
  newdata = NULL
)
```

== Arguments

#table(
  columns: 2,
  table.header([Argument], [Description]),
  [`object`], [A model formula or `terms` object, a reusable design, or a fitted model\.],
  [`data`], [A data frame, as in #link("https://rdrr.io/pkg/stats/man/model.matrix.html")[`stats::model.matrix()`]\. For a design or plan, the new data to apply it to\. A DuckDB connection \(with `table` or `sql` in `...`\) or Parquet paths give a remote design: see #link("design_remote.typ")[`design_remote`]\.],
  [`contrasts.arg, xlev`], [As in #link("https://rdrr.io/pkg/stats/man/model.matrix.html")[`stats::model.matrix()`]\.],
  [`na.action`], [Only used when reusing a design or plan: how the model frame treats missing values\. Defaults to #link("https://rdrr.io/pkg/stats/man/na.pass.html")[`stats::na.pass()`], so the output has one row per row of `data`\.],
  [`...`], [With a fitted model and `newdata`, arguments selecting the model component, such as `component = "precision"` for a beta regression\. Without `newdata`, passed to #link("https://rdrr.io/pkg/stats/man/model.matrix.html")[`stats::model.matrix()`] for non-formula objects; otherwise ignored, as in `model.matrix.default()`\. For remote data, `table`, `sql`, `engine` and `chunk_size` \(see #link("design_remote.typ")[`design_remote`]\)\.],
  [`newdata`], [A new batch of observations\. For a fitted model, `design_matrix(fit, newdata = batch)` uses its training encoding\. Use `newdata = NULL` to prepare a specification for later batches\. A formula specification learns its encoding from training data when first applied\. Reuse that resulting design for subsequent batches\. Omitting both data arguments for a fitted model returns its original matrix\.],
)

== Value

`design_matrix()` returns a design matrix object or, with explicit `newdata = NULL`, a reusable specification\.

== Details

For a formula, `as.matrix(design_matrix(f, d))` is identical to `stats::model.matrix(f, d)` in the strict sense of `identical(a, b, num.eq = FALSE, single.NA = FALSE, attrib.as.set = FALSE)`: same values bit for bit \(including signed zeros and the payloads of `NA` and `NaN`\), same names and attributes in the same order\. The C routine is a port of the one in R's stats package; its fill visits rows in chunks, optionally on several threads \(see the option `xtx.threads` in #link("xtx-package.typ")[`xtx-package`]\), and performs exactly the same floating point operations as R for every cell\.

#strong[Options\.]

- `xtx.skip_na_copy` \(default `TRUE`\): build model frames with #link("design_frame.typ")[`design_frame()`], which returns the model frame of #link("https://rdrr.io/pkg/stats/man/model.frame.html")[`stats::model.frame()`] without the copy of the whole frame that #link("https://rdrr.io/pkg/stats/man/na.omit.html")[`stats::na.omit()`] and #link("https://rdrr.io/pkg/stats/man/na.exclude.html")[`stats::na.exclude()`] make when they remove nothing, and without recoding factors that already have the levels asked for in `xlev`\. The result is identical either way; set `FALSE` to build the frame with `model.frame()` itself, for example in benchmarks against functions that call it\.
- `xtx.structured` \(default `TRUE`\): use the column-structured kernels for products with finite designs \(see #link("design_algebra.typ")[`design_algebra`]\)\.
- `xtx.dedup`: see #link("design_algebra.typ")[`design_algebra`]\.

== Reusing a design or a fitted model on new data

`design_matrix(fit, newdata = batch)` builds the model matrix of new data the way #link("https://rdrr.io/pkg/stats/man/predict.lm.html")[`stats::predict.lm()`] does, from the pieces the model stores \(see #link("design_parts.typ")[`design_parts()`]\): the terms, with the parameters of transformations such as #link("https://rdrr.io/pkg/stats/man/poly.html")[`stats::poly()`] or #link("https://rdrr.io/pkg/splines/man/ns.html")[`splines::ns()`] fixed by the training data, the factor levels and the contrasts\. `design_matrix(D, newdata = batch)` does the same from a design `D` built on training data\. Its dense matrix is identical to

```
model.matrix(tt, model.frame(tt, batch, xlev = xlevels, na.action = na.action),
             contrasts.arg = contrasts)
```

with `tt` the training terms without the response\.

`design_matrix(fit, newdata = NULL)` returns that compiled specification \(a plan\) without applying it, for repeated use: `design_matrix(plan, batch)`\. When a batch is a data frame whose variables have the same types, levels and column names as the first one, the plan skips #link("https://rdrr.io/pkg/stats/man/model.frame.html")[`stats::model.frame()`] and the contrast set-up, which makes small calls much faster; anything else takes the general path, so results never depend on the path taken\. Plans of fitted models are cached: two models with identical terms, levels and contrasts share one\.

== Examples

```r
d <- data.frame(y = rnorm(10), x = rnorm(10), f = gl(2, 5))
X <- as.matrix(design_matrix(y ~ x * f, d))
identical(X, model.matrix(y ~ x * f, d))

D <- design_matrix(y ~ x * f, d)
D %*% c(1, 2, 3, 4)
D_new <- design_matrix(D, d[1:3, ])
as.matrix(D_new)
```
