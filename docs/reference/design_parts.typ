#import "/.calepin/calepin.typ" as calepin
#show: calepin.document

#set document(title: "Terms, factor levels and contrasts of a fitted model")
#metadata((
  title: "design_parts",
  pagetitle: "Terms, factor levels and contrasts of a fitted model — design_parts",
  aliases: ("design_parts", "design_parts.afex_aov", "design_parts.betareg", "design_parts.clm", "design_parts.clmm", "design_parts.coxph", "design_parts.default", "design_parts.glm", "design_parts.glmmTMB", "design_parts.gls", "design_parts.hurdle", "design_parts.iv_robust", "design_parts.ivreg", "design_parts.lm", "design_parts.lm_robust", "design_parts.lme", "design_parts.merMod", "design_parts.multinom", "design_parts.polr", "design_parts.survreg", "design_parts.svyolr", "design_parts.zeroinfl"),
  "source-files": ("R/parts.R",),
)) <website-metadata>

#calepin.setup(echo: true, eval: true, results: "verbatim", warning: false, message: false)

#calepin.chunk(echo: false, results: "hide")[
```r
library(xtx)
```
]

= Terms, factor levels and contrasts of a fitted model

== Description

`design_parts()` extracts from a fitted model the pieces that `design_matrix(model, newdata = )` needs to build the model matrix of new data: the terms, the levels of the factors and the contrasts\. Each method reads the pieces that the model's own `predict()` method uses\.

== Usage

```r
#| eval: false
design_parts(model, ...)

# Default S3 method
design_parts(model, ...)

# S3 method for class 'lm'
design_parts(model, ...)

# S3 method for class 'glm'
design_parts(model, ...)

# S3 method for class 'polr'
design_parts(model, ...)

# S3 method for class 'multinom'
design_parts(model, ...)

# S3 method for class 'svyolr'
design_parts(model, ...)

# S3 method for class 'clm'
design_parts(model, ...)

# S3 method for class 'clmm'
design_parts(model, ...)

# S3 method for class 'survreg'
design_parts(model, ...)

# S3 method for class 'lm_robust'
design_parts(model, ...)

# S3 method for class 'iv_robust'
design_parts(model, ...)

# S3 method for class 'ivreg'
design_parts(model, ...)

# S3 method for class 'coxph'
design_parts(model, ...)

# S3 method for class 'lme'
design_parts(model, ...)

# S3 method for class 'gls'
design_parts(model, ...)

# S3 method for class 'merMod'
design_parts(model, ...)

# S3 method for class 'glmmTMB'
design_parts(model, component = c("cond", "zi", "disp"), ...)

# S3 method for class 'betareg'
design_parts(model, component = c("mean", "precision"), ...)

# S3 method for class 'zeroinfl'
design_parts(model, component = c("count", "zero"), ...)

# S3 method for class 'hurdle'
design_parts(model, component = c("count", "zero"), ...)

# S3 method for class 'afex_aov'
design_parts(model, ...)
```

== Arguments

#table(
  columns: 2,
  table.header([Argument], [Description]),
  [`model`], [A fitted model\.],
  [`...`], [Passed to methods\.],
  [`component`], [For models with several linear predictors, which one\. See Details\.],
)

== Value

A list with elements `terms`, `xlevels` and `contrasts`, and optionally `drop_intercept` \(see the section #emph[Methods for other classes]\)\.

== Details

Methods are provided for models from #strong[stats] \(`lm`, `glm`, and the classes that inherit from them, such as `MASS::rlm()`, `MASS::glm.nb()` and `survey::svyglm()`\), #strong[MASS] \(`polr`; `glmmPQL` through `lme`\), #strong[nnet] \(`multinom`\), #strong[nlme] \(`lme`, `gls`\), #strong[lme4] \(`merMod`\), #strong[glmmTMB], #strong[ordinal] \(`clm`, `clmm`\), #strong[survival] \(`coxph`, `survreg`, and the #strong[survey] classes that inherit from them\), #strong[survey] \(`svyolr`\), #strong[betareg], #strong[estimatr] \(`lm_robust`, `iv_robust`\), #strong[ivreg] and #strong[AER] \(`ivreg`\), #strong[pscl] \(`zeroinfl`, `hurdle`\) and #strong[afex] \(`afex_aov`\)\.

The default method, for other classes, takes the terms from #link("https://rdrr.io/pkg/stats/man/terms.html")[`stats::terms()`]; the levels from an `xlevels` element of the model, or else from #link("https://rdrr.io/pkg/stats/man/.getXlevels.html")[`stats::.getXlevels()`] on #link("https://rdrr.io/pkg/stats/man/model.frame.html")[`stats::model.frame()`]; and the contrasts from a `contrasts` element, or else from the `"contrasts"` attribute of #link("https://rdrr.io/pkg/stats/man/model.matrix.html")[`stats::model.matrix()`]\.

Some models have several linear predictors\. For these, `component` selects one: `"mean"` or `"precision"` for `betareg`; `"count"` or `"zero"` for `zeroinfl` and `hurdle`; `"cond"`, `"zi"` or `"disp"` for `glmmTMB`\. For instrumental-variable models \(`ivreg`, `iv_robust`\), the parts are those of the regressors\. For `coxph`, `strata()` and `cluster()` terms are removed and the intercept column is dropped, as in `model.matrix()` for these models\.

== See Also

#link("design_matrix.typ")[`design_matrix()`], which compiles these pieces into a reusable plan\.

== Methods for other classes

`design_parts()` is an S3 generic: a package \(or a user\) adds support for its own model class by writing a method `design_parts.<class>()` and registering it \(with `S3method(xtx::design_parts, <class>)` in the NAMESPACE of a package that only suggests #strong[xtx], or `S3method(design_parts, <class>)` when it imports it\)\. #link("design_matrix.typ")[`design_matrix()`], and therefore `design_matrix(fit, newdata = ...)`, then work for that class\. A method returns a list with elements:

- `terms`: a `terms` object, with the response or without\. Its `"predvars"` attribute must keep the parameters of data-dependent transformations such as #link("https://rdrr.io/pkg/stats/man/poly.html")[`stats::poly()`], #link("https://rdrr.io/pkg/splines/man/ns.html")[`splines::ns()`] or #link("https://rdrr.io/r/base/scale.html")[`scale()`], as the `terms` of a model frame do; its environment is where the variables and functions of the formula are looked up\.
- `xlevels`: a named list of character vectors, the levels of the factors \(as in the `xlev` argument of #link("https://rdrr.io/pkg/stats/man/model.frame.html")[`stats::model.frame()`]\), or `NULL`\. Levels of variables that are not in `terms` are ignored\.
- `contrasts`: a named list as in the `contrasts.arg` argument of #link("https://rdrr.io/pkg/stats/man/model.matrix.html")[`stats::model.matrix()`] \(contrast names or matrices\), or `NULL`\.
- `drop_intercept` \(optional\): `TRUE` removes the intercept column after the matrix is built, so that the factors are still coded as in a model with an intercept \(as `model.matrix()` does for `coxph`\)\.

The method reads these pieces from the fitted object, as its `predict()` method does; it does not need the training data\. Before writing one, check whether the default method already works: it does for most classes that store `terms`, `xlevels` and `contrasts`, or have #link("https://rdrr.io/pkg/stats/man/terms.html")[`stats::terms()`], #link("https://rdrr.io/pkg/stats/man/model.frame.html")[`stats::model.frame()`] and #link("https://rdrr.io/pkg/stats/man/model.matrix.html")[`stats::model.matrix()`] methods\.

== Examples

```r
fit <- lm(mpg ~ poly(hp, 2) + factor(cyl), data = mtcars)
parts <- design_parts(fit)
names(parts)
parts$xlevels

# a method for a new class: a fit that stores its formula and levels
toy <- structure(list(formula = mpg ~ wt + factor(gear),
                      levels = list(`factor(gear)` = c("3", "4", "5"))),
                 class = "toyfit")
design_parts.toyfit <- function(model, ...) {
    list(terms = terms(model$formula), xlevels = model$levels,
         contrasts = list(`factor(gear)` = "contr.treatment"))
}
# (in a package: S3method(xtx::design_parts, toyfit))
.S3method("design_parts", "toyfit", design_parts.toyfit)
p <- design_matrix(toy, newdata = NULL)
as.matrix(design_matrix(p, head(mtcars)))
```
