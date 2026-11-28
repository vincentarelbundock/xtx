#import "/.calepin/calepin.typ" as calepin
#show: calepin.document

#set document(title: [Prediction])
#metadata((title: "Prediction")) <website-metadata>

#import "/benchmarks/bench.typ": bench, bench-table, bench-code

#calepin.setup(echo: true, eval: true, results: "verbatim", warning: false, message: false)

#calepin.chunk(echo: false, results: "hide")[
```r
options(width = 200)
```
]

#title()

These recipes compute predictions from a fitted model and their standard errors. Both build the design for the new data from the model's own encoding, so they rely on #link("algebra.html#new-data")[reusing a design on new data]. Each one checks a small example against base R, then times both on 200,000 rows and 119 columns, on one thread.

```r
library(xtx)
d <- data.frame(
    y = c(0.3, -1.2, 0.8, 1.5, -0.4, 0.9, -2.1, 0.2),
    x = c(-1.5, 0.2, 1.1, -0.3, 2.0, 0.7, -0.9, 1.4),
    f = factor(c("a", "b", "c", "a", "b", "c", "a", "b")),
    cl = c(1, 1, 2, 2, 3, 3, 4, 4))
```

The benchmarks run on top of the #link("../benchmarks.html#setup")[benchmark setup], in a fresh R process each.

= Predictions

The predictions of a linear model are `X b`, the design matrix times the coefficients.

```r
fast_predict <- function(fit, newdata) {
    D <- design_matrix(fit, newdata = newdata)
    drop(D %*% coef(fit))
}

fit <- lm(y ~ x + f, d)
a <- predict(fit, d)
b <- fast_predict(fit, d)
all.equal(a, b, check.attributes = FALSE)
```

For repeated batches, prepare the specification once and reuse it:

#bench("predict")

= Standard errors

The standard errors of the predictions are the square roots of the diagonal of `X V X'`, where `V` is the variance-covariance matrix of the coefficients.

```r
base_predict_se <- function(fit, newdata) {
    predict(fit, newdata, se.fit = TRUE)$se.fit
}
fast_predict_se <- function(fit, newdata) {
    D <- design_matrix(fit, newdata = newdata)
    sqrt(quadratic_diagonal(D, vcov(fit)))
}

fit <- lm(y ~ x + f, d)
a <- base_predict_se(fit, d)
b <- fast_predict_se(fit, d)
all.equal(a, b, check.attributes = FALSE)
```

For this full-rank linear model, `predict(..., se.fit = TRUE)` returns predictions and standard errors. The helper computes only the standard errors and reuses a prepared design specification.

#bench("predict_se")
