#import "/.calepin/calepin.typ" as calepin
#show: calepin.document

#set document(title: [Benchmarks])
#metadata((title: "Benchmarks")) <website-metadata>

#import "/benchmarks/bench.typ": bench, bench-table, bench-code

#calepin.setup(echo: true, eval: true, results: "verbatim", warning: false, message: false)

#title()

How much time does it save to skip the dense model matrix? This page times six common operations twice: once with base R, which builds the full matrix first, and once with `xtx`, which works from the design matrix instead. The data have 200,000 rows, the model matrix has 119 columns, and everything runs on one thread. The answer depends on what the columns are, so there are two models: one with factors, one with numeric variables only.

Each operation is timed in two situations:

- *From a data frame.* The matrix is built and then used once, as when predicting on a new batch of data.
- *From a prebuilt design matrix.* The matrix already exists and only the operation is timed, as when the same matrix is used again and again. Base R starts from its own prebuilt dense matrix.

A third panel times model fits: `lm()` and `glm()` against `lm_xtx()` and `glm_xtx()` from `xtx.fit`, and, for the model with factors, predictions on new data, their standard errors and leverage from a fitted `lm()`, against the same quantities computed on the design matrix.

Shorter bars are better.

#calepin.chunk(echo: false, results: "hide")[
```r
library(ggplot2)
ms <- function(x) {
    x <- trimws(x)
    as.numeric(sub("[^0-9.]+$", "", x)) * c(ns = 1e-6, "\u00b5s" = 1e-3, us = 1e-3, ms = 1, s = 1e3, m = 6e4)[sub("^[0-9.]+", "", x)]
}
starts <- c("From a data frame", "From a prebuilt design matrix", "Model fits, from a data frame")
# One bar per operation, starting point and package
bar_chart <- function(res) {
    res$operation <- factor(res$operation, rev(unique(res$operation)))
    res$start <- factor(res$start, starts)
    res$package <- factor(res$package, c("xtx, design matrix", "Base R, dense matrix"))
    res$label <- paste(format(signif(res$ms, 3), big.mark = ",", drop0trailing = TRUE, trim = TRUE), "ms")
    ggplot(res, aes(ms, operation, fill = package)) +
        geom_col(position = position_dodge(width = 0.75), width = 0.7) +
        geom_text(aes(label = label), position = position_dodge(width = 0.75), hjust = -0.1, size = 4.2, colour = "grey30") +
        facet_wrap(~ start, scales = "free", ncol = 1) +
        scale_x_continuous(expand = expansion(mult = c(0, 0.32)), labels = scales::label_comma()) +
        scale_fill_manual(values = c("xtx, design matrix" = "#18756d", "Base R, dense matrix" = "#ad481e"),
                          breaks = c("Base R, dense matrix", "xtx, design matrix")) +
        labs(x = "Median time, milliseconds", y = NULL, fill = NULL) +
        theme_minimal(base_size = 16) +
        theme(legend.position = "top", legend.justification = "left",
              panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
              panel.spacing = unit(1.5, "lines"),
              strip.text = element_text(face = "bold", hjust = 0),
              axis.text.y = element_text(family = "mono", colour = "grey20"),
              axis.title = element_text(colour = "grey35"))
}
```
]

= Numeric and factors

The model is `y ~ x * f1 + f2`, with a factor of 50 levels and one of 20. All but two of its 119 columns are factor dummies and their interactions with `x`, so each row has only a few non-zero entries, which `xtx` skips.

#calepin.chunk(
  echo: false,
  results: "render",
  fig-width: "100%",
  fig-device-width: 10,
  fig-device-height: 13,
  fig-alt-text: "Bar chart of the median time in milliseconds of six operations on a model with factors (building the matrix, matrix product, crossprod with a vector, crossprod, weighted group sums, quadratic diagonal), for base R on the dense model matrix and for xtx on the design matrix, from a data frame and from a prebuilt design matrix; and of model fits and post-estimation (lm(), glm(), predictions, their standard errors, leverage) for base R and for xtx and xtx.fit.",
)[
```r
# Each results file has base R then xtx, from a data frame, then (but for
# the construction) the same two from a prebuilt matrix
ops <- c(construction = "Build the matrix", matmul = "%*%", crossprod_w = "crossprod(D, w)",
         crossprod = "crossprod(D)", colsums = "colSums(w * D, group = )", quadratic = "quadratic_diagonal()")
res <- do.call(rbind, lapply(names(ops), function(f) {
    r <- read.csv(file.path("benchmarks", "results", paste0(f, ".csv")), colClasses = "character")
    data.frame(operation = ops[[f]], ms = ms(r$median),
               package = rep(c("Base R, dense matrix", "xtx, design matrix"), length.out = nrow(r)),
               start = rep(starts[1:2], each = 2)[seq_len(nrow(r))])
}))
# Model fits and what follows them: base R (first row of each file) against
# xtx and xtx.fit (last row); leverage.csv has a middle row, base R's
# formula on the prebuilt matrix, left out here
fits <- c(lm = "lm()", glm = "glm(, binomial)", predict = "predict(fit, newdata)",
          predict_se = "predict(, se.fit = TRUE)", leverage = "hatvalues()")
res <- rbind(res, do.call(rbind, lapply(names(fits), function(f) {
    r <- read.csv(file.path("benchmarks", "results", paste0(f, ".csv")), colClasses = "character")
    data.frame(operation = fits[[f]], ms = ms(r$median[c(1, nrow(r))]),
               package = c("Base R, dense matrix", "xtx, design matrix"), start = starts[3])
})))
bar_chart(res)
```
]

For `crossprod(D)`, the design matrix adds up the products of each row's few non-zero entries, instead of all `p^2 / 2` products of the dense row. For `quadratic_diagonal()`, base R builds an `n` by `p` intermediate, where the kernel touches only the non-zero entries of each row. Examples and benchmarks of threads, reuse on new data, sparse matrices and repeated rows are in the #link("recipes/algebra.html")[Algebra recipes].

= Numeric only

The model has 118 numeric variables and an intercept: as many columns, but no zeros to skip. For the matrix operations, `xtx` has nothing to gain: the two are close on every operation, each a little ahead on some. The model fits are another matter: `lm_xtx()` and `glm_xtx()` remain several times faster than `lm()` and `glm()` (the third panel), because they factor the small `p` by `p` cross-product instead of the `n` by `p` matrix. `glm_xtx()` does that once per iteration.

#calepin.chunk(
  echo: false,
  results: "render",
  fig-width: "100%",
  fig-device-width: 10,
  fig-device-height: 11,
  fig-alt-text: "Bar chart of the median time in milliseconds of the same six operations on a model with numeric variables only, for base R on the dense model matrix and for xtx on the design matrix, from a data frame and from a prebuilt design matrix; and of lm() and glm() fits against lm_xtx() and glm_xtx().",
)[
```r
bar_chart(read.csv(file.path("benchmarks", "results", "numeric.csv")))
```
]

Two things still set `xtx` apart:

- *The dense matrix is always available.* `as.matrix(design_matrix())` returns the same matrix as `model.matrix()`, bit for bit, in about the same time (the first pair of bars). Code that needs the dense matrix of a numeric model can build it at no extra cost.
- *The same code runs out of core.* A design can be built on a DuckDB table or on Parquet files that do not fit in memory, where there is no dense matrix to fall back on: see #link("out-of-core.html")[Out-of-core data].

Other pages have their own benchmarks: #link("fit/fixed_effects.html")[fixed effects], #link("fit/lm_xtx.html#least-squares")[least squares], #link("fit/leverage.html")[leverage], #link("fit/standard_errors.html")[standard errors] and #link("recipes/prediction.html")[predictions].

= Setup

The benchmarks of this site run on this data, with one thread unless stated otherwise. `X` and `DX` are prebuilt dense and design matrices, and `same()` and `close()` are the checks for exact and approximate results. Each benchmark runs in a fresh R process on top of this setup; the results shown were saved by `make bench` (see #link("https://github.com/vincentarelbundock/xtx/tree/main/website/benchmarks")[`website/benchmarks`]).

The benchmark tables report minimum and median elapsed time and `mem_alloc`, the total bytes allocated on the R heap during one evaluation. Reproducing the memory column requires an R build with `capabilities("profmem")` enabled. For fair comparisons with functions that call `model.matrix()`, the benchmarks set `options(xtx.skip_na_copy = FALSE)`, so that the model frame is built in exactly the same way (by default, `xtx` skips a copy of the frame that `na.omit()` makes even when no row has a missing value).

= Code

These are the scripts themselves, read from `website/benchmarks/` when the site is built. `run.R` runs each benchmark in a fresh R process, on top of `setup.R`, and saves its results as a CSV file, which the charts above read.

`setup.R`:

#bench-code("setup")

`construction.R`:

#bench-code("construction")

`matmul.R`:

#bench-code("matmul")

`crossprod_w.R`:

#bench-code("crossprod_w")

`crossprod.R`:

#bench-code("crossprod")

`colsums.R`:

#bench-code("colsums")

`quadratic.R`:

#bench-code("quadratic")

`lm.R`:

#bench-code("lm")

`glm.R`:

#bench-code("glm")

`predict.R`:

#bench-code("predict")

`predict_se.R`:

#bench-code("predict_se")

`leverage.R`:

#bench-code("leverage")

`numeric.R`:

#bench-code("numeric")

`run.R`:

#bench-code("run")
