#import "/.calepin/calepin.typ" as calepin
#show: calepin.document

#set document(title: [xtx])
#metadata((title: "Introduction")) <website-metadata>

#import "/benchmarks/bench.typ": bench-table

#calepin.setup(echo: true, eval: true, results: "verbatim", warning: false, message: false)


#calepin.chunk(echo: false, results: "hide")[
```r
options(width = 200)
```
]

#title()

The `xtx` package builds model matrices that are identical, bit for bit, to those returned by `stats::model.matrix()`, and computes products, sums and quadratic forms of them without forming the dense matrix. It is an extremely fast way to compute design matrices and compute standard `X'X`-style transformations, both in and out of memory.

- *In or out of memory.* The same calls run on data frames, DuckDB tables and Parquet files, with the same results in memory and in the database.
- *Speed.* Products, least squares inputs, predictions and robust variances are computed directly from the factor codes and numeric variables, which often saves an order of magnitude or more in time.
- *Memory.* The design stores only the codes, the numeric variables and the contrast matrices, which can take much less memory than the dense matrix.

Use `xtx` when the model matrix is only an intermediate input to something else. When the dense matrix itself is what you need, there is no step to skip, and the gains are small. If you already have the dense matrix `X`, use it directly.

= Installation

```r
#| eval: false
# install.packages("remotes")
remotes::install_github("vincentarelbundock/xtx")
```

= Example

Build a design with an intercept, a numeric predictor, and a two-level factor. Multiplying it by a coefficient vector gives one linear prediction per row.

```r
library(xtx)

d <- data.frame(x = c(1, 2, 3), group = factor(c("a", "b", "a")))
D <- design_matrix(~ x + group, d)
beta <- c("(Intercept)" = 10, x = 2, groupb = -1)
drop(D %*% beta)
```

More examples: #link("recipes/algebra.html")[Algebra] (products, `crossprod()`, `qr()`, sparse matrices, new data), #link("recipes/prediction.html")[Prediction], #link("fit/lm_xtx.html")[Model fitting], #link("out-of-core.html")[Out-of-core data] (DuckDB and Parquet), and #link("recipes/package-developers.html")[Estimators] (building estimators from the primitives).

= Performance

The plot compares prediction time on new batches of rows for two numeric predictors and for a numeric predictor interacting with a 20-level factor. Both methods reuse the fitted model's terms, factor coding, and coefficients. Timing starts after one-time setup and includes multiplying by the coefficients. The package uses one thread.

#calepin.chunk(
  echo: false,
  results: "render",
  fig-width: "100%",
  fig-device-width: 10,
  fig-device-height: 5,
  fig-alt-text: "Median prediction time by batch size for base R and design_matrix().",
)[
```r
library(ggplot2)
times <- read.csv("benchmarks/results/batch_timing.csv")
ggplot(times, aes(rows, ms, colour = method)) +
    geom_line(linewidth = 1) + geom_point(size = 2.5) +
    facet_wrap(~ model) +
    expand_limits(y = 0) +
    scale_x_continuous(labels = scales::label_comma()) +
    scale_colour_manual(values = c("model.matrix()" = "#ad481e", "design_matrix()" = "#18756d")) +
    labs(title = "Predicting on new batches of rows",
         x = "Rows per batch", y = "Median prediction time (ms)", colour = NULL) +
    theme_minimal(base_size = 14) + theme(legend.position = "bottom")
```
]

The table covers time and memory for other tasks, on data with 200,000 rows and 119 columns, on one thread:

#bench-table("summary", code-column: false)

When you do need the dense matrix, `xtx` returns the identical matrix at about the same speed on one thread and about twice as fast on four threads, but it takes the same memory. The #link("operations.html")[Operations] page and the #link("fit/fixed_effects.html")[fixed effects benchmark] give the details.

= Guarantees

#table(
  columns: (auto, 1fr),
  align: (left, left),
  [*Exact*], [Model matrices pass `identical(a, b, num.eq = FALSE, single.NA = FALSE, attrib.as.set = FALSE)` against `model.matrix()`: same bits (signed zeros and `NA`/`NaN` payloads included), attributes in the same order, and the same warnings and errors, for any number of threads.],
  [*New data*], [A design built from a fitted model uses the same columns, factor coding, and transformations that the model uses for prediction.],
  [*Algebra*], [Products and sums computed without the dense matrix agree with the dense computation up to floating point rounding, with `NA` and `NaN` in the same places.],
  [*Sparse*], [Sparse model matrices hold the same non-zero values as the dense ones, bit for bit.],
)

The test suite checks these on hundreds of random formulas and data sets.

= Model fitting

The companion package `xtx.fit` (#link("https://github.com/vincentarelbundock/xtx.fit")[GitHub]) fits linear and generalized linear models on these designs, with or without absorbed fixed effects, in memory or in DuckDB. `lm_xtx()` and `glm_xtx()` give the results of `lm()`, `glm()` and `fixest` without the dense model matrix. It is built on the exported functions of `xtx`. Install it with `remotes::install_github("vincentarelbundock/xtx.fit")`.

= Setup

The benchmarks on this site ran on:

#raw(read("/benchmarks/results/machine.txt"), block: true)
