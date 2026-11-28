#import "/.calepin/calepin.typ" as calepin
#show: calepin.document

#set document(title: [Out-of-core data])
#metadata((title: "Out-of-core data")) <website-metadata>

#calepin.setup(echo: true, eval: true, results: "verbatim", warning: false, message: false)

#calepin.chunk(echo: false, results: "hide")[
```r
options(width = 200)
```
]

#title()

When the data do not fit in memory, they can stay in Parquet files or in a DuckDB table. Here `d` is a data frame and `"d.parquet"` a file with the same rows. The code does not change: only the data argument does. The rows are never read into R: products with the design run in the database, and only small results come back. This needs the DBI and duckdb packages.

#calepin.setup(fenced-chunks: false)
#table(
  columns: (1fr, 1fr),
  [*In memory*], [*Out of core*],
  [```r
data <- d
```],
  [```r
data <- "d.parquet"
```],
  [```r
f <- y ~ x + g
D <- design_matrix(f, data)
crossprod(D)
fit <- lm_xtx(f, data)
coef(fit)
```],
  [```r
f <- y ~ x + g
D <- design_matrix(f, data)
crossprod(D)
fit <- lm_xtx(f, data)
coef(fit)
```],
)
#calepin.setup(fenced-chunks: true)

The same code, run on both, gives the same results:

```r
library(xtx)
library(xtx.fit)
set.seed(1)
d <- data.frame(x = rnorm(1000), g = sample(c("a", "b", "c"), 1000, TRUE))
d$y <- d$x + (d$g == "b") + rnorm(1000)

# the same rows, in a Parquet file
con <- DBI::dbConnect(duckdb::duckdb())
DBI::dbWriteTable(con, "d", d)
pq <- file.path(tempdir(), "d.parquet")
invisible(DBI::dbExecute(con, sprintf("COPY d TO '%s' (FORMAT PARQUET)", pq)))

run <- function(data) {
    f <- y ~ x + g
    D <- design_matrix(f, data)
    fit <- lm_xtx(f, data)
    list(crossprod = crossprod(D), coef = coef(fit), vcov = vcov(fit))
}
all.equal(run(d), run(pq))
run(pq)$coef
```

A DuckDB table takes a connection and the name of the table (or a `sql` query) in place of the file:

```r
fit <- lm_xtx(y ~ x + g, con, table = "d")
all.equal(coef(fit), run(d)$coef)
```

#calepin.chunk(echo: false, results: "hide")[
```r
DBI::dbDisconnect(con)
```
]

The recipes show more on a database: #link("recipes/algebra.html#duckdb-and-parquet")[designs, algebra and columns], #link("fit/lm_xtx.html#on-a-database")[linear models], #link("fit/glm_xtx.html#on-a-database")[generalized linear models] and #link("fit/standard_errors.html#robust-standard-errors-on-a-database")[robust standard errors].

= Differences

- `model.frame()` is not available for a remote fit, and `model` is ignored: the data stay in the database (`fit$data$design` is the remote design of the fit).
- Results with one value per row are tables in the database, not vectors: `fitted()` (`.fitted`, `.resid`), `predict()` (`.fitted`; `se.fit` as a column `.se.fit` with attributes `df` and `residual.scale`; `.nonestim` instead of the attribute `"non-estim"`), and `residuals()`, a remote design.
- These tables have no row order: rows are matched to the data through its columns, which every table keeps.
- `compress = "auto"` decides from other evidence than in memory; the results are the same either way.
- Remote designs keep only the algebra promise: results agree with the in-memory computation up to rounding. `as.matrix()` refuses to collect more than 100,000 rows unless `force = TRUE`. Only DuckDB connections are supported. With infinite values in the data, `NaN` and `Inf` can land in other places than in the dense computation. See `?design_remote` for the details, and `?xtx.fit::remote_fits` for the fits on remote data.  #calepin.chunk(echo: false, results: "hide")[ ```r DBI::dbDisconnect(con) ``` ]
