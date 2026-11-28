#import "/.calepin/calepin.typ" as calepin
#show: calepin.document

#set document(title: [Fixed effects])
#metadata((title: "Fixed effects")) <website-metadata>

#import "/benchmarks/bench.typ": bench-code
#import "/benchmarks/fixed_effects/structure.typ": dgp-figure

#calepin.setup(echo: false, eval: true, results: "verbatim", warning: false, message: false)

// Building this page never runs the benchmark: the chunks read
// results/fixed_effects.csv and results/fixed_effects_scale.csv (if
// present), written by `Rscript website/benchmarks/run.R fixed_effects`.
// Older results files lack the reference, kill_s and ref_kill_s columns:
// the chunks fall back to the run's thresholds.
#calepin.chunk(echo: false, results: "hide")[
```r
options(width = 200)
library(ggplot2)
res <- read.csv("../benchmarks/results/fixed_effects.csv")
machine <- readLines("../benchmarks/results/fixed_effects_machine.txt")
# The data generating processes, by their codes in the scripts (dgp.R)
processes <- c(balanced = "Balanced panel", random = "Random matching",
               round_robin = "Round robin", low_mobility = "Low mobility")
process <- function(x) factor(processes[x$dgp], processes)
res$process <- process(res)
res <- res[res$k == max(res$k), ]
pkgs <- c("lm_xtx" = "#18756d", "lm_xtx (DuckDB)" = "#7fbfb3",
          "feols" = "#ad481e", "dbreg (DuckDB)" = "#c79a2c")
# (lfe::felm() was dropped from the benchmark: older result files list it)
res <- res[res$package %in% names(pkgs), ]
res$package <- factor(res$package, names(pkgs))
# A placeholder instead of an empty plot, when no fit finished
no_data <- function(what) {
    ggplot() + annotate("text", x = 0, y = 0, label = what, size = 5) + theme_void()
}
# Lines by number of observations: one panel per process (two by two), one
# line per package (dashed: the data stay in DuckDB); log-log
panel_plot <- function(x, y, ylab) {
    x$value <- x[[y]]
    dashed <- ifelse(grepl("DuckDB", names(pkgs)), "22", "solid")
    ggplot(x, aes(n, value, colour = package, linetype = package)) +
        geom_line(linewidth = 0.9) + geom_point(size = 2.2) +
        facet_wrap(~ process, ncol = 2) +
        scale_x_log10(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
        scale_y_log10(labels = scales::label_number(drop0trailing = TRUE, big.mark = ",")) +
        scale_colour_manual(values = pkgs, drop = FALSE) +
        scale_linetype_manual(values = setNames(dashed, names(pkgs)), drop = FALSE) +
        labs(x = "Observations", y = ylab, colour = NULL, linetype = NULL) +
        theme_minimal(base_size = 13) +
        theme(legend.position = "top", legend.justification = "left",
              legend.key.width = unit(2.2, "lines"),
              panel.grid.minor = element_blank(),
              panel.spacing = unit(1.3, "lines"),
              strip.text = element_text(face = "bold", hjust = 0),
              axis.title = element_text(colour = "grey35"))
}
provisional <- any(res$provisional)
# The thresholds of the run, saved with the results (older results files
# predate that line and used the defaults, 120 s and 1800 s)
limits <- regmatches(machine, regexec("^Limits: .*after ([0-9.e+]+) s .*after ([0-9.e+]+) s", machine))
limits <- Filter(length, limits)
limit <- if (length(limits)) as.numeric(limits[[1]][2]) else 120
kill <- if (length(limits)) as.numeric(limits[[1]][3]) else 1800
duration <- function(s) {
    if (s >= 60 && s %% 60 == 0) {
        m <- s / 60
        if (m >= 60 && m %% 60 == 0) paste(m / 60, if (m == 60) "hour" else "hours")
        else paste(m, if (m == 1) "minute" else "minutes")
    } else paste(format(s), "seconds")
}
durations <- function(s) vapply(s, duration, "")
# Results saved before the cached data and separate references (no kill_s
# column) were measured with the previous protocol, and every cell was
# stopped after `kill` seconds, references included
old_protocol <- is.null(res$kill_s)
if (old_protocol) res$kill_s <- kill
if (is.null(res$ref_kill_s)) res$ref_kill_s <- kill
# Cells stopped sooner than the others (known not to converge)
short <- unique(res[res$kill_s < kill, c("package", "process", "kill_s")])
fmt_n <- function(n) format(n, big.mark = ",", scientific = FALSE, trim = TRUE)
# The at-scale grid, saved in its own file once it has run
scale_file <- "../benchmarks/results/fixed_effects_scale.csv"
sc <- if (file.exists(scale_file)) read.csv(scale_file)
if (!is.null(sc)) {
    sc$process <- droplevels(process(sc))
    sc$package <- factor(sc$package, names(pkgs))
}
scale_machine <- if (file.exists("../benchmarks/results/fixed_effects_scale_machine.txt")) readLines("../benchmarks/results/fixed_effects_scale_machine.txt")
```
]

#title()

= Usage

Fixed effects go after a `|` in the formula, as in `fixest`. The same model, fitted with both packages on a dataset that comes with `fixest`:

#calepin.chunk(echo: true)[
```r
library(fixest)
library(xtx.fit)
data(trade)

fit1 <- feols(log(Euros) ~ log(dist_km) | Origin + Destination + Year, trade)
fit2 <- lm_xtx(log(Euros) ~ log(dist_km) | Origin + Destination + Year, trade)

coef(fit1)
coef(fit2)
all.equal(c(vcov(fit1)), c(vcov(fit2)))
```
]

See #link("../reference-fit/lm_xtx.html")[`lm_xtx()`] for the options (clustered standard errors, predictions, the values of the fixed effects) and #link("glm_xtx.html")[`glm_xtx()`], which takes the same formulas, for generalized linear models.

= Benchmark

How fast is `lm_xtx()` with fixed effects, compared with the other R packages that fit linear models with high-dimensional fixed effects? This page replicates the fixed-effects OLS benchmark of #link("https://github.com/py-econometrics/pyfixest")[pyfixest], itself adapted from Kyle Butts' #link("https://github.com/kylebutts/fixest_benchmarks")[fixest_benchmarks], for R packages only.

#calepin.chunk(echo: false, results: "typst")[
```r
if (provisional) cat("*Provisional results.* These numbers come from a reduced grid (", paste(fmt_n(sort(unique(res$n))), collapse = ", "), " observations), run while validating the scripts. The full grid, up to 1 million observations, will replace them.", if (old_protocol) " They were also measured with an earlier version of the protocol below: each process simulated its own data, the DuckDB fits read an in-database table instead of the Parquet file, and the reference fit ran in the `feols()` process, after the timed fits.", "\n", sep = "")
```
]

= Data generating processes

Every dataset is simulated from the same model:
$ y = sum_(j=1)^k x_j beta_j + alpha_"indiv" + delta_"year" + gamma_"group" + epsilon. $

The balanced panel has no group effect. The processes differ in how individuals are matched with years and groups.

== Balanced panel

States observed every year. Each state has one row per year, without gaps: every state appears in every year.

#figure(
  dgp-figure("balanced"),
  caption: [A balanced panel of states and years: every state is observed in every year.],
)

- States or countries observed every year.
- Stores with sales recorded every week.
- The subjects of an experiment, measured in every round.

#calepin.chunk(echo: true, eval: false)[
```r
n_indiv <- 10000; n_year <- 10       # indiv_id: the state
d <- data.frame(indiv_id = rep(1:n_indiv, each = n_year), year = rep(1:n_year, n_indiv))
d$x1 <- rnorm(nrow(d))
d$y <- d$x1 + rnorm(n_indiv)[d$indiv_id] + rnorm(n_year)[d$year] + rnorm(nrow(d))
lm_xtx(y ~ x1 | indiv_id + year, d)
```
]

== Random matching

Shoppers and the products they buy. Every week, each shopper buys one product, picked at random among all the products, whatever they bought the week before. Over time, a shopper goes through many different products.

#figure(
  dgp-figure("random"),
  caption: [Random matching of shoppers and products: each shopper buys a product picked at random every week, without any pattern.],
)

- Shoppers and the products they buy, or viewers and the films they rate.
- Freelancers or temporary workers and the many clients they work for.
- Callers and the call-centre agents who answer them.

#calepin.chunk(echo: true, eval: false)[
```r
n_indiv <- 10000; n_year <- 10; n_group <- 435   # indiv_id: the shopper, year: the week, group_id: the product
d <- data.frame(indiv_id = rep(1:n_indiv, each = n_year), year = rep(1:n_year, n_indiv))
d$group_id <- sample(n_group, nrow(d), replace = TRUE)
d$x1 <- rnorm(nrow(d))
d$y <- d$x1 + rnorm(n_indiv)[d$indiv_id] + rnorm(n_year)[d$year] + rnorm(n_group)[d$group_id] + rnorm(nrow(d))
lm_xtx(y ~ x1 | indiv_id + year + group_id, d)
```
]

== Round robin

Students and their teachers. Every year, each student moves up to the next teacher in line, as a class moves up to the teacher of the grade above. Over time, a student has a run of consecutive teachers. Students who start with different teachers follow the same line, a few teachers ahead or behind. In the simulated data the line closes into a loop: after the last teacher comes the first.

#figure(
  dgp-figure("round-robin"),
  caption: [Round robin of students and teachers: each student (row) has a run of consecutive teachers, one a year, and starts a few teachers after the student above. The last student goes on with the first teachers.],
)

- Students passed from the teachers of one grade to those of the next.
- Medical residents rotating through wards, or staff on a fixed rota.
- Surveys that rotate overlapping cohorts.

#calepin.chunk(echo: true, eval: false)[
```r
n_indiv <- 10000; n_year <- 10; n_group <- 435   # indiv_id: the student, group_id: the teacher
d <- data.frame(indiv_id = rep(1:n_indiv, each = n_year), year = rep(1:n_year, n_indiv))
d$group_id <- rep_len(1:n_group, nrow(d))
d$x1 <- rnorm(nrow(d))
d$y <- d$x1 + rnorm(n_indiv)[d$indiv_id] + rnorm(n_year)[d$year] + rnorm(n_group)[d$group_id] + rnorm(nrow(d))
lm_xtx(y ~ x1 | indiv_id + year + group_id, d)
```
]

== Low mobility

Workers and their firms. Workers stay at the same firm from one year to the next and change only now and then: each year, a worker leaves for another firm, picked at random, with a small probability.

#figure(
  dgp-figure("low-mobility"),
  caption: [Low mobility of workers between firms: workers spend most of their years at one firm (dark cells) and change only now and then.],
)

- Workers and firms in matched employer-employee data.
- Patients and their physicians.
- Students and their schools.

#calepin.chunk(echo: true, eval: false)[
```r
n_indiv <- 10000; n_year <- 10; n_group <- 435   # indiv_id: the worker, group_id: the firm
d <- data.frame(indiv_id = rep(1:n_indiv, each = n_year), year = rep(1:n_year, n_indiv))
# a new spell in a worker's first year, and when they move (8% of years)
spell <- cumsum(d$year == 1 | runif(nrow(d)) < 0.08)
d$group_id <- sample(n_group, max(spell), replace = TRUE)[spell]
d$x1 <- rnorm(nrow(d))
d$y <- d$x1 + rnorm(n_indiv)[d$indiv_id] + rnorm(n_year)[d$year] + rnorm(n_group)[d$group_id] + rnorm(nrow(d))
lm_xtx(y ~ x1 | indiv_id + year + group_id, d)
```
]

= Packages

In memory:

- `xtx.fit::lm_xtx()`
- `fixest::feols()`

In database (DuckDB):

- `xtx.fit::lm_xtx()`
- `dbreg::dbreg()`

= Algorithms

- `fixest` and `pyfixest` demean the data by alternating projections: they subtract the means of each fixed effect in turn from the outcome and from every regressor, sweeping over the rows until nothing changes, then run least squares on the demeaned data. The number of sweeps depends on how individuals are matched with groups: one on a balanced panel, a few with random matching, many with round robin.
- `xtx.fit` never demeans. A fixed number of passes over the rows computes sums (cross-products of the regressors, their sums within each level, and counts of the levels that occur together), then a linear system over the levels of the fixed effects is solved: the largest fixed effect is eliminated exactly, the others by a direct factorization (up to 1,000 levels) or by conjugate gradients. `options(xtx.fe_precond = "jacobi")` or `"ichol"` forces a preconditioner, and `xtx.fe_form_max` sets when the reduced system is formed (`0` never, `Inf` always).

The time of the first grows with the number of sweeps, so with the data generating process; the time of the second depends on it much less. See #link("../architecture.html")[Architecture] for details.

= Time

Median time of one fit with ten regressors, in seconds, by number of observations. Each panel is a data generating process, each line a package (dashed when the data stay in DuckDB). Both axes are on log scales. A line that stops early means that the fit did not finish at the next size (see Protocol).

#calepin.chunk(
  echo: false,
  results: "render",
  fig-width: "100%",
  fig-device-width: 10,
  fig-device-height: 8.5,
  fig-alt-text: "Median fitting time with ten regressors by number of observations, log-log scale, for lm_xtx and feols in memory and lm_xtx and dbreg on DuckDB, one panel per data generating process: balanced panel, random matching, round robin, low mobility.",
)[
```r
x <- res[!is.na(res$seconds), ]
if (nrow(x) == 0) no_data("No timed fits.") else panel_plot(x, "seconds", "Median seconds")
```
]

= Memory

Extra peak memory of one fit with ten regressors, in MB, by number of observations, with the same panels and lines.

#calepin.chunk(
  echo: false,
  results: "render",
  fig-width: "100%",
  fig-device-width: 10,
  fig-device-height: 8.5,
  fig-alt-text: "Extra peak resident memory of one fit with ten regressors by number of observations, log-log scale, for lm_xtx and feols in memory and lm_xtx and dbreg on DuckDB, one panel per data generating process.",
)[
```r
x <- res[!is.na(res$mem_mb), ]
x$mem <- pmax(x$mem_mb, 0.1)
if (nrow(x) == 0) no_data("No memory measurements.") else panel_plot(x, "mem", "Extra peak memory, MB")
```
]

= At scale

The balanced panel and random matching with ten regressors, at 10 and 50 million observations, for three fits only: `lm_xtx()` in memory, `lm_xtx()` on DuckDB reading the Parquet file, and `feols()`. The protocol is the same, with a time limit of its own, and without reference coefficients.

#calepin.chunk(echo: false, results: "typst")[
```r
if (is.null(sc)) {
    cat("_Not run yet._ These results will come with the full run of the benchmark.\n")
} else {
    if (any(sc$provisional)) cat("*Provisional results* (reduced grid: ", paste(fmt_n(sort(unique(sc$n))), collapse = ", "), " observations).\n\n", sep = "")
    m <- sc[sc$status %in% c("skipped", "killed", "error"), ]
    if (nrow(m) > 0) {
        cat("Missing points: ", paste(sprintf("%s, %s, %s rows: %s", m$package, tolower(m$process), fmt_n(m$n),
            ifelse(m$status == "killed", paste("stopped after", durations(m$kill_s)),
                   ifelse(m$status == "skipped", "not run (too slow at a smaller size)", "failed"))), collapse = "; "),
            ".\n", sep = "")
    }
}
```
]

#calepin.chunk(
  echo: false,
  results: "render",
  fig-width: "100%",
  fig-device-width: 10,
  fig-device-height: 6,
  fig-alt-text: "Median fitting time and extra peak memory at 10 and 50 million observations, log scales, for lm_xtx in memory, lm_xtx on DuckDB and feols, on the balanced panel and with random matching (columns).",
)[
```r
# No figure until the at-scale grid has run
if (!is.null(sc)) {
    x <- sc[!is.na(sc$seconds), ]
    long <- if (nrow(x) > 0) rbind(
        data.frame(x[c("package", "n", "process")], what = "Median seconds", value = x$seconds),
        data.frame(x[c("package", "n", "process")], what = "Extra peak memory, MB", value = pmax(x$mem_mb, 0.1)))
    if (!is.null(long)) long$what <- factor(long$what, c("Median seconds", "Extra peak memory, MB"))
    print(if (nrow(x) == 0) no_data("No fit finished at scale.") else
    ggplot(long, aes(n, value, colour = package)) +
        geom_line(linewidth = 0.8) + geom_point(size = 2) +
        facet_grid(what ~ process, scales = "free_y", switch = "y") +
        scale_x_log10(labels = scales::label_number(scale_cut = scales::cut_short_scale()),
                      breaks = sort(unique(x$n))) +
        scale_y_log10(labels = scales::label_number(drop0trailing = TRUE, big.mark = ",")) +
        scale_colour_manual(values = pkgs[levels(droplevels(x$package))]) +
        labs(x = "Observations (log scale)", y = NULL, colour = NULL) +
        theme_minimal(base_size = 13) +
        theme(legend.position = "bottom", strip.placement = "outside"))
}
```
]

#calepin.chunk(echo: false, results: "verbatim")[
```r
if (!is.null(scale_machine)) cat(scale_machine, sep = "\n")
```
]

= Accuracy

Largest absolute difference between each package's coefficients and those of `feols()` at a tight tolerance, over every fit.

#calepin.chunk(echo: false, results: "typst")[
```r
# One row per package in the results, whether or not it has a comparison
present <- levels(droplevels(res$package))
acc <- data.frame(Package = present)
acc$diff <- vapply(present, function(p) {
    d <- res$max_coef_diff[res$package == p & !is.na(res$max_coef_diff)]
    if (length(d)) format(max(d), digits = 2)
    else if (old_protocol && !"feols" %in% present) "not available (feols reference not run)"
    else "not available (no reference for these fits)"
}, "")
names(acc)[2] <- "Largest coefficient difference"
tinytable::tt(acc)
```
]

= Protocol

- *Data.* 1,000 to 1 million observations with ten regressors: for $n$ observations, $n slash 10$ individuals over 10 years and $n slash 230$ groups. Each dataset is simulated once and saved as a Parquet file.
- *Models.* `y ~ x1 + ... + x10 | indiv_id + year` on the balanced panel, with `+ group_id` on the other processes; iid standard errors.
- *One thread.* Every package and the BLAS run on one thread, and the whole run is pinned to one CPU core.
- *Time.* Each fit runs in a fresh R process: one warmup fit, then the median of three. In-memory fits read the data before the clock starts; DuckDB fits query the Parquet file, so their time includes reading it.
- *Memory.* Peak resident memory during the warmup fit, over its level just before. The data frame of the in-memory fits is not counted.
#calepin.chunk(echo: false, results: "typst")[
```r
stopped_sooner <- if (nrow(short) > 0) {
    paste0(" (", paste(sprintf("%s for %s on %s", durations(short$kill_s), short$package, tolower(short$process)),
                       collapse = "; "), ")")
}
cat("- *Limits.* A fit that takes more than ", duration(limit), " is timed once, and the package is not run at larger sizes. A fit is stopped after ",
    duration(kill), stopped_sooner, ". `dbreg()` runs with its limit of alternating projections raised from 100 to 10,000.\n", sep = "")
cat("- *Accuracy.* Coefficients are compared with those of `feols()` at a tight tolerance (`fixef.tol = 1e-11`), fitted outside the timed runs.\n", sep = "")
```
]

#calepin.chunk(echo: false, results: "verbatim")[
```r
cat(machine, sep = "\n")
```
]

= Code

The benchmark is `website/benchmarks/fixed_effects.R`. It simulates each dataset once with `fixed_effects/simulate.R` and this function:

#bench-code("fixed_effects/dgp")

Each fit runs in a fresh R process with this script (the reference fits run separately, in `fixed_effects/reference.R`):

#bench-code("fixed_effects/cell")
