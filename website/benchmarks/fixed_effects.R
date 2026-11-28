# Fixed-effects OLS: lm_xtx() against other R packages, after the
# benchmark of pyfixest (itself adapted from Kyle Butts' fixest_benchmarks).
# lm_xtx() is in the package xtx.fit, which must be installed (with xtx).
# The full run is one command, pinned to one core:
#
#   taskset -c 2 Rscript website/benchmarks/run.R fixed_effects
#
# Two grids, with k = 10 regressors:
# * main (results/fixed_effects.csv): every package, the four data
#   generating processes of fixed_effects/dgp.R (balanced, random,
#   round_robin, low_mobility), n = 1e3 to 1e6;
# * at scale (results/fixed_effects_scale.csv): lm_xtx() in memory,
#   lm_xtx() on DuckDB reading the Parquet file, and feols(); the balanced
#   and random processes, n = 1e7 and 5e7, no reference coefficients. Memory: the 5e7 dataset is about 5 GB in R
#   (as much Parquet on disk), and a fit in memory needs two to three times
#   that: allow 32 GB of RAM, or set FE_BENCH_SCALE_N="1e7".
#
# Data: each (process, n) dataset is simulated once (fixed_effects/simulate.R,
# with the data generating process of fixed_effects/dgp.R and the same seed
# for every cell) and saved as Parquet in a cache folder, FE_BENCH_CACHE
# (default tools::R_user_dir("xtx", "cache")/fixed_effects). The
# file name holds the process, n, the numbers of regressors and a hash of
# dgp.R and simulate.R: editing either makes new files (delete the folder
# to reclaim the space). One file holds the responses for every k (y_<k>),
# so each fit reads exactly the data it would have simulated itself.
#
# Cells: each (package, process, n, k) runs in a fresh R
# process (fixed_effects/cell.R), which reads the Parquet file (into a data
# frame, or as a DuckDB view for the DuckDB fits): one warmup call, then the
# median of three timed calls. Memory is the extra peak resident memory of
# the warmup call, so it counts allocations made in C, C++ and DuckDB. Once
# a package's fit takes more than `limit` seconds at some size, or fails
# after that long (dbreg's alternating projections can stop without
# converging), it is not run at larger sizes; a cell that runs longer than
# `kill` seconds is stopped, sooner for the cells in `kill_short` (known not
# to converge). The `kill_s` column gives each cell's limit.
#
# References (main grid only): after the timed cells, feols() at a tight
# tolerance (fixef.tol = 1e-11, fixef.iter = 1e5) runs once per (process, n,
# k) for which some fit finished, in its own process
# (fixed_effects/reference.R), stopped after FE_BENCH_REF_KILL seconds. It
# is cached next to the data (key: data file, fixed effects, k, fixest
# version, hash of dgp.R and reference.R); a reference that failed is cached
# too, one that was stopped is retried by a run with a longer
# FE_BENCH_REF_KILL. The `reference` column gives its status (ok, error,
# killed, or skipped: a smaller size was stopped); max_coef_diff compares
# each fit with it.
#
# BENCH_SIZE=quick (make bench-quick) runs a small preset instead: n = 1e3,
# 1e4, 1e5, lm_xtx and feols, main grid only, marked provisional;
# run.R then saves the results outside results/ (BENCH_RESULTS). Environment
# variables also shrink the grid; the results are then marked provisional:
#
#   FE_BENCH_N="1e3,1e4,1e5" taskset -c 2 Rscript website/benchmarks/run.R fixed_effects
#
# FE_BENCH_PACKAGES selects packages (comma separated, names as in `full`),
# FE_BENCH_DGPS the processes of the main grid, FE_BENCH_K the numbers of
# regressors,
# FE_BENCH_GRIDS the grids ("main,scale" by default), FE_BENCH_SCALE_N the
# at-scale sizes. FE_BENCH_LIMIT and FE_BENCH_KILL set `limit` and `kill`
# in seconds for the main grid (default 120 and 1800), FE_BENCH_KILL_SHORT
# the stop of the `kill_short` cells (default 600), FE_BENCH_SCALE_LIMIT and
# FE_BENCH_SCALE_KILL those of the at-scale grid (default 600 and 3600); all
# are saved in results/fixed_effects_machine.txt (and _scale_machine.txt).
processes <- c("balanced", "random", "round_robin", "low_mobility")
full <- list(
    n = c(1e3, 1e4, 1e5, 1e6),
    k = 10,
    package = c("lm_xtx", "lm_xtx (DuckDB)", "feols", "dbreg (DuckDB)"))
preset_quick <- list(n = c(1e3, 1e4, 1e5), k = 10, package = c("lm_xtx", "feols"))
scale <- list(
    dgp = c("balanced", "random"),
    n = c(1e7, 5e7),
    k = 10,
    package = c("lm_xtx", "lm_xtx (DuckDB)", "feols"))
# Cells known not to converge (package = processes): stopped after
# FE_BENCH_KILL_SHORT seconds
kill_short <- list("dbreg (DuckDB)" = "round_robin")

quick <- identical(Sys.getenv("BENCH_SIZE", "full"), "quick")
results <- Sys.getenv("BENCH_RESULTS", "results")
grid_env <- function(name, default, numeric = TRUE) {
    v <- Sys.getenv(name)
    if (v == "") return(default)
    v <- trimws(strsplit(v, ",")[[1]])
    if (numeric) as.numeric(v) else v
}
default <- if (quick) preset_quick else full
sizes <- grid_env("FE_BENCH_N", default$n)
ks <- grid_env("FE_BENCH_K", default$k)
packages <- grid_env("FE_BENCH_PACKAGES", default$package, numeric = FALSE)
dgps <- grid_env("FE_BENCH_DGPS", processes, numeric = FALSE)
grids <- grid_env("FE_BENCH_GRIDS", if (quick) "main" else c("main", "scale"), numeric = FALSE)
scale_sizes <- grid_env("FE_BENCH_SCALE_N", scale$n)
limit <- as.numeric(Sys.getenv("FE_BENCH_LIMIT", "120"))
kill <- as.numeric(Sys.getenv("FE_BENCH_KILL", "1800"))
kill_short_s <- as.numeric(Sys.getenv("FE_BENCH_KILL_SHORT", "600"))
scale_limit <- as.numeric(Sys.getenv("FE_BENCH_SCALE_LIMIT", "600"))
scale_kill <- as.numeric(Sys.getenv("FE_BENCH_SCALE_KILL", "3600"))
ref_kill <- as.numeric(Sys.getenv("FE_BENCH_REF_KILL", "1800"))
provisional <- quick || !(setequal(sizes, full$n) && setequal(ks, full$k) && setequal(packages, full$package) &&
                          setequal(dgps, processes))
scale_provisional <- quick || !setequal(scale_sizes, scale$n)
stopifnot(all(grids %in% c("main", "scale")), all(packages %in% full$package), all(dgps %in% processes),
          !anyNA(c(sizes, ks, scale_sizes)))

cache <- Sys.getenv("FE_BENCH_CACHE", file.path(tools::R_user_dir("xtx", "cache"), "fixed_effects"))
dir.create(cache, showWarnings = FALSE, recursive = TRUE)
cache <- normalizePath(cache)
# A short hash of some files: part of the cache keys
hash <- function(files) {
    tf <- tempfile()
    writeLines(unname(tools::md5sum(files)), tf)
    substr(unname(tools::md5sum(tf)), 1, 8)
}
data_hash <- hash(c("fixed_effects/dgp.R", "fixed_effects/simulate.R"))
ref_hash <- hash(c("fixed_effects/dgp.R", "fixed_effects/reference.R"))
fixest_version <- tryCatch(as.character(packageVersion("fixest")), error = function(e) "none")

rscript <- file.path(R.home("bin"), "Rscript")
one_thread <- c("OMP_NUM_THREADS=1", "OPENBLAS_NUM_THREADS=1")
tmp <- tempfile(fileext = ".csv")

# The Parquet file of one process and size, with responses for `with_k`:
# simulated in its own process the first time it is needed, NA if that
# fails (the cells that need it are then errors).
data_file <- function(dgp, n, with_k) {
    path <- file.path(cache, sprintf("data_%s_n%.0f_k%s_%s.parquet", dgp, n,
                                     paste(sort(with_k), collapse = "-"), data_hash))
    if (!file.exists(path)) {
        message(sprintf("simulating %s, n = %g", dgp, n))
        status <- system2(rscript, c("fixed_effects/simulate.R", dgp, n, paste(with_k, collapse = ","),
                                     shQuote(path)), env = one_thread)
        if (status != 0 || !file.exists(path)) {
            message(sprintf("could not simulate %s, n = %g", dgp, n))
            return(NA_character_)
        }
    }
    path
}

# Time one grid; `cell_kill(package, dgp)` gives a cell's stop in seconds.
run_grid <- function(dgps, ks, packages, sizes, limit, cell_kill) {
    rows <- list()
    for (dgp in dgps) for (k in ks) for (package in packages) {
        too_slow <- FALSE
        cell_kill_s <- cell_kill(package, dgp)
        for (n in sort(sizes)) {
            cell <- data.frame(package = package, dgp = dgp, n = n, k = k,
                               seconds = NA_real_, reps = NA_integer_, mem_mb = NA_real_,
                               coef = NA_character_, status = NA_character_, kill_s = cell_kill_s)
            path <- if (too_slow) NA_character_ else data_file(dgp, n, ks)
            if (too_slow) {
                cell$status <- "skipped"
            } else if (is.na(path)) {
                # No data: an error, and the larger sizes are skipped
                cell$status <- "error"
                too_slow <- TRUE
            } else {
                message(sprintf("%s, %s, k = %d, n = %g", package, dgp, k, n))
                unlink(tmp)
                wall <- Sys.time()
                status <- system2(rscript,
                    c("fixed_effects/cell.R", shQuote(package), shQuote(path), dgp, k, limit, tmp),
                    env = one_thread, timeout = cell_kill_s)
                wall <- as.numeric(difftime(Sys.time(), wall, units = "secs"))
                killed <- identical(status, 124L)
                out <- if (file.exists(tmp)) read.csv(tmp, colClasses = c(coef = "character", fit_status = "character"))
                if (!is.null(out) && out$fit_status == "ok") {
                    cell[c("n", "seconds", "reps", "mem_mb", "coef")] <- out[c("n", "seconds", "reps", "mem_mb", "coef")]
                    cell$status <- if (out$reps == 0) "single" else "ok"
                    fit_seconds <- out$fit_seconds
                } else {
                    if (!is.null(out)) cell$n <- out$n
                    cell$status <- if (killed) "killed" else "error"
                    # The time to failure of the fit, or of the whole process if
                    # it died before saying so
                    fit_seconds <- if (!is.null(out)) out$fit_seconds else wall
                }
                # A fit that failed (e.g. did not converge) after more than `limit`
                # seconds would only fail more slowly at larger sizes.
                too_slow <- cell$status %in% c("single", "killed") || isTRUE(cell$seconds > limit) ||
                    (cell$status == "error" && isTRUE(fit_seconds > limit))
            }
            rows[[length(rows) + 1]] <- cell
        }
    }
    res <- do.call(rbind, rows)
    res$seconds <- signif(res$seconds, 4)
    res
}

# The reference coefficients of one (process, n, k), from the cache or from
# reference.R: list(status, coef, kill_s).
reference <- function(path, dgp, k) {
    if (is.na(path)) return(list(status = "error", coef = NA_character_, kill_s = ref_kill))
    key <- file.path(cache, sprintf("ref_%s_k%d_fixest%s_%s.csv",
                                    sub("\\.parquet$", "", basename(path)), k, fixest_version, ref_hash))
    if (file.exists(key)) {
        r <- read.csv(key, colClasses = c(status = "character", coef = "character"))
        if (r$status != "killed" || r$kill_s >= ref_kill) return(as.list(r))
    }
    unlink(tmp)
    status <- system2(rscript, c("fixed_effects/reference.R", shQuote(path), dgp, k, tmp),
                      env = one_thread, timeout = ref_kill)
    r <- if (identical(status, 124L)) {
        data.frame(status = "killed", coef = NA_character_, kill_s = ref_kill)
    } else if (file.exists(tmp)) {
        cbind(read.csv(tmp, colClasses = c(status = "character", coef = "character"))[c("status", "coef")],
              kill_s = ref_kill)
    }
    # A process that died without a result (out of memory, say) is not cached
    if (is.null(r)) return(list(status = "error", coef = NA_character_, kill_s = ref_kill))
    write.csv(r, key, row.names = FALSE)
    as.list(r)
}

# Largest absolute difference between coefficients, stored as "b1;b2;..."
split_num <- function(s) if (is.na(s)) NA_real_ else as.numeric(strsplit(s, ";")[[1]])
coef_diff <- function(b, r) {
    if (is.na(b) || is.na(r)) return(NA_real_)
    signif(max(abs(split_num(b) - split_num(r))), 3)
}

deps <- list("lm_xtx" = c("xtx", "xtx.fit"), "lm_xtx (DuckDB)" = c("xtx", "xtx.fit", "duckdb"),
             "feols" = "fixest", "dbreg (DuckDB)" = c("dbreg", "duckdb"))
cpu <- tryCatch(trimws(sub(".*:", "", grep("^model name", readLines("/proc/cpuinfo"), value = TRUE)[1])),
                error = function(e) "unknown")
# The machine, the package versions and the thresholds, next to the results.
# Only the packages that ran are looked up (and duckdb, which reads the
# data), and a failed lookup never stops the results from being saved.
write_machine <- function(file, packages, lines) {
    used <- unique(c(unlist(deps[packages]), "duckdb"))
    versions <- vapply(used, function(p) tryCatch(as.character(packageVersion(p)),
                                                  error = function(e) "not installed"), "")
    tryCatch(writeLines(c(
        paste("CPU:", cpu, "(one core, one thread)"),
        R.version.string,
        paste("Packages:", paste(names(versions), versions, collapse = ", ")),
        lines),
        file.path(results, file)),
        error = function(e) warning("could not write ", file, ": ", conditionMessage(e)))
}
run_note <- function(provisional) {
    paste("Run:", format(Sys.Date()),
          if (quick) "(provisional: quick preset, BENCH_SIZE=quick)" else if (provisional) "(provisional: reduced grid)")
}

main <- NULL
if ("main" %in% grids) {
    res <- run_grid(dgps, ks, packages, sizes, limit,
                    function(package, dgp) if (dgp %in% kill_short[[package]]) min(kill, kill_short_s) else kill)

    # References, outside the timed cells, for every (process, n, k) with a
    # finished fit; once one is stopped, the larger sizes are skipped
    res$ref_coef <- NA_character_
    res$reference <- NA_character_
    res$ref_kill_s <- NA_real_
    for (dgp in dgps) for (k in ks) {
        stopped <- FALSE
        for (n in sort(sizes)) {
            i <- which(res$dgp == dgp & res$k == k & res$n == n)
            if (!any(res$status[i] %in% c("ok", "single"))) next
            if (stopped) {
                res$reference[i] <- "skipped"
                next
            }
            message(sprintf("reference: %s, k = %d, n = %g", dgp, k, n))
            r <- reference(data_file(dgp, n, ks), dgp, k)
            res$ref_coef[i] <- r$coef
            res$reference[i] <- r$status
            res$ref_kill_s[i] <- r$kill_s
            stopped <- r$status == "killed"
        }
    }
    res$max_coef_diff <- mapply(coef_diff, res$coef, res$ref_coef, USE.NAMES = FALSE)
    res$provisional <- provisional

    write_machine("fixed_effects_machine.txt", packages, c(
        sprintf("Limits: larger sizes skipped after %g s (FE_BENCH_LIMIT); cells stopped after %g s (FE_BENCH_KILL)",
                limit, kill),
        sprintf("Stopped sooner: %s after %g s (FE_BENCH_KILL_SHORT); references stopped after %g s (FE_BENCH_REF_KILL)",
                paste(names(kill_short), "on", unlist(kill_short), collapse = ", "),
                kill_short_s, ref_kill),
        run_note(provisional)))
    main <- res[, c("package", "dgp", "n", "k", "seconds", "reps", "mem_mb",
                    "max_coef_diff", "status", "provisional", "reference", "kill_s", "ref_kill_s")]
}

scale_packages <- intersect(scale$package, packages)
if ("scale" %in% grids && length(scale_sizes) > 0 && length(scale_packages) > 0) {
    sc <- run_grid(scale$dgp, scale$k, scale_packages, scale_sizes, scale_limit,
                   function(package, dgp) scale_kill)
    sc$provisional <- scale_provisional
    write.csv(sc[, c("package", "dgp", "n", "k", "seconds", "reps", "mem_mb", "status",
                     "provisional", "kill_s")],
              file.path(results, "fixed_effects_scale.csv"), row.names = FALSE)
    write_machine("fixed_effects_scale_machine.txt", scale_packages, c(
        sprintf("Limits: larger sizes skipped after %g s (FE_BENCH_SCALE_LIMIT); cells stopped after %g s (FE_BENCH_SCALE_KILL)",
                scale_limit, scale_kill),
        run_note(scale_provisional)))
}

# The main grid, saved by run.R as fixed_effects.csv (nothing if it did not run)
main
