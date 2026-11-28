# One cell of the fixed-effects benchmark (see ../fixed_effects.R), in a fresh
# R process:
#
#   Rscript cell.R <package> <data.parquet> <dgp> <k> <limit> <out.csv>
#
# Reads the data that simulate.R saved (dgp.R): into a data frame for the
# in-memory packages, or as a DuckDB view that queries the Parquet file for
# the DuckDB ones. Then fits the model once (warmup), then three timed times,
# and writes one row to <out.csv>: the median time in seconds, the extra peak
# resident memory of the warmup call in MB, and the coefficients. When the
# warmup takes more than <limit> seconds, its time is reported alone (reps =
# 0) instead of three more runs. If the warmup fails, the row says so
# (fit_status = "error") with the seconds it took to fail. The reference
# coefficients are computed separately (reference.R).
args <- commandArgs(trailingOnly = TRUE)
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(here, "dgp.R"))
package <- args[1]
path <- args[2]
dgp_type <- args[3]
k <- as.integer(args[4])
limit <- as.numeric(args[5])
out <- args[6]

xs <- paste0("x", seq_len(k))
f <- fe_formula(dgp_type, k)

# The data: DuckDB runs on one thread, as every package below
remote <- package %in% c("lm_xtx (DuckDB)", "dbreg (DuckDB)")
if (remote) {
    # A view, so that every query reads the Parquet file: the rows are never
    # loaded into the database or into R
    suppressMessages({
        con <- DBI::dbConnect(duckdb::duckdb())
        DBI::dbExecute(con, "SET threads = 1")
        DBI::dbExecute(con, paste("CREATE VIEW d AS", fe_sql(path, dgp_type, k)))
    })
    nrows <- DBI::dbGetQuery(con, "SELECT count(*) AS n FROM d")$n
} else {
    d <- fe_read(path, dgp_type, k)
    nrows <- nrow(d)
}

# Every package on one thread
suppressPackageStartupMessages(switch(package,
    "lm_xtx" = ,
    "lm_xtx (DuckDB)" = {
        library(xtx)
        library(xtx.fit)
        options(xtx.threads = 1)
    },
    "feols" = {
        library(fixest)
        setFixest_nthreads(1)
        data.table::setDTthreads(1)
    },
    "dbreg (DuckDB)" = {
        library(dbreg)
        # The default (100 iterations of alternating projections) stops with
        # an error on round robin.
        options(dbreg.ap_max_iter = 1e4)
    },
    stop("unknown package: ", package)))

fit <- switch(package,
    "lm_xtx" = function() lm_xtx(f, d),
    "lm_xtx (DuckDB)" = function() lm_xtx(f, con, table = "d"),
    "feols" = function() feols(f, d, vcov = "iid", notes = FALSE),
    "dbreg (DuckDB)" = function() dbreg(f, conn = con, table = "d", vcov = "iid"))

# Peak resident memory of this process (Linux): reset the high-water mark,
# then read it after the call.
status_mb <- function(field) {
    line <- grep(paste0("^", field, ":"), readLines("/proc/self/status"), value = TRUE)
    as.numeric(gsub("\\D", "", line)) / 1024
}

# One row of results, written to a temporary file then renamed, so that a
# process stopped while writing never leaves a truncated <out.csv>.
save_row <- function(...) {
    row <- data.frame(package = package, dgp = dgp_type, n = nrows, k = k, ...)
    part <- paste0(out, ".part")
    write.csv(row, part, row.names = FALSE)
    invisible(file.rename(part, out))
}

invisible(gc())
writeLines("5", "/proc/self/clear_refs")
rss0 <- status_mb("VmRSS")
start <- Sys.time()
t0 <- tryCatch(system.time(model <- fit())[["elapsed"]], error = function(e) {
    message("fit failed: ", conditionMessage(e))
    NULL
})
if (is.null(t0)) {
    save_row(seconds = NA_real_, reps = NA_integer_, mem_mb = NA_real_,
             coef = NA_character_, fit_status = "error",
             fit_seconds = as.numeric(difftime(Sys.time(), start, units = "secs")))
    quit(save = "no", status = 1)
}
mem_mb <- status_mb("VmHWM") - rss0
b <- coef(model)[xs]

if (t0 > limit) {
    seconds <- t0
    reps <- 0L
} else {
    rm(model)
    times <- replicate(3, system.time(fit())[["elapsed"]])
    seconds <- median(times)
    reps <- 3L
}

save_row(seconds = seconds, reps = reps, mem_mb = round(mem_mb, 1),
         coef = paste(sprintf("%.17g", b), collapse = ";"), fit_status = "ok", fit_seconds = t0)
