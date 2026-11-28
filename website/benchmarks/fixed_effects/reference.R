# Reference coefficients of the fixed-effects benchmark (see
# ../fixed_effects.R): feols() at a tight tolerance, in a fresh R process,
# outside the timed cells:
#
#   Rscript reference.R <data.parquet> <dgp> <k> <out.csv>
#
# Writes one row to <out.csv>: status "ok" with the coefficients, or "error"
# if the fit failed (the driver records a process it stops as "killed").
args <- commandArgs(trailingOnly = TRUE)
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(here, "dgp.R"))
path <- args[1]
dgp_type <- args[2]
k <- as.integer(args[3])
out <- args[4]

d <- fe_read(path, dgp_type, k)
f <- fe_formula(dgp_type, k)
xs <- paste0("x", seq_len(k))
suppressPackageStartupMessages(library(fixest))
setFixest_nthreads(1)
data.table::setDTthreads(1)

start <- Sys.time()
ref <- tryCatch(
    coef(feols(f, d, vcov = "iid", notes = FALSE, fixef.tol = 1e-11, fixef.iter = 1e5))[xs],
    error = function(e) {
        message("reference fit failed: ", conditionMessage(e))
        NULL
    })
ok <- !is.null(ref) && length(ref) == length(xs) && all(is.finite(ref))
row <- data.frame(status = if (ok) "ok" else "error",
                  coef = if (ok) paste(sprintf("%.17g", ref), collapse = ";") else NA_character_,
                  seconds = round(as.numeric(difftime(Sys.time(), start, units = "secs")), 2))
part <- paste0(out, ".part")
write.csv(row, part, row.names = FALSE)
invisible(file.rename(part, out))
