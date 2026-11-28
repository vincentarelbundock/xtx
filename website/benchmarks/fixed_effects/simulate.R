# Simulate one dataset of the fixed-effects benchmark (see ../fixed_effects.R)
# and save it as Parquet, in a fresh R process:
#
#   Rscript simulate.R <dgp> <n> <ks> <out.parquet>
#
# <ks> lists the numbers of regressors, comma separated: the file has
# x1..x_max(ks) and one response y_<k> per k, each the response of
# dgp(n, <dgp>, k) (see dgp.R), so that every fit reads exactly the data it
# would have simulated itself. The file is written under a temporary name,
# then renamed.
args <- commandArgs(trailingOnly = TRUE)
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(here, "dgp.R"))
type <- args[1]
n <- as.numeric(args[2])
ks <- sort(unique(as.integer(strsplit(args[3], ",")[[1]])))
out <- args[4]

kmax <- max(ks)
d <- dgp(n, type, kmax)
d[[paste0("y_", kmax)]] <- d$y
for (k in setdiff(ks, kmax)) {
    dk <- dgp(n, type, k)
    # The regressors and the ids do not depend on k (only the response does)
    for (v in c(paste0("x", seq_len(k)), fe_ids(type))) stopifnot(identical(dk[[v]], d[[v]]))
    d[[paste0("y_", k)]] <- dk$y
    rm(dk)
}
d <- d[c(paste0("y_", ks), paste0("x", seq_len(kmax)), fe_ids(type))]

con <- DBI::dbConnect(duckdb::duckdb())
duckdb::duckdb_register(con, "d", d)
part <- paste0(out, ".part")
unlink(part)
invisible(DBI::dbExecute(con, sprintf("COPY d TO '%s' (FORMAT parquet)", gsub("'", "''", part))))
DBI::dbDisconnect(con, shutdown = TRUE)
stopifnot(file.rename(part, out))
