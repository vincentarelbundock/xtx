# The data of the fixed-effects benchmark (see ../fixed_effects.R), shared by
# simulate.R (which saves each dataset once, as Parquet), cell.R (timed fits)
# and reference.R (tight-tolerance coefficients).
#
# Four data generating processes, adapted from Kyle Butts' fixest_benchmarks
# (as used by pyfixest). Each is a panel of individuals over 10 years, with
# k standard normal regressors; in all but the first, every row also belongs
# to a group, with about 23 individuals per group and year:
#
# * "balanced": individuals and years only (states observed every year).
# * "random": each row's group is drawn at random (shoppers and the
#   products they buy).
# * "round_robin": groups are assigned in turn over the rows, so every
#   individual moves to the next group every year (students moving up to
#   the next teacher).
# * "low_mobility": individuals keep their group and move to another one,
#   drawn at random, in 8% of years (workers and their firms).
#
# The model fitted has individual and year effects, and group effects when
# the process has groups (fe_formula()).
fe_processes <- c("balanced", "random", "round_robin", "low_mobility")

dgp <- function(n, type, k, seed = 42) {
    type <- match.arg(type, fe_processes)
    set.seed(seed)
    nb_year <- 10
    nb_indiv <- round(n / nb_year)
    nb_group <- round(nb_indiv / 23)
    n <- nb_indiv * nb_year
    indiv_id <- rep(seq_len(nb_indiv), each = nb_year)
    year <- rep(seq_len(nb_year), nb_indiv)
    group_id <- switch(type,
        balanced = NULL,
        random = sample.int(nb_group, n, TRUE),
        round_robin = rep_len(seq_len(nb_group), n),
        low_mobility = {
            # a new spell in an individual's first year, and when they move
            spell <- cumsum(year == 1 | runif(n) < 0.08)
            sample.int(nb_group, max(spell), TRUE)[spell]
        })
    # x1..xk ~ N(0, 1), column by column (same draws as an n by k matrix)
    d <- data.frame(y = numeric(n))
    y <- numeric(n)
    for (j in seq_len(k)) {
        x <- rnorm(n)
        d[[paste0("x", j)]] <- x
        y <- y + x / j
    }
    y <- y + rnorm(nb_indiv)[indiv_id] + rnorm(nb_year)[year]
    if (!is.null(group_id)) y <- y + rnorm(nb_group)[group_id]
    d$y <- y + rnorm(n)
    d$indiv_id <- indiv_id
    d$year <- year
    d$group_id <- group_id
    d
}

# The fixed effects of a process: individuals and years, and its groups
fe_ids <- function(type) c("indiv_id", "year", if (type != "balanced") "group_id")

# A Parquet file of simulate.R holds the data of one process and size for
# several k: x1..x_max(k) and the ids are the same draws whatever k (they
# come first in the random stream), and the response for k regressors is
# column y_<k>. This query gives exactly dgp(n, type, k): y, x1..xk, then
# the ids, in that order.
fe_sql <- function(path, type, k) {
    sprintf("SELECT y_%d AS y, %s, %s FROM read_parquet('%s')",
            as.integer(k), paste0("x", seq_len(k), collapse = ", "),
            paste(fe_ids(type), collapse = ", "), gsub("'", "''", path))
}

# The data for k regressors as a data frame (in-memory fits)
fe_read <- function(path, type, k) {
    con <- DBI::dbConnect(duckdb::duckdb())
    on.exit(DBI::dbDisconnect(con, shutdown = TRUE))
    DBI::dbExecute(con, "SET threads = 1")
    DBI::dbGetQuery(con, fe_sql(path, type, k))
}

fe_formula <- function(type, k) {
    xs <- paste0("x", seq_len(k))
    f <- as.formula(paste("y ~", paste(xs, collapse = " + "), "|", paste(fe_ids(type), collapse = " + ")))
    environment(f) <- globalenv()
    f
}
