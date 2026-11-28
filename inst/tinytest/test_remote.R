## Remote designs on DuckDB: results agree with the design in memory.
if (!requireNamespace("duckdb", quietly = TRUE) || !requireNamespace("DBI", quietly = TRUE))
    exit_file("duckdb is not installed")
source("helpers.R")
library(xtx)

con <- test_con()
scratch <- tempfile("xtx_remote")
dir.create(scratch)

## data as the database returns it (ordered factors become plain factors,
## row names are dropped), with an id to match rows of predictions
remote_data <- function(n, na = TRUE) {
    d <- fuzz_data(n, special = FALSE, na = na)
    d$m <- NULL
    rownames(d) <- NULL
    d$id <- seq_len(n)
    d$w <- runif(n)
    d$g <- sample(c(3, 1, 2, NA), n, TRUE, prob = c(0.3, 0.3, 0.3, 0.1))
    d
}

same <- function(a, b) isTRUE(all.equal(unname(a), unname(b), tolerance = 1e-7))

## the rows of `dd` that the design in memory keeps
kept_rows <- function(fo, dd, spec = NULL) {
    mf <- if (is.null(spec)) {
        design_frame(fo, dd)
    } else {
        design_frame(spec$terms, dd, xlev = spec$xlevels, na.action = stats::na.omit)
    }
    as.integer(rownames(mf))
}

## compare a remote design with the design in memory, on every operation
check_remote <- function(Dr, M, dd, rows, info) {
    expect_equal(dim(Dr), dim(M), info = info)
    expect_identical(colnames(Dr), colnames(M), info = info)
    expect_true(same(crossprod(Dr), crossprod(M)), info = paste(info, "X'X"))
    expect_true(same(crossprod(Dr, "w"), crossprod(M, dd$w[rows])), info = paste(info, "X'w"))
    w <- dd$w[rows]
    expect_true(same(crossprod(Dr, "w" * Dr), crossprod(M, w * M)), info = paste(info, "X'WX"))
    expect_true(same(crossprod(sqrt(2) * Dr), crossprod(sqrt(2) * M)), info = paste(info, "(sX)'(sX)"))
    expect_true(same(crossprod("w" * Dr, "x1"), crossprod(M, w * dd$x1[rows])), info = paste(info, "X'Wx"))
    expect_true(same(colSums(Dr), colSums(M)), info = paste(info, "colSums"))
    g <- dd$g[rows]
    wg <- if (anyNA(g)) "missing values" else NA
    a <- capture(function() colSums(Dr, group = "g"))
    b <- capture(function() colSums(M, group = g))
    expect_true(same(a$value, b$value), info = paste(info, "colSums(group)"))
    expect_identical(rownames(a$value), rownames(b$value), info = paste(info, "group names"))
    expect_identical(a$warnings, b$warnings, info = paste(info, "group warning"))
    a <- suppressWarnings(colMeans("w" * Dr, group = "g"))
    b <- suppressWarnings(colMeans(dd$w[rows] * M, group = g))
    expect_true(same(a, b), info = paste(info, "colMeans(group, weights)"))
    beta <- rnorm(ncol(M))
    p <- as.data.frame(Dr %*% beta, keep = "id")
    p <- p[order(p$id), ]
    expect_identical(p$id, dd$id[rows], info = paste(info, "prediction rows"))
    expect_true(same(p$.fitted, drop(M %*% beta)), info = paste(info, "predictions"))
}

set.seed(20261125)
dd <- remote_data(400)
DBI::dbWriteTable(con, "d", dd)
dd <- DBI::dbReadTable(con, "d")
dd <- dd[order(dd$id), ]
rownames(dd) <- NULL

## random formulas: computed in SQL, or with variables evaluated in R
checked <- c(sql = 0, evaluated = 0)
for (i in seq_len(if (full_tests) 40 else 30)) {
    fo <- fuzz_formula()
    ## (a remote design needs a plan, which refuses the response on the right)
    if ("m" %in% all.vars(fo) || "y" %in% all.vars(fo[[3L]])) next
    ctr <- fuzz_contrasts(dd)
    M <- capture(function() design_matrix(fo, dd, contrasts.arg = ctr))
    for (engine in "auto") {
        info <- paste(deparse(fo), engine)
        Dr <- capture(function() suppressMessages(
            design_matrix(fo, con, table = "d", contrasts.arg = ctr, engine = engine, chunk_size = 37)))
        if (inherits(M$value, "captured_error")) {
            expect_inherits(Dr$value, "captured_error", info = info)
            next
        }
        if (inherits(Dr$value, "captured_error")) {
            expect_true(FALSE, info = paste(info, Dr$value$msg))
            next
        }
        Dr <- Dr$value
        how <- if (is.null(Dr$evaluated)) "sql" else "evaluated"
        checked[how] <- checked[how] + 1
        check_remote(Dr, M$value, dd, kept_rows(fo, dd), info)
        expect_true(same(crossprod(Dr, "y"), crossprod(M$value, dd$y[kept_rows(fo, dd)])),
                    info = paste(info, "X'y"))
    }
}
expect_true(all(checked > 5))

## fit in R, score in the database: a specification from part of the data
for (i in seq_len(if (full_tests) 25 else 10)) {
    fo <- fuzz_formula(response = FALSE)
    if ("m" %in% all.vars(fo)) next
    train <- dd[1:120, ]
    spec <- capture(function() design_matrix(fo, train, contrasts.arg = fuzz_contrasts(train)))
    if (inherits(spec$value, "captured_error")) next
    spec <- spec$value
    plan <- design_matrix(spec, newdata = NULL)
    ## rows with missing values kept (na.pass, the default) or dropped
    for (na in c("na.pass", "na.omit")) {
        M <- capture(function() design_matrix(spec, dd, na.action = get(na)))
        rows <- if (na == "na.pass") seq_len(nrow(dd)) else kept_rows(NULL, dd, plan)
        for (engine in "auto") {
            info <- paste("spec", deparse(fo), engine, na)
            Dr <- capture(function() suppressMessages(
                design_matrix(spec, con, table = "d", engine = engine, chunk_size = 50, na.action = get(na))))
            if (inherits(M$value, "captured_error")) {
                ## new levels: the same error as in memory
                expect_inherits(Dr$value, "captured_error", info = info)
                if (inherits(Dr$value, "captured_error")) expect_identical(Dr$value$msg, M$value$msg, info = info)
                next
            }
            if (inherits(Dr$value, "captured_error")) {
                expect_true(FALSE, info = paste(info, Dr$value$msg))
                next
            }
            check_remote(Dr$value, M$value, dd, rows, info)
        }
    }
}

## na.pass explicitly: NULL predictions, NA in the columns of terms with a
## missing value, the other columns computed over all rows
train <- data.frame(x = rnorm(20), f = rep(c("a", "b"), 10), z = runif(20))
spec <- design_matrix(~ x + f + z, train)
nd <- data.frame(x = c(1, NA, 3, NaN, 5), f = c("a", "b", NA, "a", "b"), z = 1:5, g = c(1, 1, 2, 2, 2))
DBI::dbWriteTable(con, "nd", nd)
M <- design_matrix(spec, nd)
for (engine in "auto") {
    Dr <- suppressMessages(design_matrix(spec, con, table = "nd", engine = engine))
    expect_equal(nrow(Dr), 5)
    expect_true(same(crossprod(Dr), crossprod(M)), info = engine)
    expect_identical(is.na(crossprod(Dr)), is.na(crossprod(M)), info = engine)
    expect_true(same(colSums(Dr, group = "g"), colSums(M, group = nd$g)), info = engine)
    expect_identical(is.na(colSums(Dr, group = "g")), is.na(colSums(M, group = nd$g)), info = engine)
    p <- as.data.frame(Dr %*% seq_len(4), keep = "z")
    p <- p[order(p$z), ]
    expect_identical(is.na(p$.fitted), c(FALSE, TRUE, TRUE, TRUE, FALSE), info = engine)
    expect_true(same(p$.fitted, drop(M %*% seq_len(4))), info = engine)
}

## a column of another type than in training: the error of predict.lm(),
## before anything is built
train <- data.frame(x = rnorm(20), f = rep(c("a", "b"), 10))
spec <- design_matrix(~ x + f, train)
DBI::dbWriteTable(con, "badtype", data.frame(x = as.character(1:200), f = "a"))
bad <- DBI::dbReadTable(con, "badtype")
msg <- tryCatch(design_matrix(spec, bad), error = conditionMessage)
expect_true(grepl("fitted with type", msg))
for (engine in c("auto", "sql")) {
    expect_error(design_matrix(spec, con, table = "badtype", engine = engine), msg, fixed = TRUE)
}
## a factor stored as text is fine (model.frame() converts it)
spec2 <- design_matrix(~ x + f, transform(train, f = factor(f)))
DBI::dbWriteTable(con, "texttype", data.frame(x = 1:3, f = c("a", "b", "a")))
expect_inherits(design_matrix(spec2, con, table = "texttype"), "xtx_remote")
## a variable computed in R: checked on each batch
spec3 <- design_matrix(~ I(x) + f, train)
msg3 <- tryCatch(design_matrix(spec3, bad), error = conditionMessage)
expect_true(grepl("fitted with type", msg3))
expect_error(suppressMessages(crossprod(design_matrix(spec3, con, table = "badtype"))), msg3, fixed = TRUE)

## new levels, explicitly
train <- data.frame(x = rnorm(10), ch = rep(c("u", "v"), 5))
spec <- design_matrix(~ x + ch, train)
DBI::dbWriteTable(con, "newlev", data.frame(x = 1:4, ch = c("u", "w", "a", NA)))
msg <- tryCatch(design_matrix(spec, DBI::dbReadTable(con, "newlev"), na.action = stats::na.omit),
                error = conditionMessage)
expect_error(design_matrix(spec, con, table = "newlev"), msg, fixed = TRUE)

## NaN and NA rows, transformations giving NaN
na <- data.frame(y = c(1, 2, NA, 4, 5, 6, 7, 8), x = c(1, NaN, 3, 4, -1, 6, 7, 8),
                 f = c("a", "b", "a", NA, "b", "a", "b", "a"), w = 1, g = 1, id = 1:8)
DBI::dbWriteTable(con, "na", na)
fo <- y ~ log(x) * f
M <- suppressWarnings(design_matrix(fo, na))
Dr <- design_matrix(fo, con, table = "na")
expect_identical(Dr$engine, "sql")
expect_equal(nrow(Dr), nrow(M))
expect_true(same(crossprod(Dr), crossprod(M)))
expect_true(same(crossprod(Dr, "y"), crossprod(M, na$y[c(1, 6, 7, 8)])))
expect_error(crossprod(Dr, weights = "w"), "only `x` and `y`")
expect_true(same(crossprod(Dr, "w" * Dr), crossprod(M)))
expect_error(Dr + Dr, "not supported")
expect_error(Dr * Dr, "two remote designs")
expect_error(1:3 * Dr, "only be scaled")
expect_error(crossprod(Dr, suppressMessages(design_matrix(fo, con, table = "na"))), "same design")
## a column with missing values in the rows used gives NA, as in memory
DBI::dbWriteTable(con, "na2", transform(na, w = c(1, 1, 1, 1, 1, NA, 1, 1)))
Dr2 <- design_matrix(fo, con, table = "na2")
expect_true(all(is.na(crossprod(Dr2, "w"))))

## log(x), log(x, base), log2, log10, log1p in formula terms (y12r): the
## values of model.matrix(), also at 0 (-Inf), below (NaN), for NA, NaN,
## Inf, -Inf, for any base, and near 0 for log1p
lx <- c(NA, NaN, Inf, -Inf, 0, -0.5, -1, -2, 1e-10, -1e-12, 3e-17, 1, 2.5, 1e300)
DBI::dbWriteTable(con, "lg", data.frame(y = seq_along(lx), x = lx))
lgf <- list(y ~ log(x), y ~ log(x, 2), y ~ log(x, 10), y ~ log(x, 3), y ~ log(base = 0.5, x),
            y ~ log(x, 1), y ~ log(x, 0), y ~ log(x, Inf), y ~ log(x, 1e-3), y ~ log2(x), y ~ log10(x),
            y ~ log1p(x), y ~ log1p(2 * x) + log(x + 1, 2), y ~ log(abs(x), 0.25) + x)
for (fo in lgf) {
    M <- suppressWarnings(design_matrix(fo, data.frame(y = seq_along(lx), x = lx)))
    Dr <- design_matrix(fo, con, table = "lg")
    info <- deparse(fo)
    expect_identical(Dr$engine, "sql", info = info)
    expect_true(isTRUE(all.equal(unname(as.matrix(Dr)), unname(as.matrix(M)), tolerance = 1e-15)), info = info)
}
## a base that is not a constant is evaluated in R
expect_message(design_matrix(y ~ log(x, y), con, table = "lg"), "Evaluated in R", fixed = TRUE)

expect_message(design_matrix(y ~ poly(z, 2) + f1, con, table = "d"),
               "Evaluated in R: term `poly(z, 2)` cannot be translated", fixed = TRUE)
expect_error(design_matrix(y ~ poly(z, 2) + f1, con, table = "d", engine = "sql"),
             "SQL engine cannot compute this design", fixed = TRUE)
expect_silent(D <- design_matrix(y ~ x1 * f1 + I(x2^2) + log(z, 10), con, table = "d"))
expect_identical(D$engine, "sql")
expect_stdout(print(summary(D, sql = TRUE)), "GROUP BY")
expect_silent(ex <- summary(D))
expect_inherits(ex, "summary.xtx_remote")
expect_identical(ex$sql, rep(TRUE, nrow(ex)))
expect_stdout(print(D), "all terms pushed down")
expect_identical(xtx:::xtx_fmt_count(1e6), "1,000,000")
Dc <- suppressMessages(design_matrix(y ~ poly(z, 2) + f1, con, table = "d"))
out <- capture.output(print(Dc))
expect_identical(sum(grepl("cannot be translated", out, fixed = TRUE)), 0L)
expect_identical(sum(grepl("^engine: sql, on columns computed in R.*poly\\(z, 2\\)", out)), 1L)

## sql= source, as.matrix, reuse of a remote design
D <- design_matrix(y ~ x1 + f1, con, sql = "SELECT * FROM d WHERE id <= 200")
M <- design_matrix(y ~ x1 + f1, dd[dd$id <= 200, ])
expect_true(same(crossprod(D), crossprod(M)))
expect_equal(unname(as.matrix(D)), unname(as.matrix(M)))
op <- options(xtx.remote_max_rows = 10)
expect_error(as.matrix(D), "force = TRUE")
expect_equal(unname(as.matrix(D, force = TRUE)), unname(as.matrix(M)))
options(op)
expect_equal(as.matrix(design_matrix(D, dd[1:5, ])), as.matrix(design_matrix(M, dd[1:5, ])))
D2 <- design_matrix(D, con, table = "d", na.action = stats::na.omit)
expect_true(same(crossprod(D2), crossprod(design_matrix(M, dd, na.action = stats::na.omit))))

## predictions: a remote column (design_column), stored by col_compute(); scaled
## predictions; scaled as.matrix()
beta <- seq_len(ncol(D))
col <- D %*% beta
expect_inherits(col, "xtx_column")
expect_equal(ncol(col), 1)
expect_equal(length(col), nrow(D))
expect_equal(dim(col), c(nrow(D), 1))
p <- as.data.frame(col_compute(col), keep = "id")
p <- p[order(p$id), ]
expect_identical(names(p), c("id", ".fitted"))
expect_true(same(p$.fitted, drop(M %*% beta)))
ps <- as.data.frame((DBI::SQL("2 * x1") * D / "w") %*% beta, keep = c("id", "x1", "w"))
ps <- ps[order(ps$id), ]
expect_true(same(ps$.fitted, drop(M %*% beta) * 2 * ps$x1 / ps$w))
expect_equal(unname(as.matrix(-D)), unname(as.matrix(-M)))

## several prediction columns; a fitted model scored with `newdata`
fit <- lm(y ~ x1 + f1 + splines::ns(z, 2), dd)
B <- cbind(a = coef(fit), b = -coef(fit))
ref <- as.matrix(design_matrix(fit, newdata = dd)) %*% B
for (engine in "auto") {
    Df <- suppressMessages(design_matrix(fit, newdata = con, table = "d", engine = engine))
    p <- as.data.frame(Df %*% B, keep = "id")
    p <- p[order(p$id), ]
    expect_true(all(c("id", "a", "b") %in% names(p)))
    expect_true(same(as.matrix(p[c("a", "b")]), ref), info = engine)
}

## Parquet files, with a private connection
pq <- file.path(scratch, "d.parquet")
DBI::dbExecute(con, sprintf("COPY (SELECT * FROM d) TO '%s' (FORMAT PARQUET)", pq))
dp <- DBI::dbGetQuery(con, sprintf("SELECT * FROM read_parquet('%s') ORDER BY id", pq))
for (fo in list(y ~ x1 * f2 + ch + lg, y ~ splines::ns(z, 3) + f1)) {
    M <- design_matrix(fo, dp)
    for (engine in "auto") {
        Dp <- suppressMessages(design_matrix(fo, pq, engine = engine, chunk_size = 64))
        check_remote(Dp, M, dp, kept_rows(fo, dp), paste("parquet", deparse(fo), engine))
    }
}

## SQL expressions as weights and groups; sandwich meats (with the
## residuals of a remote fit: tests of xtx.fit)
set.seed(11)
ds <- data.frame(y = rnorm(300), x = rnorm(300), f = sample(letters[1:4], 300, TRUE),
                 cl = sample(1:12, 300, TRUE))
ds$y[5] <- NA
DBI::dbWriteTable(con, "ds", ds)
dk <- ds[!is.na(ds$y), ]
X <- model.matrix(y ~ x * f, dk)
for (engine in "sql") {
    D <- design_matrix(y ~ x * f, con, table = "ds", engine = engine, chunk_size = 64)
    sw <- DBI::SQL("1 + abs(x)")
    expect_true(same(crossprod(D, sw * D), crossprod(X, (1 + abs(dk$x)) * X)), info = engine)
    expect_true(same(crossprod(sw * D, "y"), crossprod(X, (1 + abs(dk$x)) * dk$y)), info = engine)
    expect_true(same(crossprod(D, DBI::SQL("y * x")), crossprod(X, dk$y * dk$x)), info = engine)
    expect_true(same(colSums(D, group = DBI::SQL("cl % 3")), rowsum(X, dk$cl %% 3)), info = engine)
    expect_true(same(colMeans(DBI::SQL("2 * x") * D, group = "cl"),
                     rowsum(2 * dk$x * X, dk$cl) / as.vector(table(dk$cl))), info = engine)
}
expect_error(DBI::SQL("nope + 1") * design_matrix(y ~ x, con, table = "ds"), "invalid SQL expression")

## errors
expect_error(design_matrix(y ~ x1, con), "exactly one of")
expect_error(design_matrix(y ~ nope, con, table = "d"))
expect_error(crossprod(design_matrix(y ~ x1, con, table = "d"), "nope"), "no column")

test_con_done(con)
unlink(scratch, recursive = TRUE)

## factor(x), as.factor(x), as.character(x) of a text, ENUM or integer
## column are factor variables on that column, computed in the database
## (no evaluation in R), with the levels R gives them: numeric order for
## factor() of an integer, the order of the type for an ENUM, sorted as
## text for as.character()
if (requireNamespace("duckdb", quietly = TRUE) && requireNamespace("DBI", quietly = TRUE)) {
    con <- DBI::dbConnect(duckdb::duckdb())
    set.seed(31)
    n <- 400
    df <- data.frame(x = rnorm(n), xi = sample(c(2L, 9L, 10L, 11L, -3L), n, TRUE), ch = sample(c("b", "a", "c"), n, TRUE),
                     en = factor(sample(c("z", "m", "a"), n, TRUE), levels = c("z", "m", "a")), lg = sample(c(TRUE, FALSE), n, TRUE))
    df$xi[7] <- NA
    DBI::dbWriteTable(con, "df", df)
    for (fo in list(~ factor(xi), ~ x * as.factor(xi), ~ as.character(xi), ~ factor(ch) + x, ~ as.character(ch),
                    ~ factor(en), ~ as.character(en) + factor(xi))) {
        info <- deparse(fo)
        msgs <- character()
        D <- withCallingHandlers(design_matrix(fo, con, table = "df"),
                                 message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
        M <- model.matrix(fo, df)
        expect_false(any(grepl("Evaluated in R", msgs)), info = info)
        expect_identical(colnames(D), colnames(M), info = info)
        expect_equal(unname(as.matrix(crossprod(D))), unname(crossprod(M)), tolerance = 1e-12, info = info)
        expect_identical(D$plan$xlevels, stats::.getXlevels(stats::terms(fo), stats::model.frame(fo, df)), info = info)
    }
    ## other columns and other calls are evaluated in R, as before
    for (fo in list(~ factor(lg), ~ factor(round(x)), ~ factor(xi, levels = c(2, 9)))) {
        msgs <- character()
        D <- withCallingHandlers(design_matrix(fo, con, table = "df"),
                                 message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
        expect_true(any(grepl("Evaluated in R", msgs)), info = deparse(fo))
        expect_identical(colnames(D), colnames(model.matrix(fo, df)), info = deparse(fo))
    }
    ## the specification of such a design applies to new remote data
    D <- design_matrix(~ factor(xi) + x, con, table = "df")
    D2 <- design_matrix(D$plan, con, table = "df", na.action = stats::na.omit)
    expect_equal(unname(as.matrix(crossprod(D2))), unname(as.matrix(crossprod(D))), tolerance = 1e-12)
    DBI::dbDisconnect(con)
}
