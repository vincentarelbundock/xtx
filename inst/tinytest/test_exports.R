## The exported layer for package developers: design_parts(), plan_of()
## and design_patterns(). Only exported functions are called here (and base
## R / stats), in memory and on a DuckDB table.
source("helpers.R")
library(xtx)

expect_true(all(c("design_parts", "design_patterns") %in% getNamespaceExports("xtx")))

same <- function(a, b, tol = 1e-9) isTRUE(all.equal(unname(a), unname(b), tolerance = tol))

## ---- design_patterns() ----

set.seed(17)
n <- 600
d <- data.frame(f = sample(c("a", "b", "c"), n, TRUE), x = sample(0:4, n, TRUE) / 2,
                g = sample(1:3, n, TRUE), w = sample(c(0.5, 1, 2), n, TRUE), xc = rnorm(n))
d$y <- d$x + match(d$f, c("a", "b", "c")) + rnorm(n)
fo <- y ~ f * x

## The patterns as plain matrices, sorted by the rows of the design and the
## keys (the order of the patterns is not part of the contract).
sorted <- function(p) {
    X <- as.matrix(p$x)
    X <- matrix(as.double(X), nrow(X))
    K <- if (!is.null(p$keys)) matrix(as.double(p$keys), nrow(p$keys))
    M <- cbind(X, K)
    o <- do.call(order, unname(as.data.frame(M)))
    list(x = X[o, , drop = FALSE], z = p$z[o], w = p$w[o], within = p$within, rows = p$rows,
         keys = if (!is.null(K)) K[o, , drop = FALSE])
}

## The reference in base R: group the rows of model.matrix() (and keys).
reference <- function(X, z, w = NULL, K = NULL) {
    if (is.null(w)) w <- rep(1, length(z))
    M <- cbind(unname(X), K)
    key <- do.call(paste, c(lapply(seq_len(ncol(M)), function(j) sprintf("%a", as.double(M[, j]))), sep = "\r"))
    idx <- match(key, unique(key))
    sw <- drop(rowsum(w, idx, reorder = TRUE))
    m <- drop(rowsum(w * z, idx, reorder = TRUE)) / sw
    pos <- w > 0
    within <- sum(w[pos] * (z[pos] - m[idx[pos]])^2)
    keep <- sw > 0
    U <- M[!duplicated(key), , drop = FALSE][keep, , drop = FALSE]
    p <- list(x = U[, seq_len(ncol(X)), drop = FALSE], z = m[keep], w = sw[keep], within = within,
              rows = length(z), keys = if (!is.null(K)) U[, -seq_len(ncol(X)), drop = FALSE])
    sorted(p)
}

same_patterns <- function(a, b, tol = 1e-9) {
    same(a$x, b$x, tol) && same(a$z, b$z, tol) && same(a$w, b$w, tol) &&
        same(a$within, b$within, tol) && same(a$rows, b$rows) &&
        (is.null(a$keys) && is.null(b$keys) || same(a$keys, b$keys, tol))
}

D <- design_matrix(fo, d)
X <- model.matrix(fo, d)

## the public list
p <- design_patterns(D, d$y, compress = TRUE)
expect_identical(names(p), c("x", "z", "w", "within", "rows", "keys"))
expect_inherits(p$x, "xtx_matrix")
expect_identical(colnames(as.matrix(p$x)), colnames(X))
expect_true(nrow(p$x) < n)
expect_true(is.double(p$w) && length(p$w) == nrow(p$x))
expect_null(p$keys)
expect_equal(sum(p$w), n)
expect_equal(p$rows, n)
expect_true(same_patterns(sorted(p), reference(X, d$y)))

## weighted, and with keys
pw <- design_patterns(D, d$y, d$w, compress = TRUE)
expect_true(same_patterns(sorted(pw), reference(X, d$y, d$w)))
pk <- design_patterns(D, d$y, d$w, keys = cbind(g = d$g), compress = TRUE)
expect_identical(colnames(pk$keys), "g")
expect_true(same_patterns(sorted(pk), reference(X, d$y, d$w, cbind(d$g))))
expect_true(nrow(pk$x) > nrow(pw$x))

## zero weights: the rows do not enter the statistics
w0 <- d$w
w0[d$f == "a" & d$x == 0] <- 0
p0 <- design_patterns(D, d$y, w0, compress = TRUE)
expect_true(all(p0$w > 0))
expect_true(same_patterns(sorted(p0), reference(X, d$y, w0)))

## lm on the patterns, with their weights, is lm on the rows
l <- lm(fo, d)
lp <- lm(p$z ~ as.matrix(p$x) - 1, weights = p$w)
expect_true(same(coef(lp), coef(l)))
expect_true(same(sum(weights(lp) * residuals(lp)^2) + p$within, sum(residuals(l)^2)))
lw <- lm(fo, d, weights = w)
lpw <- lm(pw$z ~ as.matrix(pw$x) - 1, weights = pw$w)
expect_true(same(coef(lpw), coef(lw)))
expect_true(same(sum(weights(lpw) * residuals(lpw)^2) + pw$within, sum(weights(lw) * residuals(lw)^2)))

## when compression does not apply: NULL
expect_null(design_patterns(D, d$y, compress = FALSE))
Dc <- design_matrix(y ~ xc, d)
expect_null(design_patterns(Dc, d$y, compress = TRUE))

## errors
expect_error(design_patterns(as.matrix(D), d$y), "design")
expect_error(design_patterns(D, d$y, compress = "yes"), "compress")
expect_error(design_patterns(D, d$y[-1]), "one value per row")
yna <- d$y
yna[3] <- NA
expect_error(design_patterns(D, yna, compress = TRUE), "missing or infinite")
wneg <- d$w
wneg[2] <- -1
expect_error(design_patterns(D, d$y, wneg, compress = TRUE), "non-negative")

## the same patterns and statistics on a DuckDB table
if (requireNamespace("duckdb", quietly = TRUE) && requireNamespace("DBI", quietly = TRUE)) {
    con <- test_con()
    DBI::dbWriteTable(con, "d", d)
    d$yna <- yna
    d$wneg <- wneg
    d$w0 <- w0
    DBI::dbWriteTable(con, "d2", d)
    R <- suppressMessages(design_matrix(fo, con, table = "d"))
    R2 <- suppressMessages(design_matrix(fo, con, table = "d2"))

    pr <- design_patterns(R, "y", compress = TRUE)
    expect_identical(names(pr), names(p))
    expect_inherits(pr$x, "xtx_matrix")
    expect_identical(colnames(as.matrix(pr$x)), colnames(X))
    expect_true(same_patterns(sorted(pr), sorted(p)))

    prw <- design_patterns(R, "y", "w", compress = TRUE)
    expect_true(same_patterns(sorted(prw), sorted(pw)))
    prk <- design_patterns(R, "y", "w", keys = "g", compress = TRUE)
    expect_identical(colnames(prk$keys), "g")
    expect_true(same_patterns(sorted(prk), sorted(pk)))
    pr0 <- design_patterns(R2, "y", "w0", compress = TRUE)
    expect_true(same_patterns(sorted(pr0), sorted(p0)))

    ## SQL expressions, as in memory
    prs <- design_patterns(R, DBI::SQL("2 * y"), DBI::SQL("w + 1"), keys = DBI::SQL("g + 1"), compress = TRUE)
    pms <- design_patterns(D, 2 * d$y, d$w + 1, keys = d$g + 1, compress = TRUE)
    expect_true(same_patterns(sorted(prs), sorted(pms)))

    ## design_column()s of the design, and products with it, as in memory
    b <- coef(lm(fo, d))
    prc <- design_patterns(R, design_column(R, "y") - R %*% b, 2 * design_column(R, "w"),
                           keys = design_column(R, "g"), compress = TRUE)
    pmc <- design_patterns(D, d$y - D %*% b, 2 * d$w, keys = d$g, compress = TRUE)
    expect_true(same_patterns(sorted(prc), sorted(pmc)))
    expect_true(same_patterns(sorted(pmc), reference(X, d$y - drop(X %*% b), 2 * d$w, cbind(d$g))))

    ## lm on the remote patterns is lm on the rows
    lpr <- lm(prw$z ~ as.matrix(prw$x) - 1, weights = prw$w)
    expect_true(same(coef(lpr), coef(lw)))

    ## NULL and errors, as in memory
    expect_null(design_patterns(R, "y", compress = FALSE))
    Rc <- suppressMessages(design_matrix(y ~ xc, con, table = "d"))
    expect_null(design_patterns(Rc, "y", compress = TRUE))
    expect_error(design_patterns(R, "nope", compress = TRUE), "no column")
    expect_error(design_patterns(R2, "yna", compress = TRUE), "missing or infinite")
    expect_error(design_patterns(R2, "y", "wneg", compress = TRUE), "non-negative")
    expect_error(design_patterns(R, "f", compress = TRUE), "not numeric")
    test_con_done(con)
}

## ---- plan_of() and design_parts() ----

set.seed(3)
tr <- data.frame(y = rnorm(50), x = runif(50, 1, 2), f = factor(sample(c("u", "v", "w"), 50, TRUE)))
tr$yb <- as.integer(tr$y > 0)
nd <- data.frame(x = c(1.2, 1.9, 1.5), f = factor(c("w", "u", "w"), levels = levels(tr$f)))
fo2 <- y ~ poly(x, 2) + f

## on a formula: the matrix of model.matrix() on the training terms
pf <- plan_of(fo2, tr)
expect_inherits(pf, "xtx_plan")
mf <- model.frame(pf$terms, nd, xlev = pf$xlevels)
expect_true(strict_identical(as.matrix(design_matrix(pf, nd)),
                             model.matrix(pf$terms, mf, contrasts.arg = pf$contrasts)))
expect_identical(plan_of(terms(fo2), tr)$xlevels, pf$xlevels)

## on lm and glm fits: the columns predict() uses
fl <- lm(fo2, tr)
parts <- design_parts(fl)
expect_true(all(c("terms", "xlevels", "contrasts") %in% names(parts)))
expect_inherits(parts$terms, "terms")
expect_identical(parts$xlevels, fl$xlevels)
pl <- plan_of(fl)
expect_inherits(pl, "xtx_plan")
expect_true(same(drop(as.matrix(design_matrix(pl, nd)) %*% coef(fl)), predict(fl, nd), 1e-12))
fg <- glm(yb ~ poly(x, 2) + f, binomial, tr)
expect_identical(design_parts(fg)$xlevels, fg$xlevels)
pg <- plan_of(fg)
expect_true(same(drop(as.matrix(design_matrix(pg, nd)) %*% coef(fg)), predict(fg, nd), 1e-12))
## the same plan as design_matrix(fit, newdata = NULL)
expect_true(strict_identical(as.matrix(design_matrix(pg, nd)),
                             as.matrix(design_matrix(design_matrix(fg, newdata = NULL), nd))))

## a toy class with its own design_parts() method works with plan_of()
toy <- structure(list(formula = y ~ log(x) + f, levels = list(f = levels(tr$f))), class = "toyfit")
design_parts.toyfit <- function(model, ...) {
    list(terms = terms(model$formula), xlevels = model$levels, contrasts = list(f = "contr.sum"))
}
.S3method("design_parts", "toyfit", design_parts.toyfit)
expect_identical(design_parts(toy)$contrasts, list(f = "contr.sum"))
pt <- plan_of(toy)
expect_inherits(pt, "xtx_plan")
Xt <- as.matrix(design_matrix(pt, nd))
expect_true(same(c(Xt), c(model.matrix(~ log(x) + f, nd, contrasts.arg = list(f = "contr.sum")))))
expect_identical(colnames(Xt), c("(Intercept)", "log(x)", "f1", "f2"))
## and so does design_matrix(fit, newdata = ...)
expect_true(same(c(as.matrix(design_matrix(toy, newdata = nd))), c(Xt)))
## a method that returns something else is caught
design_parts.badfit <- function(model, ...) list(terms = "y ~ x", xlevels = NULL, contrasts = NULL)
.S3method("design_parts", "badfit", design_parts.badfit)
expect_error(plan_of(structure(list(), class = "badfit")), "terms")

## design_relation(): the rows of a remote design and columns over them as
## what a query reads (connection, source, filter, expressions)
if (requireNamespace("duckdb", quietly = TRUE) && requireNamespace("DBI", quietly = TRUE)) {
    con <- DBI::dbConnect(duckdb::duckdb())
    set.seed(5)
    dr <- data.frame(y = rnorm(30), x = rnorm(30), id = rep(1:6, 5), g = rep(c("a", "b", "c"), 10))
    dr$x[4] <- NA
    DBI::dbWriteTable(con, "dr", dr)
    D <- design_matrix(~ x, con, table = "dr")
    r <- design_relation(D, list(xb = D %*% c(1, 2), key = design_key(D, c("id", "g"))))
    expect_identical(names(r$sql), c("xb", "key"))
    expect_identical(lengths(r$sql), c(xb = 1L, key = 2L))
    got <- DBI::dbGetQuery(r$con, sprintf("SELECT %s AS id, SUM(%s) AS s FROM (%s) AS s WHERE %s GROUP BY 1 ORDER BY 1",
                                          r$sql$key[1L], r$sql$xb, r$select, r$filter))
    ok <- !is.na(dr$x)
    expect_equal(got$s, as.vector(tapply(1 + 2 * dr$x[ok], dr$id[ok], sum)), tolerance = 1e-12)
    expect_error(design_relation(D, list(a = 1)), "must be a column")
    expect_error(design_relation(design_matrix(~ x, dr), list(a = 1)))
    DBI::dbDisconnect(con)
}
