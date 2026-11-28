source("helpers.R")
library(xtx)

# In memory, qr() of a design with numeric variables far from zero factors
# the X'X of the centered design (as remote designs do), so that R, hence
# vcov(), keeps the accuracy of lm() on badly conditioned designs. The
# columns kept and the NA pattern still come from the uncentered X'X, as in
# lm(). (The fits on these designs are tested in xtx.fit.)
same <- function(a, b, tol = 1e-7) isTRUE(all.equal(unname(a), unname(b), tolerance = tol))

set.seed(42)
n <- 3000
d <- data.frame(
    xd = sample(1:5, n, TRUE),
    xc = rnorm(n),
    xe = rnorm(n),
    f1 = sample(c("b", "a", "c", "d"), n, TRUE),
    f2 = factor(sample(c("u", "v", "w"), n, TRUE), levels = c("w", "u", "v", "unused")),
    lg = sample(c(TRUE, FALSE), n, TRUE),
    w = runif(n, 0.5, 2))
d$y <- d$xd + d$xc + as.integer(factor(d$f1)) + rnorm(n)
d$xs <- 100 * d$xc + 1e7  # far from zero, not aliased by lm()
d$xt <- d$xc + 1e7        # so close to aliased that lm() drops it
d$xd7 <- d$xd + 1e7
d$xb <- 1e7 + 10 * d$xe   # ~1e7 with a small variance
d$f3 <- d$f1              # aliased with f1
d$y[c(3, 50)] <- NA
d$f1[7] <- NA
d$xc[9] <- NaN
d$w[11] <- 0

# qr() directly: rank, pivot and the factor of lm()'s QR, up to sign
d0 <- na.omit(d)
for (fo in list(y ~ xs + f1, y ~ xb * xs, y ~ xs * f2 + f1, y ~ xt + f1, y ~ xb * f2 + xt)) {
    info <- deparse(fo)
    D <- design_matrix(fo, d0)
    X <- as.matrix(D)
    q <- qr(D)
    q0 <- qr(X)
    expect_identical(q$rank, q0$rank, info = paste(info, "rank"))
    expect_identical(q$pivot, q0$pivot, info = paste(info, "pivot"))
    k <- q$rank
    i <- seq_len(k)
    expect_true(same(abs(qr.R(q)[i, i]), abs(qr.R(q0)[i, i])), info = paste(info, "R"))
    V <- chol2inv(qr.R(q)[i, i])
    V0 <- chol2inv(qr.R(q0)[i, i])
    expect_true(same(V, V0), info = paste(info, "(R'R)^-1"))
    expect_true(same(qr.coef(q, d0$y), qr.coef(q0, d0$y)), info = paste(info, "coef"))
    expect_true(same(qr.fitted(q, d0$y), qr.fitted(q0, d0$y)), info = paste(info, "fitted"))
}

# variables computed from a numeric variable far from zero (I(x^2) is a
# column of the model frame, not rebuilt from the shifted x)
set.seed(3)
dc <- data.frame(x = runif(5000) + 100, f = gl(5, 1000))
dc$y <- 1 + dc$x + rnorm(5000)
for (fo in list(y ~ f + x + I(x^2), y ~ f * log(x) + x)) {
    D <- design_matrix(fo, dc)
    expect_true(same(qr.coef(qr(D), dc$y), qr.coef(qr(as.matrix(D)), dc$y), 1e-8), info = deparse(fo))
}

# (design.matrix#r6g7) the shifted design is X T on every row: T is found
# on prototypes built, like the shifted design, from the frame's columns as
# stored; such designs are now centered (no guard)
ns <- asNamespace("xtx")
for (fo in list(y ~ f + x + I(x^2), y ~ f * log(x) + x, y ~ I(x^2) + f, y ~ log(x) * I(x^2))) {
    info <- paste(deparse(fo), "X T")
    D <- design_matrix(fo, dc)
    ctr <- ns$xtx_p_center(D)
    expect_false(is.null(ctr), info = info)
    X <- as.matrix(D)
    err <- abs(X %*% ctr$T - as.matrix(ctr$x))
    scale <- abs(X) %*% abs(ctr$T)
    expect_true(all(err <= 1e-10 * pmax(scale, 1)), info = info)
    ## (badly conditioned: qr() centers; I(x^2) alone is not bad enough)
    if (!identical(fo, y ~ I(x^2) + f)) expect_false(is.null(qr(D)$centered), info = paste(info, "centered"))
}

# centering leaves the user's random numbers alone
set.seed(1)
invisible(qr(design_matrix(y ~ xb * xs, d0)))
a <- runif(1)
set.seed(1)
expect_identical(runif(1), a)

# well-conditioned designs are not centered (no extra pass)
q <- qr(design_matrix(y ~ xc * f1 + xd, d0))
expect_null(q$centered)
q <- qr(design_matrix(y ~ xb * xs, d0))
expect_false(is.null(q$centered))

# design_center(): the transform T (X~ = X T) is exact in every entry. With
# x and z near 1e7, the column of x:z has an entry near 1e14 (the product
# of the shifts) and entries near 1e7 (minus each shift): products of the
# shifts, not their least squares estimates, which would give the small
# ones nine digits (and shift the fitted values of coefficients mapped
# back with them).
set.seed(7)
dt <- data.frame(x = 1e7 + rnorm(500), z = 1e7 + 10 * rnorm(500), f = gl(2, 250))
Xt <- design_matrix(~ x * z + f, dt)
Tm <- design_center(Xt)$T
cn <- colnames(Xt)
sx <- -Tm[1L, match("x", cn)]
sz <- -Tm[1L, match("z", cn)]
j <- match("x:z", cn)
expect_true(abs(sx - 1e7) < 1 && abs(sz - 1e7) < 10)
expect_equal(Tm[match("x", cn), j], -sz, tolerance = 1e-13)
expect_equal(Tm[match("z", cn), j], -sx, tolerance = 1e-13)
expect_equal(Tm[1L, j], sx * sz, tolerance = 1e-13)
