## quadratic_diagonal() and grouped column sums through the structured
## kernels: the dense columns (intercept, numeric variables) and the terms
## with factors are handled separately (xtx_struct_quad(), and the kernel
## of design_group_sums() for colSums(, group = )). Results agree with the
## dense computation and with the row kernels up to rounding, and a row's
## result does not depend on the other rows.
library(xtx)
set.seed(11)
n <- 3000
d <- data.frame(y = rnorm(n), x = rnorm(n), z = runif(n),
                f1 = factor(sample(letters[1:6], n, TRUE)),
                f2 = factor(sample(LETTERS[1:4], n, TRUE)),
                f3 = factor(sample(c("u", "v", "w"), n, TRUE)))
forms <- list(y ~ x + z, y ~ 0 + x + z, y ~ f1, y ~ 0 + f1, y ~ x * f1 + f2, y ~ x * z * f1,
              y ~ f1:f2:f3 + x, y ~ poly(x, 3) * f1 + z:f2, y ~ f1 * f2 * f3 + poly(z, 2):f1)
close_rel <- function(a, b) max(abs(a - b) / (1 + abs(b))) < 1e-11
w <- runif(n, 0.5, 2)
g <- sample(c("a", "b", "c", NA), n, TRUE)
rows_kernel <- function(code) {
    op <- options(xtx.structured = FALSE)
    on.exit(options(op))
    code()
}
op <- options(xtx.dedup = FALSE)
for (contr in c("contr.treatment", "contr.sum")) {
    oc <- options(contrasts = c(contr, "contr.poly"))
    for (f in forms) {
        X <- model.matrix(f, d)
        D <- design_matrix(f, d)
        p <- ncol(X)
        V <- matrix(rnorm(p * p), p)        # any finite matrix: not symmetric
        ref <- unname(rowSums((X %*% V) * X))
        st <- xtx:::xtx_struct_single(xtx:::as_jac(D))
        expect_false(is.null(st))
        expect_true(close_rel(xtx:::xtx_struct_quad(st$design, V), ref))
        expect_true(close_rel(quadratic_diagonal(D, V), ref))
        expect_true(close_rel(quadratic_diagonal(D, V), rows_kernel(function() quadratic_diagonal(D, V))))
        expect_true(close_rel(quadratic_diagonal(w * D, V), w^2 * ref))
        ## grouped sums, with and without weights
        expect_true(close_rel(suppressWarnings(colSums(D, group = g)), suppressWarnings(rowsum(X, g))))
        expect_true(close_rel(suppressWarnings(colSums(w * D, group = g)), suppressWarnings(rowsum(w * X, g))))
        expect_identical(dimnames(suppressWarnings(colSums(w * D, group = g))),
                         dimnames(rows_kernel(function() suppressWarnings(colSums(w * D, group = g)))))
        expect_true(close_rel(suppressWarnings(colMeans(D, group = g)),
                              suppressWarnings(rowsum(X, g)) / as.vector(table(g, useNA = "ifany"))))
    }
    options(oc)
}

## a row's result is the same bit for bit whatever its position and the
## other rows, with any number of threads
f <- y ~ poly(x, 3) * f1 + z:f2 + z
D <- design_matrix(f, d)
V <- crossprod(matrix(rnorm(ncol(D)^2), ncol(D)))
q <- quadratic_diagonal(D, V)
i <- sample(n, 700)
expect_identical(quadratic_diagonal(D[i, ], V), q[i])
dm <- cbind(d, matrix(rnorm(n * 20), n, dimnames = list(NULL, paste0("v", 1:20))))
many <- reformulate(c(paste0("v", 1:20), "f1"), "y")   # enough dense columns for the BLAS
Dm <- design_matrix(many, dm)
Vm <- crossprod(matrix(rnorm(ncol(Dm)^2), ncol(Dm))) / ncol(Dm)
qm <- quadratic_diagonal(Dm, Vm)
expect_identical(quadratic_diagonal(Dm[i, ], Vm), qm[i])
expect_true(close_rel(qm, rowSums((as.matrix(Dm) %*% Vm) * as.matrix(Dm))))
ot <- options(xtx.threads = 1L)
q1 <- quadratic_diagonal(Dm, Vm)
options(xtx.threads = 3L)
expect_identical(quadratic_diagonal(Dm, Vm), q1)
options(ot)

## with groups of identical rows: the same results bit for bit
dd <- d[sample(60, n, TRUE), ]
Dd <- design_matrix(f, dd)
options(xtx.dedup = TRUE)
q_on <- quadratic_diagonal(Dd, V)
expect_true(is.list(Dd$cache$groups))
options(xtx.dedup = FALSE)
expect_identical(q_on, quadratic_diagonal(design_matrix(f, dd), V))

## tables too large for the cap: NULL, and the row kernel answers
st <- xtx:::xtx_struct_single(xtx:::as_jac(D))
expect_null(xtx:::xtx_struct_quad(st$design, V, cap = 10))

## non-finite values are left to the row kernel: NaN where the dense
## computation has them
dn <- d
dn$x[5] <- Inf
Dn <- design_matrix(y ~ x * f1, dn)
Xn <- model.matrix(y ~ x * f1, dn)
Vn <- diag(ncol(Xn))
expect_identical(is.nan(quadratic_diagonal(Dn, Vn)), unname(is.nan(rowSums((Xn %*% Vn) * Xn))))
expect_equal(suppressWarnings(colSums(Dn, group = g)), suppressWarnings(rowsum(Xn, g)), check.attributes = FALSE)
options(op)
