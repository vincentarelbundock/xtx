source("helpers.R")
library(xtx)

set.seed(3)
d <- data.frame(x = rnorm(40), f = gl(4, 10), g = sample(c(3, 1, 2), 40, TRUE))
D <- design_matrix(~ x * f, d)
X <- model.matrix(~ x * f, d)
w <- runif(40)

# rowsum(): the groups, their order and the row names of base::rowsum()
expect_equal(colSums(D, group = d$g), rowsum(X, d$g))
expect_identical(dimnames(colSums(D, group = d$g)), dimnames(rowsum(X, d$g)))
expect_equal(colSums(D, group = d$g, reorder = FALSE), rowsum(X, d$g, reorder = FALSE))
expect_identical(rownames(colSums(D, group = d$g, reorder = FALSE)), rownames(rowsum(X, d$g, reorder = FALSE)))
fg <- factor(sample(c("b", "a", "c"), 40, TRUE), levels = c("c", "b", "a", "unused"))
expect_identical(dimnames(colSums(D, group = fg)), dimnames(rowsum(X, fg)))
expect_equal(colSums(D, group = fg), rowsum(X, fg))
cg <- sample(c("u", "v"), 40, TRUE)
expect_equal(colSums(D, group = cg), rowsum(X, cg))

# missing groups: a row of their own, with base R's warning
ng <- d$g
ng[c(2, 7)] <- NA
expect_warning(colSums(D, group = ng), "missing values for 'group'")
expect_equal(suppressWarnings(colSums(D, group = ng)), suppressWarnings(rowsum(X, ng)))
expect_identical(rownames(suppressWarnings(colSums(D, group = ng))), rownames(suppressWarnings(rowsum(X, ng))))

# weights through arithmetic, on one design and on a combination of two
expect_equal(colSums(w * D, group = d$g), rowsum(w * X, d$g))
D1 <- design_matrix(~ x * f, transform(d, x = x + 1))
X1 <- model.matrix(~ x * f, transform(d, x = x + 1))
expect_equal(colSums(w * D1 - D, group = d$g), rowsum(w * X1 - X, d$g), check.attributes = FALSE)

# group means
expect_equal(colSums((1 / tabulate(d$f)[d$f]) * D, group = d$f), rowsum(X, d$f) / as.vector(table(d$f)))

# colSums() and friends, with the names of the dense versions
expect_equal(colSums(D), colSums(X))
expect_equal(colMeans(D), colMeans(X))
expect_equal(rowSums(D), rowSums(X))
expect_equal(rowMeans(D), rowMeans(X))
expect_identical(names(colSums(D)), colnames(X))
expect_identical(names(rowSums(D)), rownames(X))
expect_equal(colMeans(w * D1 - D), colMeans(w * X1 - X), check.attributes = FALSE)

# zero rows and zero columns
D0 <- design_matrix(~ x * f, d[0, ])
expect_identical(colSums(D0), colSums(model.matrix(~ x * f, d[0, ])))
expect_identical(dim(colSums(D0, group = integer(0))), c(0L, ncol(D0)))

# Dense fallback, including combinations, delegates grouping and NA handling to base.
dense <- design_matrix(lm(x ~ f, d))
fg[1L] <- NA
for (J in list(dense, -dense, w * dense, dense - w * dense)) {
    for (reorder in c(FALSE, TRUE)) {
        expect_identical(suppressWarnings(colSums(J, group = fg, reorder = reorder)),
                         suppressWarnings(rowsum(as.matrix(J), fg, reorder)))
    }
}
expect_equal(colSums(-D, group = d$g), rowsum(-X, d$g))
expect_identical(dim(colSums(-D0, group = integer(0))), c(0L, ncol(D0)))
for (X in list(matrix(.Machine$integer.max, 3L), matrix(TRUE, 3L))) {
    dense <- xtx:::xtx_dense_design(X)
    expect_identical(unname(colSums(dense, group = c(1L, 1L, 1L))), matrix(sum(as.double(X))))
}

# Grouped column summaries keep rowsum ordering and omit unused levels.
for (J in list(D, w * D1 - D, dense)) {
    group <- rep(fg, length.out = nrow(J))
    ref <- suppressWarnings(rowsum(as.matrix(J) * 1.0, group))
    counts <- suppressWarnings(rowsum(rep(1, nrow(J)), group))
    expect_equal(suppressWarnings(colSums(J, group = group)), ref)
    expect_equal(suppressWarnings(colMeans(J, group = group)), ref / as.vector(counts))
}
expect_identical(dim(colMeans(D0, group = integer())), c(0L, ncol(D0)))

expect_silent(rowsum(D, d$g))
expect_identical(rowsum(D, d$g, reorder = FALSE),
                 colSums(D, group = d$g, reorder = FALSE))

# ---- the bound of the structured grouped sums (xtx_colsums_struct_ok()):
# relative to the G x p output, so that a grouping with millions of levels
# (a fixed effect at 5 million levels) keeps the structured kernel; the
# check allocates nothing
local({
    dat <- as.data.frame(matrix(rnorm(50), 5, 10))
    dat$y <- rnorm(5)
    D <- design_matrix(reformulate(paste0("V", 1:10), "y"), dat)
    ok <- getFromNamespace("xtx_colsums_struct_ok", "xtx")
    expect_true(ok(D, 1e6))
    expect_true(ok(D, 5e6))
    expect_true(ok(D, c(5e6, 10)))
})
