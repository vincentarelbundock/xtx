# A purely numeric design: 118 numeric columns and an intercept, as many
# columns as the design with factors of the other benchmarks. Returns one
# row per operation, starting point and package: the median time in
# milliseconds.
p <- 118
num <- as.data.frame(matrix(rnorm(n * p), n, p))
names(num) <- paste0("x", seq_len(p))
fn <- reformulate(names(num), "y")
num$y <- rnorm(n)
Xn <- model.matrix(fn, num)
Dn <- design_matrix(fn, num)

beta <- rnorm(ncol(Xn))
w <- rnorm(n)
g <- sample(1:20, n, TRUE)
V <- crossprod(matrix(rnorm(ncol(Xn)^2), ncol(Xn))) / ncol(Xn)
base_diagonal <- function(X, V) rowSums((X %*% V) * X)

# Each operation: base R then xtx from the data frame, then (but for the
# first) base R then xtx from the prebuilt matrix
time <- function(operation, ..., check = close, iterations = 10) {
    res <- bench::mark(..., iterations = iters(iterations), check = check)
    data.frame(operation = operation,
               start = rep(c("From a data frame", "From a prebuilt design matrix"), each = 2)[seq_len(nrow(res))],
               package = rep(c("Base R, dense matrix", "xtx, design matrix"), length.out = nrow(res)),
               ms = signif(as.numeric(res$median) * 1000, 3))
}
rbind(
    time("Build the dense matrix",
         model.matrix(fn, num),
         as.matrix(design_matrix(fn, num)), check = same),
    time("%*%",
         model.matrix(fn, num) %*% beta,
         design_matrix(fn, num) %*% beta,
         Xn %*% beta,
         Dn %*% beta),
    time("crossprod(D, w)",
         crossprod(model.matrix(fn, num), w),
         crossprod(design_matrix(fn, num), w),
         crossprod(Xn, w),
         crossprod(Dn, w)),
    time("crossprod(D)",
         crossprod(model.matrix(fn, num)),
         crossprod(design_matrix(fn, num)),
         crossprod(Xn),
         crossprod(Dn)),
    time("colSums(w * D, group = )",
         rowsum(w * model.matrix(fn, num), g),
         colSums(w * design_matrix(fn, num), group = g),
         rowsum(w * Xn, g),
         colSums(w * Dn, group = g)),
    time("quadratic_diagonal()",
         base_diagonal(model.matrix(fn, num), V),
         quadratic_diagonal(design_matrix(fn, num), V),
         base_diagonal(Xn, V),
         quadratic_diagonal(Dn, V), iterations = 5))
