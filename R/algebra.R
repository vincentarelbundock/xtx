## the S3 class, registered for S4 methods (colSums() and friends, as())
#' @importFrom methods setOldClass
NULL
setOldClass("xtx_matrix")

#' Arithmetic on design matrices
#'
#' Scaling rows and adding or subtracting design matrices gives another
#' `design_matrix`, still without forming the dense matrix. Row `i` of the
#' result is `s1[i] * x1_i + s2[i] * x2_i`, where `x1_i` and `x2_i` are rows of
#' at most two designs with the same rows and columns.
#'
#' * `s * D` and `D * s`: scale the rows of `D` by `s`, a number or a numeric
#'   vector with one value per row.
#' * `D / s`: the same with `1 / s`.
#' * `-D`: the same with `-1`.
#' * `D1 + D2` and `D1 - D2`: add or subtract two design matrices with the same
#'   rows and columns.
#'
#' This covers the Jacobians used for delta-method standard errors: the
#' Jacobian of predictions from a model with a linear predictor is
#' `dmu * D`, where `dmu` is the derivative of the inverse link, and the
#' Jacobian of a comparison between two counterfactual data sets is
#' `dmu1 * D1 - dmu0 * D0`. It also covers weighted designs, such as
#' `w * D` in `crossprod(D, w * D)`.
#'
#' Scales are combined before any product is taken: `as.matrix(s2 * (s1 * D))`
#' is `(s2 * s1) * as.matrix(D)`. A result that would combine more than two
#' designs is an error.
#'
#' @param e1,e2 A `design_matrix`, or a number or numeric vector for `*`
#'   and `/`.
#' @return A `design_matrix`.
#' @examples
#' d <- data.frame(x = rnorm(50), f = gl(5, 10))
#' p <- design_matrix(~ x * f, d)
#' X0 <- design_matrix(p, d)
#' X1 <- design_matrix(p, transform(d, x = x + 1))
#' J <- X1 - X0
#' V <- diag(10) / 100
#' quadratic_diagonal(J, V)
#' @name design_arithmetic
NULL

#' @rdname design_arithmetic
#' @export
Ops.xtx_matrix <- function(e1, e2) {
    op <- .Generic
    if (nargs() == 1L) {
        if (op == "+") return(e1)
        if (op == "-") return(xtx_scale(e1, -1))
    } else if (op %in% c("*", "/")) {
        d1 <- inherits(e1, "xtx_matrix")
        d2 <- inherits(e2, "xtx_matrix")
        if (d1 && d2)
            stop(sprintf("`%s` between two design matrices is not supported.", op), call. = FALSE)
        if (op == "/" && d2)
            stop("Only `design / number` is supported, with the design first.", call. = FALSE)
        s <- if (d1) e2 else e1
        if (op == "/") s <- 1 / s
        return(xtx_scale(if (d1) e1 else e2, s))
    } else if (op %in% c("+", "-")) {
        if (!inherits(e1, "xtx_matrix") || !inherits(e2, "xtx_matrix"))
            stop(sprintf("`%s` needs two design matrices.", op), call. = FALSE)
        t2 <- xtx_terms(e2)
        if (op == "-") t2 <- lapply(t2, function(t) { t$scale <- xtx_times(t$scale, -1); t })
        return(xtx_from_terms(c(xtx_terms(e1), t2)))
    }
    stop(sprintf("`%s` is not supported for design matrices.", op), call. = FALSE)
}

## Is `x` a combination of scaled designs (rather than one plain design)?
xtx_is_comb <- function(x) !is.null(x[["design"]])

## A design matrix as a list of terms, each a plain design and its row scale
## (NULL for 1).
xtx_terms <- function(x) {
    if (!xtx_is_comb(x)) return(list(list(design = x, scale = NULL)))
    out <- list(list(design = x$design, scale = x$scale))
    if (!is.null(x$design2)) out[[2L]] <- list(design = x$design2, scale = x$scale2)
    out
}

## The error for a scale (or the weights of a fit, design_problem() of xtx.fit) that
## is not one value per row of a design of n rows.
xtx_scale_length_msg <- function(n) {
    sprintf("A design matrix can only be scaled by a number or a numeric vector of length %d (one value per row).", n)
}

## scale * s, where a NULL scale is 1 and s is a number or a vector of length n
xtx_times <- function(scale, s) {
    s <- as.double(s)
    if (is.null(scale)) s else scale * s
}

xtx_scale <- function(x, s) {
    n <- nrow(x)
    if (!is.numeric(s) || !(length(s) %in% c(1L, n)) || !is.null(dim(s)))
        stop(xtx_scale_length_msg(n), call. = FALSE)
    terms <- lapply(xtx_terms(x), function(t) { t$scale <- xtx_times(t$scale, s); t })
    xtx_from_terms(terms)
}

xtx_from_terms <- function(terms) {
    if (length(terms) > 2L)
        stop("A design matrix can combine at most two designs.", call. = FALSE)
    t1 <- terms[[1L]]
    if (length(terms) == 1L && is.null(t1$scale)) return(t1$design)
    t2 <- if (length(terms) == 2L) terms[[2L]] else list(design = NULL, scale = NULL)
    xtx_combine(t1$design, t1$scale, t2$design, t2$scale)
}

## Row i is scale[i] * x1_i + scale2[i] * x2_i, for plain designs x1 and x2
## (x2 optional) with the same rows and columns.
xtx_combine <- function(design, scale = NULL, design2 = NULL, scale2 = NULL) {
    if (!is.null(design2) &&
        (nrow(design2) != nrow(design) || !identical(colnames(design2), colnames(design))))
        stop("The design matrices must have the same rows and columns.", call. = FALSE)
    structure(list(design = design, scale = scale,
                   design2 = design2, scale2 = scale2,
                   cache = new.env(parent = emptyenv())),
              class = "xtx_matrix")
}

## the dense matrix of a combination
xtx_comb_dense <- function(x) {
    if (jac_kernel_ok(x)) {
        J <- do.call(.Call, c(list(C_xtx_jac_dense), jac_args(x), list(design_threads())))
        ## with a NaN anywhere, use R's arithmetic, whose NaN payloads
        ## (NA or NaN) the kernel cannot guarantee to reproduce
        if (is.null(attr(J, "xtx_nan"))) {
            dimnames(J) <- dimnames(x$design)
            return(J)
        }
    }
    X <- as.matrix(x$design)
    J <- if (is.null(x$scale)) X else x$scale * X
    if (!is.null(x$design2)) {
        X2 <- as.matrix(x$design2)
        J <- J + (if (is.null(x$scale2)) X2 else x$scale2 * X2)
    }
    attributes(J) <- attributes(X)[c("dim", "dimnames")]
    J
}

## Can the C kernels handle this object? Otherwise use dense algebra.
kernel_ok <- function(d) {
    if (is.null(d$setup)) return(FALSE)
    s <- d$setup
    if (s$nterms == 0L) return(TRUE)
    used <- (s$factors != 0) & (s$columns != 0)
    max(colSums(used)) <= 32L
}

as_jac <- function(x) {
    if (xtx_is_comb(x)) return(x)
    J <- xtx_combine(x)
    J$cache <- x$cache # a design's groups are those of its plain combination
    J
}

jac_args <- function(J) {
    d2 <- J$design2
    list(J$design$setup, J$design$vars, J$scale,
         if (is.null(d2)) NULL else d2$setup,
         if (is.null(d2)) NULL else d2$vars,
         J$scale2)
}

jac_kernel_ok <- function(J) {
    kernel_ok(J$design) && (is.null(J$design2) || kernel_ok(J$design2))
}

#' Algebra on design matrices without forming the dense matrix
#'
#' Products with a `design_matrix` are computed directly from
#' the factor codes and numeric variables, with the usual base R operators:
#'
#' * `x %*% beta`: `X %*% beta`, an `nrow(x)` by `ncol(beta)` matrix.
#' * `crossprod(x, w)`: `t(X) %*% w`, an `ncol(x)` by `ncol(w)` matrix.
#' * `crossprod(x)` and `crossprod(x, y)` with `y` another `design_matrix`
#'   with the same rows: `t(X) %*% X` and `t(X) %*% Y`, ordinary base R
#'   matrices. Weighted crossproducts are `crossprod(x, w * x)` or
#'   `crossprod(sqrt(w) * x)`. The result can go straight to [chol()],
#'   [chol2inv()], [eigen()] or [solve()]; see also [qr()][design_qr] for
#'   least squares.
#' * `colSums(x, group = g)` and `colMeans(x, group = g)` summarize each
#'   design column within groups of rows. Group order, row names and missing
#'   group handling follow [base::rowsum()]. For weighted sums, use
#'   `colSums(w * x, group = g)`. For design matrices, `rowsum(x, g)` is
#'   deprecated in favor of `colSums(x, group = g)`.
#' * `colSums(x)`, `colMeans(x)`, `rowSums(x)` and `rowMeans(x)`: as for the
#'   dense matrix (without `na.rm` or `dims`).
#' * `quadratic_diagonal(x, V)`: `rowSums((X %*% V) * X)`, the diagonal of
#'   `X V X'`. For a Jacobian and a coefficient covariance matrix, these are
#'   the variances of unit-level estimates.
#'
#' The results agree with the same computation on the dense matrix up to
#' floating point rounding (the order of the sums differs), with `NA` and
#' `NaN` in the same places. Only `design %*% beta`, `crossprod(design, w)`
#' and `crossprod(design, design)` are supported, with a design first.
#'
#' `beta` and `V` must be finite: the kernels skip structural zeros, so
#' `0 * Inf` terms that the dense product would turn into `NaN` are not
#' computed.
#'
#' Designs that the kernels cannot handle (non-formula objects, terms with
#' more than 32 variables) fall back to dense algebra.
#'
#' **Identical rows.** Many data sets have far fewer distinct covariate
#' patterns than rows. The kernels can group rows whose inputs (factor
#' codes, numeric values, scales and group) have the same bits, compute
#' once per pattern, and expand. Row-wise results (`%*%`, `quadratic_diagonal()`)
#' are then identical bit for bit to computing every row; `crossprod()` and
#' grouped column summaries first sum the weights within each pattern, which changes
#' only the rounding. The groups are found once per object and cached.
#' The option `xtx.dedup` controls this:
#' `"auto"` (default) groups when the work per row is large enough to pay for
#' hashing it and at most half the rows are distinct, `TRUE` always groups and
#' `FALSE` never does.
#'
#' @param x A `design_matrix`.
#' @param group A vector or factor with one value per row, as in
#'   [base::rowsum()].
#' @param reorder As in [base::rowsum()]: if `TRUE`, the groups are sorted;
#'   otherwise they are in order of first appearance.
#' @param na.rm,dims Only the defaults are supported.
#' @param V A square numeric matrix with `ncol(x)` rows.
#' @return `rowsum()`: a matrix with one row per group and the columns of
#'   `x`. `colSums()` and `colMeans()`: a vector with one value per column, or with
#'   `group`, a matrix with one row per observed group (sorted by default).
#'   `rowSums()`, `rowMeans()` and `quadratic_diagonal()`: a vector with one
#'   value per row (for a remote design, a column: see [design_column]).
#' @examples
#' d <- data.frame(x = rnorm(100), f = gl(4, 25), g = gl(2, 1, 100))
#' D <- design_matrix(~ x * f, d)
#' X <- model.matrix(~ x * f, d)
#' b <- rnorm(ncol(D))
#' all.equal(D %*% b, X %*% b, check.attributes = FALSE)
#' crossprod(D, rep(1, 100))
#' colMeans(D, group = d$g)
#' colSums(runif(100) * D, group = d$g)
#' colMeans(D)
#' quadratic_diagonal(D, diag(ncol(D)))
#' @name design_algebra
NULL

## X %*% B for a design matrix: an n x K matrix.
xtx_mult <- function(x, beta) {
    p <- ncol(x)
    if (is.matrix(beta)) {
        checkmate::assert_matrix(beta, mode = "numeric", nrows = p, any.missing = FALSE)
    } else {
        checkmate::assert_numeric(beta, len = p, any.missing = FALSE)
    }
    checkmate::assert_true(all(is.finite(beta)))
    J <- as_jac(x)
    B <- beta
    if (!is.matrix(B)) B <- matrix(B, ncol = 1L)
    storage.mode(B) <- "double"
    out <- NULL
    if (!is.null(st <- xtx_struct_single(J)) && xtx_in_range(B, 1e-100, 1e100)) {
        out <- xtx_struct_mult(st$design, B)
        if (!is.null(st$scale)) out <- st$scale * out
        ## the scale is applied after the product: if that order overflows,
        ## use the row kernel on the scaled design
        if (!xtx_all_finite(out)) out <- NULL
    }
    if (!is.null(out)) {
        ## done
    } else if (!jac_kernel_ok(J)) {
        out <- as.matrix(J) %*% B
    } else {
        g <- xtx_dedup_groups(J, cost = xtx_row_cost(J) * ncol(B))
        out <- do.call(.Call, c(list(C_xtx_jac_mult), jac_args(J),
                                list(B, if (is.null(g)) NULL else g$first, design_threads())))
        if (!is.null(g)) out <- out[g$map, , drop = FALSE]
    }
    unname(out)
}

## t(X) %*% W: a p x K matrix.
xtx_crossprod <- function(x, w) {
    n <- nrow(x)
    if (is.matrix(w)) {
        checkmate::assert_matrix(w, mode = "numeric", nrows = n)
    } else {
        checkmate::assert_numeric(w, len = n)
    }
    J <- as_jac(x)
    W <- w
    if (!is.double(W)) storage.mode(W) <- "double"
    st <- xtx_struct_single(J)
    out <- NULL
    if (!is.null(st) && xtx_in_range(W, 1e-100, 1e100)) {
        ## the scale moves to the weights (safe within the bounds); a
        ## vector becomes a one-column matrix after the product with the
        ## scale (one copy, not two)
        sW <- if (is.null(st$scale)) W else st$scale * W
        if (!is.matrix(sW)) dim(sW) <- c(n, 1L)
        out <- xtx_struct_xty(st$design, sW)
        if (!xtx_all_finite(out)) out <- NULL
    }
    if (is.null(out) && !is.matrix(W)) dim(W) <- c(n, 1L)
    out <- if (!is.null(out)) {
        out
    } else if (jac_kernel_ok(J)) {
        xtx_crossprod_kernel(J, W)
    } else {
        crossprod(as.matrix(J), W)
    }
    dimnames(out) <- list(colnames(x), colnames(w))
    out
}

## With group codes, accumulate each column of W separately for each group.
xtx_crossprod_kernel <- function(J, W, codes = NULL, G = 1L) {
    g <- xtx_dedup_groups(J, codes, cost = xtx_row_cost(J) * ncol(W))
    first <- NULL
    if (!is.null(g)) {
        W <- rowsum(W, g$map, reorder = TRUE)
        first <- g$first
        if (!is.null(codes)) codes <- codes[first]
    }
    do.call(.Call, c(list(C_xtx_jac_crossprod), jac_args(J),
                     list(W, codes, G, first, design_threads())))
}

#' @rdname design_algebra
#' @export
rowsum.xtx_matrix <- function(x, group, reorder = TRUE, na.rm = FALSE, ...) {
    xtx_defaults_only(na.rm = na.rm)
    xtx_group_summary(x, group, reorder)
}

xtx_group_summary <- function(x, group, reorder = TRUE, mean = FALSE) {
    n <- nrow(x)
    checkmate::assert_atomic_vector(group, len = n)
    checkmate::assert_flag(reorder)
    ## groups as in rowsum.default()
    if (anyNA(group)) warning("missing values for 'group'", call. = FALSE)
    ugroup <- unique(group)
    if (reorder) ugroup <- sort(ugroup, na.last = TRUE, method = "quick")
    codes <- match(group, ugroup)
    J <- as_jac(x)
    if (!jac_kernel_ok(J)) {
        X <- as.matrix(x)
        if (is.integer(X) || is.logical(X)) storage.mode(X) <- "double"
        out <- base::rowsum(X, codes)
    } else {
        ## a design with scaled rows: the scale is a weight, so that rows with
        ## the same inputs but different weights can still be grouped
        w <- 1
        if (is.null(J$design2) && !is.null(J$scale)) {
            w <- J$scale
            J <- as_jac(J$design)
        }
        G <- length(ugroup)
        out <- NULL
        ## A plain design within the bounds of the structured kernels: the
        ## sums by column and group in one pass over the columns (the kernel
        ## of design_group_sums()), whatever the number of entries per row.
        ## The row kernel below visits each entry of each row, which costs
        ## several times more when rows are dense (numeric variables).
        if (is.null(J$design2) && is.null(J$scale) && xtx_colsums_struct_ok(J$design, G) &&
            xtx_safe_scale(w)) {
            ws <- if (identical(w, 1)) NULL else rep_len(as.double(w), n)
            S <- xtx_gsums_memory(J$design, ws, list(list(G = G, codes = codes)))$gs[[1L]]$S
            if (xtx_all_finite(S)) out <- S
        }
        if (is.null(out)) {
            W <- matrix(w, nrow = n, ncol = 1L)
            out <- t(xtx_crossprod_kernel(J, W, codes, G))
        }
    }
    if (mean) out <- out / tabulate(codes, nbins = length(ugroup))
    dimnames(out) <- list(as.character(ugroup), colnames(x))
    out
}

xtx_defaults_only <- function(na.rm = FALSE, dims = 1) {
    if (!isFALSE(na.rm)) stop("`na.rm = TRUE` is not supported for design matrices.", call. = FALSE)
    if (!identical(dims, 1) && !identical(dims, 1L))
        stop("Only `dims = 1` is supported for design matrices.", call. = FALSE)
}

xtx_col_sums <- function(x) {
    out <- xtx_crossprod(x, rep(1, nrow(x)))[, 1L]
    names(out) <- colnames(x)
    out
}

xtx_row_sums <- function(x) {
    out <- xtx_mult(x, rep(1, ncol(x)))[, 1L]
    names(out) <- rownames(x)
    out
}

#' @rdname design_algebra
#' @importFrom methods setMethod
#' @export
setMethod("colSums", "xtx_matrix", function(x, na.rm = FALSE, dims = 1, ..., group = NULL, reorder = TRUE) {
    xtx_defaults_only(na.rm, dims)
    if (!is.null(group)) return(xtx_group_summary(x, group, reorder))
    xtx_col_sums(x)
})

#' @rdname design_algebra
#' @export
setMethod("colMeans", "xtx_matrix", function(x, na.rm = FALSE, dims = 1, ..., group = NULL, reorder = TRUE) {
    xtx_defaults_only(na.rm, dims)
    if (!is.null(group)) return(xtx_group_summary(x, group, reorder, mean = TRUE))
    xtx_col_sums(x) / nrow(x)
})

#' @rdname design_algebra
#' @export
setMethod("rowSums", "xtx_matrix", function(x, na.rm = FALSE, dims = 1, ...) {
    xtx_defaults_only(na.rm, dims)
    xtx_row_sums(x)
})

#' @rdname design_algebra
#' @export
setMethod("rowMeans", "xtx_matrix", function(x, na.rm = FALSE, dims = 1, ...) {
    xtx_defaults_only(na.rm, dims)
    xtx_row_sums(x) / ncol(x)
})

#' @rdname design_algebra
#' @export
quadratic_diagonal <- function(x, V) {
    checkmate::assert_multi_class(x, c("xtx_matrix", "xtx_remote"))
    p <- ncol(x)
    checkmate::assert_matrix(V, mode = "numeric", nrows = p, ncols = p, any.missing = FALSE)
    checkmate::assert_true(all(is.finite(V)))
    ## a vector in memory, a remote column (design_column) remotely
    if (inherits(x, "xtx_remote")) return(xtx_column_label(xtx_p_quad(x, V), ".quad"))
    xtx_p_output(x, list(.quad = xtx_p_quad(x, V)), memory = function(d) d$.quad)
}

## Groups of identical rows (same bits in every input of the row), or NULL.
##
## Option `xtx.dedup`: FALSE (never), TRUE (always) or "auto".
## Groups are computed once per object and cached (the object's inputs never
## change). In "auto" mode they are computed only when the kernel is costly
## per row compared with hashing a row (`cost`, an estimate of the work per
## row), and kept only when at most half the rows are distinct; the search
## stops as soon as that is exceeded, and the failure is cached too.
xtx_dedup_groups <- function(J, group = NULL, cost = Inf, force = FALSE) {
    mode <- getOption("xtx.dedup", "auto")
    if (isFALSE(mode) && !force) return(NULL)
    always <- force || isTRUE(mode)
    cache <- J$cache
    g <- if (!is.null(cache)) cache$groups else NULL
    if (isFALSE(g) && force) g <- NULL # an "auto" search gave up; search fully
    if (is.null(g)) {
        n <- nrow(J)
        g <- xtx_groups_from_designs(J)
        if (!is.null(g) && !always && length(g$first) > 0.5 * n) g <- FALSE
        if (is.null(g)) {
            if (!always && cost < 16) return(NULL)
            g <- do.call(.Call, c(list(C_xtx_row_groups), jac_args(J),
                                  list(NULL, if (always) 1 else 0.5)))
            if (is.null(g)) g <- FALSE else names(g) <- c("map", "first")
        }
        if (!is.null(cache)) cache$groups <- g
    }
    if (isFALSE(g)) return(NULL)
    if (!is.null(group)) {
        ## patterns within groups: combine the cached ids with the group codes
        group[is.na(group)] <- 0L
        key <- (g$map - 1) * (max(group) + 1) + group
        u <- unique(key)
        map <- match(key, u)
        first <- match(seq_along(u), map)
        g <- list(map = map, first = first)
    }
    g
}

## A Jacobian's groups from the cached groups of its designs: rows in the
## same group of every design, whose scales also have the same bits, give
## the same Jacobian row. NULL when the designs have no cached groups or a
## scale differs within a group (the caller then hashes the Jacobian).
xtx_groups_from_designs <- function(J) {
    if (!xtx_is_comb(J) || identical(J$cache, J$design$cache)) return(NULL)
    g1 <- J$design$cache$groups
    if (!is.list(g1)) return(NULL)
    map <- g1$map
    if (!is.null(J$design2)) {
        g2 <- J$design2$cache$groups
        if (!is.list(g2)) return(NULL)
        key <- (map - 1) * length(g2$first) + g2$map # exact: below 2^53
        map <- match(key, unique(key))
    }
    first <- match(seq_len(max(0L, map)), map)
    for (s in list(J$scale, J$scale2)) {
        if (length(s) > 1L &&
            !identical(s, s[first][map], num.eq = FALSE, single.NA = FALSE)) return(NULL)
    }
    list(map = map, first = first)
}

#' Groups of identical rows
#'
#' Finds the rows of a design matrix whose inputs (factor codes,
#' numeric values and scales) have the same bits, and therefore give the
#' same row in the model matrix and in every row-wise computation. Use it to
#' compute once per distinct row and expand, or to weight distinct rows by
#' their counts when averaging. The kernels in [design_algebra] use the same
#' groups (see the option `xtx.dedup`); they are found once per
#' object and cached.
#'
#' @param x A `design_matrix` built from a formula or plan.
#' @return A list: `map`, the id (1 to `U`) of each row's group, in order of
#'   first appearance; `first`, the index of the first row of each group;
#'   and `counts`, the number of rows in each group.
#' @examples
#' d <- data.frame(x = rep(1:3, 4), f = gl(2, 6))
#' g <- design_row_groups(design_matrix(~ x * f, d))
#' g$counts
#' @keywords internal
#' @noRd
design_row_groups <- function(x) {
    checkmate::assert_class(x, "xtx_matrix")
    J <- as_jac(x)
    if (!jac_kernel_ok(J))
        stop("`x` must be built from a formula or a plan.", call. = FALSE)
    g <- xtx_dedup_groups(J, force = TRUE)
    list(map = g$map, first = g$first, counts = tabulate(g$map, length(g$first)))
}

## Average number of entries per row over the first rows (a cost estimate),
## cached.
xtx_row_cost <- function(J) {
    cache <- J$cache
    if (!is.null(cache) && !is.null(cache$entries)) return(cache$entries)
    m <- do.call(.Call, c(list(C_xtx_jac_entries_mean), jac_args(J), list(256L)))
    if (!is.null(cache)) cache$entries <- m
    m
}

#' @rdname design_algebra
#' @param y For `%*%`, a numeric vector of length `ncol(x)` or a matrix with
#'   `ncol(x)` rows (the coefficients). For `crossprod()`, a numeric vector
#'   with one value per row or a matrix with `nrow(x)` rows (the weights),
#'   another `design_matrix` with the same rows, or `NULL` for `t(X) %*% X`.
#' @export
`%*%.xtx_matrix` <- function(x, y) {
    if (!inherits(x, "xtx_matrix"))
        stop("Only `design %*% beta` is supported, with the design first.", call. = FALSE)
    xtx_mult(x, y)
}

#' @rdname design_algebra
#' @param ... Not used.
#' @method crossprod xtx_matrix
#' @export
crossprod.xtx_matrix <- function(x, y = NULL, ...) {
    if (!inherits(x, "xtx_matrix"))
        stop("Only `crossprod(design, w)` and `crossprod(design, design)` are supported, with a design first.", call. = FALSE)
    if (is.null(y) || inherits(y, "xtx_matrix")) return(xtx_gram(x, y))
    xtx_crossprod(x, y)
}

## t(X) %*% Y for two design matrices (Y = X when `y` is NULL): a base
## matrix. When both are scaled rows of the same plain design, the scales
## become weights of a symmetric product, so that identical rows of the
## design can be grouped.
xtx_gram <- function(x, y = NULL) {
    n <- nrow(x)
    if (!is.null(y) && nrow(y) != n)
        stop("The design matrices must have the same number of rows.", call. = FALSE)
    nm <- list(colnames(x), colnames(if (is.null(y)) x else y))
    x0 <- x
    y0 <- y
    w <- NULL
    tx <- xtx_terms(x)
    if (is.null(y)) {
        if (length(tx) == 1L && !is.null(tx[[1L]]$scale) && xtx_safe_scale(tx[[1L]]$scale)) {
            w <- tx[[1L]]$scale^2
            x <- tx[[1L]]$design
        }
    } else {
        ty <- xtx_terms(y)
        if (length(tx) == 1L && length(ty) == 1L &&
            xtx_same_design(tx[[1L]]$design, ty[[1L]]$design) &&
            xtx_safe_scale(tx[[1L]]$scale) && xtx_safe_scale(ty[[1L]]$scale)) {
            w <- xtx_times(tx[[1L]]$scale, if (is.null(ty[[1L]]$scale)) 1 else ty[[1L]]$scale)
            x <- tx[[1L]]$design
            y <- NULL
        }
    }
    if (is.null(y) && (is.null(w) || xtx_all_finite(w)) && xtx_struct_ok(x)) {
        out <- xtx_struct_gram(x, if (is.null(w)) NULL
                                 else if (is.double(w) && length(w) == n && is.null(attributes(w))) w
                                 else rep_len(as.double(w), n))
        if (xtx_all_finite(out)) {
            dimnames(out) <- nm
            return(out)
        }
        ## an overflow in the structured sums: the row kernel on the
        ## original operands reproduces the dense result
        x <- x0
        y <- y0
        w <- NULL
    }
    A <- as_jac(x)
    B <- if (is.null(y)) NULL else as_jac(y)
    if (!jac_kernel_ok(A) || (!is.null(B) && !jac_kernel_ok(B))) {
        XA <- as.matrix(A)
        out <- if (is.null(B)) {
            if (is.null(w)) crossprod(XA) else crossprod(XA, w * XA)
        } else {
            crossprod(XA, as.matrix(B))
        }
        dimnames(out) <- nm
        return(out)
    }
    if (!is.null(w) && length(w) == 1L) w <- rep(w, n)
    first <- NULL
    if (is.null(B)) {
        g <- xtx_dedup_groups(A, cost = xtx_row_cost(A)^2)
        if (!is.null(g)) {
            w <- drop(rowsum(if (is.null(w)) rep(1, n) else w, g$map, reorder = TRUE))
            first <- g$first
        }
    }
    if (!is.null(w)) w <- as.double(w)
    Bargs <- if (is.null(B)) rep(list(NULL), 6L) else jac_args(B)
    out <- do.call(.Call, c(list(C_xtx_jac_gram), jac_args(A), Bargs, list(w, first, design_threads())))
    dimnames(out) <- nm
    out
}

## A plain design with row scales that the structured kernels can handle
## (moving the scale outside the product is safe), as list(design, scale),
## or NULL.
xtx_struct_single <- function(J) {
    if (!is.null(J$design2) || !xtx_struct_ok(J$design)) return(NULL)
    if (!xtx_safe_scale(J$scale)) return(NULL)
    list(design = J$design, scale = J$scale)
}

## Can row scales be combined or moved to the other factor of a product
## (s^2 * x * x for (s * x) * (s * x), say) without overflow or underflow
## that the product in its original order would not have? NULL is 1.
xtx_safe_scale <- function(s) {
    if (is.null(s)) return(TRUE)
    xtx_in_range(s, 1e-50, 1e50)
}

## Are all values finite, with non-zero absolute values within [lo, hi]?
## The structured kernels regroup products; with inputs within the bounds
## that xtx_struct_ok() enforces (1e-50 to 1e50 per factor), row scales
## within xtx_safe_scale() and the other operand of a product within
## [1e-100, 1e100], no intermediate leaves 1e+-300, so regrouping cannot
## overflow or underflow.
xtx_in_range <- function(x, lo, hi) {
    if (!is.double(x)) storage.mode(x) <- "double"
    .Call(C_xtx_in_range, x, lo, hi)
}

## Are all values of a numeric vector or matrix finite? (No allocation.)
xtx_all_finite <- function(x) {
    if (!is.double(x)) storage.mode(x) <- "double"
    .Call(C_xtx_all_finite, x)
}

## Are two plain designs the same object (or equal)?
xtx_same_design <- function(a, b) {
    if (!is.null(a$cache) && identical(a$cache, b$cache)) return(TRUE)
    identical(a, b)
}
