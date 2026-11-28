#' QR decomposition of a design matrix without forming the dense matrix
#'
#' `qr(x)` decomposes a `design_matrix` as `X P = Q R` without forming the
#' dense matrix `X` or the `n` by `p` matrix `Q`. It returns an object of
#' class `"xtx_qr"`, which the usual base R functions accept: `qr.coef()`,
#' `qr.fitted()`, `qr.resid()`, `qr.qty()`, `qr.qy()`, `qr.Q()`, `qr.R()` and
#' `solve()`. This gives least squares coefficients, fitted values and
#' residuals at the cost of a few passes over the rows.
#'
#' **Method.** `R` is the Cholesky factor of `X'X`, computed by
#' `crossprod(x)` without forming `X`. Columns are pivoted as in
#' `qr(as.matrix(x))` (LINPACK's limited pivoting, used by [lm()]): they
#' are taken in order, and a column whose residual norm is below `tol`
#' times its norm is moved to the end and counts as aliased. From `X'X`, a
#' residual norm is accurate only to about `sqrt(.Machine$double.eps)`
#' times the column norm, which cannot separate an aliased column from one
#' just above `tol`; the residual norms of such columns are computed from
#' the data (one product with `x` per block of columns). `Q` is
#' `X[, pivot] %*% solve(R)`, never formed. Coefficients solve the normal
#' equations and are then refined from the residuals `y - X b`, each step
#' costing one `x %*% b` and one `crossprod(x, r)`. Refinement recovers
#' much of the accuracy that the normal equations lose on badly conditioned
#' designs. When `X'X` is badly conditioned and numeric variables are far
#' from zero (`x + 1e7`), the values of `R` come from the `X'X` of the
#' design with those variables centered, when that is an exact linear
#' transform `X T` of the design (`R` is then `R~ T^-1`); the columns kept
#' and their order still come from `X'X`.
#'
#' **Agreement with `qr()`.** Results agree with those of
#' `qr(as.matrix(x))` up to floating point rounding, not bit for bit. On
#' badly conditioned designs the coefficients agree to fewer digits than
#' with a Householder QR. A column whose residual norm is within rounding
#' of `tol` can be aliased by one decomposition and not the other. Rows of `R` can
#' differ in sign (the diagonal here is positive), and so can the columns
#' of `Q` and the values of `qr.qty()`. For an aliased design, the block of
#' `R` for the aliased columns is zero.
#'
#' **What is not formed.** The complete `Q` (`n` by `n`) is never formed,
#' so:
#' * `qr.Q(q)` returns the `n` by `rank` matrix of the first `rank` columns
#'   of `Q`, which spans the columns of `X`; `complete = TRUE` is not
#'   supported.
#' * `qr.qty(q, y)` returns the first `rank` rows of `t(Q) %*% y`.
#' * `qr.qy(q, y)` takes `y` with `rank` rows and returns `Q[, 1:rank] %*% y`.
#' * `qr.fitted()` only supports `k = q$rank`.
#'
#' Non-finite entries in `x` are an error, as in [qr()]. `y` must be finite.
#' When `X'X` cannot represent the columns' scale (entries whose squares
#' overflow or underflow, around `1e154` and `1e-154`), `qr()` falls back
#' to the dense QR of `qr(as.matrix(x))`, kept in the object, and the
#' methods give its results.
#'
#' **Calling the methods.** `qr.coef()`, `qr.fitted()`, `qr.resid()`,
#' `qr.qty()`, `qr.qy()`, `qr.Q()` and `qr.R()` are S4 generics created by
#' this package (as the Matrix package does for its decompositions): they
#' dispatch when xtx is attached, or when called as
#' `xtx::qr.coef()`. Calls to `base::qr.coef()` (for example from
#' code in another package that imports base's function) do not dispatch,
#' and stop because the object is not a base QR decomposition. `qr()` and
#' `solve()` are S3 generics and always dispatch.
#'
#' @param x A `design_matrix`.
#' @param tol The tolerance for detecting aliased columns, as in [qr()].
#' @param ... Not used, except `y` for a remote design: a column of the
#'   design (a [design_column()]) whose products `X'y` are computed in the
#'   scan that computes `X'X`, and kept for `qr.coef(qr, y)` on that column.
#' @return An object of class `"xtx_qr"`: a list with `design` (`x`),
#'   `R` (the `p` by `p` upper triangular factor, columns in pivoted
#'   order), `rank`, `pivot` and `tol`. For a remote design, also the sums
#'   that its scan computed: `weight_sum` (the sum of the weights of the
#'   rows: their number for an unscaled design) and, with `y`,
#'   `response_ss` (the sum of squares of `y`). `qr.coef()` on a remote design returns the
#'   coefficients with an attribute `"moments"`: a list with `b` (the
#'   coefficients of the kept columns), `h` (their `X'y`) and `G` (their
#'   `X'X`), from which sums of squares follow without another pass over
#'   the rows (`y'y - sum(b * h)` is the residual sum of squares); these
#'   are those of the centered design when the factor is.
#' @examples
#' set.seed(1)
#' d <- data.frame(y = rnorm(100), x = rnorm(100), f = gl(4, 25))
#' D <- design_matrix(y ~ x * f, d)
#' q <- qr(D)
#' all.equal(qr.coef(q, d$y), coef(lm(y ~ x * f, d)))
#' all.equal(qr.fitted(q, d$y), unname(fitted(lm(y ~ x * f, d))))
#' @name design_qr
NULL

#' @rdname design_qr
#' @method qr xtx_matrix
#' @export
qr.xtx_matrix <- function(x, tol = 1e-07, ...) {
    checkmate::assert_number(tol, lower = 0, finite = TRUE)
    dots <- list(...)
    if (isTRUE(dots$LAPACK))
        stop("`LAPACK = TRUE` is not supported for design matrices.", call. = FALSE)
    G <- crossprod(x)
    storage.mode(G) <- "double"
    if (!all(is.finite(G))) {
        if (!all(finite_rows(x))) stop("NA/NaN/Inf in 'x'", call. = FALSE)
        ## finite entries whose squares overflow
        return(xtx_qr_dense(x, tol))
    }
    dg <- diag(G)
    ## squares near the limits of double precision lose the column's scale
    if (any(dg > 1e280 | (dg > 0 & dg < 1e-280))) return(xtx_qr_dense(x, tol))
    ch <- xtx_qr_chol(x, G, tol)
    if (is.null(ch)) return(xtx_qr_dense(x, tol))
    ch <- xtx_qr_center_memory(x, G, ch)
    R <- ch$R
    cn <- colnames(x)
    if (!is.null(cn)) colnames(R) <- cn[ch$pivot]
    out <- structure(list(design = x, R = R, rank = ch$rank, pivot = ch$pivot, tol = tol),
                     class = "xtx_qr")
    out$centered <- ch$centered
    out
}

## In memory, the factor of a badly conditioned X'X takes the values of the
## centered design, as remotely (xtx_qr_remote(), xtx_qr_centered()); the
## columns kept and their order still come from the uncentered X'X, as in
## lm(). The normal equations lose about kappa^2 * eps (kappa: the
## condition number of X with columns scaled to unit norm, see xtx_qr_coef()):
## refinement recovers the coefficients but not R, hence vcov(). Centering
## costs a shifted copy of the numeric variables and one more crossprod(),
## so it is done only when kappa^2 > 1e6, where that loss could exceed
## about 2e-10 (with a margin for the estimate of kappa from rcond());
## below, the factor is that of X'X, unchanged.
xtx_qr_center_memory <- function(x, G, ch) {
    k <- ch$rank
    if (k < 2L || !isTRUE(getOption("xtx.qr_center", TRUE))) return(ch)
    if (!(xtx_qr_kappa2(ch$R[seq_len(k), seq_len(k), drop = FALSE]) > 1e6)) return(ch)
    ctr <- xtx_qr_center_design(x)
    if (is.null(ctr)) return(ch)
    Gs <- crossprod(ctr$x)
    storage.mode(Gs) <- "double"
    if (!all(is.finite(Gs))) return(ch)
    xtx_qr_centered(ch, G, ctr$T, ctr$x, list(xtx = Gs))
}

## kappa^2 of the columns of a triangular factor scaled to unit norm (an
## estimate, from rcond()).
xtx_qr_kappa2 <- function(R11) {
    k <- ncol(R11)
    kap <- 1 / rcond(R11 / rep(sqrt(colSums(R11^2)), each = k), triangular = TRUE)
    kap^2
}

## An in-memory design (or one with rows scaled, as the weighted fits make
## it) with its numeric variables centered: list(x, T) with x = X T (see
## xtx_p_center()), or NULL. The centering of a plain design is cached in it,
## so that refits on reweighted copies (the steps of a GLM) center once;
## the cache holds the shifted numeric variables only. The sample that
## decides the shift does not change the session's random numbers.
xtx_qr_center_design <- function(x) {
    if (inherits(x, "xtx_remote")) return(NULL)
    sc <- NULL
    if (xtx_is_comb(x)) {
        if (!is.null(x$design2)) return(NULL)
        sc <- x$scale
        x <- x$design
    }
    if (is.null(x$vars) || !is.environment(x$cache)) return(NULL)
    key <- list(x$vars, x$setup, x$contrasts, colnames(x))
    hit <- x$cache$qr_center
    if (!is.null(hit) && identical(hit$key, key)) {
        ctr <- hit$value
    } else {
        ctr <- xtx_qr_keep_seed(tryCatch(xtx_p_center(x), error = function(e) NULL))
        x$cache$qr_center <- list(key = key, value = ctr)
    }
    if (is.null(ctr)) return(NULL)
    if (!is.null(sc)) ctr$x <- xtx_scale(ctr$x, sc)
    ctr
}

## The value of expr, leaving .Random.seed as it was.
xtx_qr_keep_seed <- function(expr) {
    env <- globalenv()
    seed <- get0(".Random.seed", envir = env, inherits = FALSE)
    on.exit(if (is.null(seed)) {
        if (exists(".Random.seed", envir = env, inherits = FALSE)) rm(".Random.seed", envir = env)
    } else {
        assign(".Random.seed", seed, envir = env)
    })
    expr
}

## The same factorization for a remote design, from the primitives (see
## primitives.R): X'X in the database, and the residual norms of suspect
## columns from products with the data. There is no dense fallback.
#' @method qr xtx_remote
#' @export
qr.xtx_remote <- function(x, tol = 1e-07, ..., y = NULL) {
    if (!is.null(y) && !xtx_is_column(y)) stop("`y` must be a column of the design (a `design_column()`).", call. = FALSE)
    xtx_qr_remote(x, tol, y)
}

## (`Y`, remote columns: X'Y comes from the scan that computes X'X, and is
## kept for xtx_qr_coef(), see xtx_remote_gram_xty())
##
## The columns kept and their order come from X'X (lm()'s rule, as in
## memory), but the values of the factor of the kept columns come from the
## X'X of the design with numeric variables shifted by their mean in a
## sample, when that design is X~ = X T with T unit upper triangular
## (xtx_center()): R = R~ T^-1. A numeric variable far from zero (x + 1e7)
## makes X'X nearly singular, while the shifted one is well conditioned; the
## database sums X'X with another rounding than the in-memory kernels, so
## without the shift the factor would lose the accuracy that in-memory
## designs keep. The sums are then those of the shifted design only (one
## scan), and X'X = T^-T (X~'X~) T^-1 follows from them.
xtx_qr_remote <- function(x, tol, Y = NULL) {
    checkmate::assert_number(tol, lower = 0, finite = TRUE)
    ctr <- if (isTRUE(getOption("xtx.remote_center", TRUE))) xtx_center(x)
    mc <- xs <- NULL
    if (!is.null(ctr)) {
        xs <- xtx_shifted(x, ctr$shift)
        ## (not compressed: the design of the patterns would not be shifted)
        op <- options(xtx.remote_compress = FALSE)
        mc <- tryCatch(xtx_remote_gram_xty(xs, Y), finally = options(op))
        if (!all(is.finite(mc$xtx))) mc <- NULL
    }
    if (is.null(mc)) {
        mo <- xtx_remote_gram_xty(x, Y)
    } else {
        Ti <- backsolve(ctr$T, diag(ncol(ctr$T)))
        G <- crossprod(Ti, mc$xtx %*% Ti)
        mo <- list(xtx = (G + t(G)) / 2, xty = if (!is.null(mc$xty)) crossprod(Ti, mc$xty))
        dimnames(mo$xtx) <- dimnames(mc$xtx)
    }
    G <- mo$xtx
    storage.mode(G) <- "double"
    if (!all(is.finite(G))) stop("NA/NaN/Inf in 'x'", call. = FALSE)
    ch <- xtx_qr_chol(x, G, tol)
    if (is.null(ch)) stop("The squares of some entries of the design underflow.", call. = FALSE)
    if (!is.null(mc)) ch <- xtx_qr_centered(ch, G, ctr$T, xs, mc, Y)
    R <- ch$R
    cn <- colnames(x)
    if (!is.null(cn)) colnames(R) <- cn[ch$pivot]
    structure(list(design = x, R = R, rank = ch$rank, pivot = ch$pivot, tol = tol, G = unname(G),
                   response_ss = if (is.null(mc)) mo$yty else mc$yty,
                   weight_sum = if (is.null(mc)) mo$sumw else mc$sumw,
                   centered = ch$centered, xty = mo$xty, xty_of = if (!is.null(mo$xty)) xtx_column_key(Y)),
              class = "xtx_qr")
}

## The factorization of a design for a least squares fit of the columns Y:
## qr(x), which remotely also computes X'Y in the same scan.
xtx_qr_for <- function(x, tol, Y) {
    if (inherits(x, "xtx_remote") && xtx_is_column(Y)) xtx_qr_remote(x, tol, Y) else qr(x, tol = tol)
}

## X'Y for the columns Y, from the scan of the factorization when it was
## made for them (xtx_qr_for()), else computed.
xtx_qr_xty <- function(q, D, Y) {
    ## (the response by its identity, xtx_column_key(): the same SQL text
    ## over another design or stored relation is another response)
    if (!is.null(q$xty) && xtx_is_column(Y) && identical(q$xty_of, xtx_column_key(Y))) return(q$xty)
    xtx_p_crossprod(D, Y)
}

## The factor of a design with the values of the shifted design (see
## xtx_qr_remote(); in memory, xtx_qr_center_memory()): `xs` the shifted
## design, `mc` its X'X and (remotely) X'Y, `Tm` with X~ = X T. Unchanged
## when the kept columns do not allow it (a kept column of X~ uses an
## aliased column of X).
xtx_qr_centered <- function(ch, G, Tm, xs, mc, Y = NULL) {
    if (ch$rank == 0L) return(ch)
    k <- ch$rank
    p <- length(ch$pivot)
    keep <- ch$pivot[seq_len(k)]
    al <- ch$pivot[seq_len(p) > k]
    ## The kept columns are in increasing order, so that T[keep, keep] is
    ## upper triangular and R~ T^-1 a triangular factor: lm()'s pivoting
    ## (xtx_chol_pivot()) takes columns in order and moves only the others
    ## to the end. Anything else would be a bug, not a design to fall back
    ## on silently.
    if (is.unsorted(keep)) stop("internal error: the kept columns are not in order", call. = FALSE)
    ## (not when a kept column of X~ uses an aliased column of X)
    if (length(al) && any(Tm[al, keep] != 0)) return(ch)
    Rt <- tryCatch(chol(unname(mc$xtx[keep, keep, drop = FALSE])), error = function(e) NULL)
    if (is.null(Rt)) return(ch)
    R11 <- Rt %*% backsolve(Tm[keep, keep, drop = FALSE], diag(k))
    R <- ch$R
    R[seq_len(k), seq_len(k)] <- R11
    if (length(al)) {
        R[seq_len(k), k + seq_along(al)] <- backsolve(R11, unname(G[keep, al, drop = FALSE]), transpose = TRUE)
    }
    ch$R <- R
    ch$centered <- list(design = xs, T = Tm[keep, keep, drop = FALSE], Rt = Rt, xty = mc$xty,
                        xty_of = if (!is.null(mc$xty)) xtx_column_key(Y))
    ch
}

## The remote design with numeric variables shifted (their SQL).
xtx_shifted <- function(x, shift) {
    for (nm in names(shift))
        x$vars[[nm]]$sql <- sprintf("(%s - %s)", x$vars[[nm]]$sql, xtx_sql_double(shift[[nm]]))
    x
}

## Centering of the numeric variables of a remote design: list(shift (named
## by variable), T) with X(shifted) = X T, T unit upper triangular, or NULL
## (nothing to center, or no such T); from a sample of its rows.
xtx_center <- function(D) {
    ## (once for a design: its scaled copies share the cache)
    if (is.null(D$cache$center)) D$cache$center <- list(xtx_center_sample(D))
    D$cache$center[[1L]]
}

xtx_center_sample <- function(D) {
    plan <- D$plan
    sample <- xtx_remote_sample(D)
    if (!nrow(sample)) return(NULL)
    mf <- tryCatch(design_frame(plan$terms, sample, xlev = plan$xlevels, na.action = stats::na.omit),
                   error = function(e) NULL)
    if (is.null(mf) || !nrow(mf)) return(NULL)
    if (D$p > 2000L) return(NULL)
    xtx_center_shift(mf, mf, function(fr) xtx_frame_design(D, fr), D$p)
}

## The shift of the numeric variables of a model frame whose mean (in the
## rows of `stats_frame`) is large compared with their spread, and T with
## X(shifted) = X T unit upper triangular, checked on prototype rows drawn
## from `proto_frame` (levels and logical values at random, numeric values
## near zero, which keeps the prototype design well conditioned), the
## designs built by build(frame) with p columns: list(shift, T), or NULL.
## The columns of the frames are taken as stored (I(x^2), log(x): columns,
## shifted on their own when far from zero), in the prototypes as in the
## shifted design (xtx_shifted() remotely, xtx_p_center() in memory).
## (For remote designs, xtx_center(); in memory, xtx_center_memory().)
xtx_center_shift <- function(stats_frame, proto_frame, build, p) {
    tt <- attr(stats_frame, "terms")
    resp <- if (!is.null(tt) && attr(tt, "response") > 0L) xtx_terms_vars(tt)$names[attr(tt, "response")]
    cand <- setdiff(names(stats_frame), resp)
    num <- cand[vapply(cand, function(nm) identical(stats::.MFclass(stats_frame[[nm]]), "numeric"), NA)]
    if (!length(num)) return(NULL)
    mu <- vapply(num, function(nm) mean(stats_frame[[nm]]), 0)
    sdv <- vapply(num, function(nm) stats::sd(stats_frame[[nm]]), 0)
    sdv[!is.finite(sdv)] <- 0
    ## only variables whose mean is large compared with their spread
    big <- is.finite(mu) & abs(mu) > sdv
    if (!any(big)) return(NULL)
    shift <- mu[big]
    m <- max(4L * p, 200L)
    proto <- proto_frame[sample.int(nrow(proto_frame), m, replace = TRUE), , drop = FALSE]
    for (nm in names(proto)) {
        x <- proto[[nm]]
        if (is.factor(x)) {
            proto[[nm]] <- factor(sample(levels(x), m, TRUE), levels = levels(x))
            attr(proto[[nm]], "contrasts") <- attr(x, "contrasts")
        } else if (is.logical(x)) {
            proto[[nm]] <- sample(c(FALSE, TRUE), m, TRUE)
        } else if (nm %in% num) {
            proto[[nm]] <- stats::rnorm(m)
        }
    }
    ## The designs are built from the columns of the frames as they are
    ## stored (a model frame: its "terms"), as the shifted design is: a
    ## variable such as I(x^2) or log(x) is a column of the frame, shifted
    ## on its own, not rebuilt from a shifted x (design.matrix#r6g7).
    if (!is.null(tt)) attr(proto, "terms") <- tt
    shifted <- proto
    for (nm in names(shift)) shifted[[nm]] <- proto[[nm]] - shift[[nm]]
    X0 <- tryCatch(as.matrix(build(proto)), error = function(e) NULL)
    X1 <- tryCatch(as.matrix(build(shifted)), error = function(e) NULL)
    if (is.null(X0) || is.null(X1) || ncol(X0) != p || ncol(X1) != p) return(NULL)
    q <- qr(X0)
    if (q$rank < p) return(NULL)
    Tm <- unname(qr.coef(q, X1))
    if (!all(is.finite(Tm)) || max(abs(X0 %*% Tm - X1)) > 1e-8 * max(1, abs(X1))) return(NULL)
    ## (relative to the scale of each column of T: the shift of a product
    ## of variables far from zero, such as x * z with x and z near 1e7,
    ## puts entries near 1e14 and near 1 in T)
    ts <- matrix(pmax(1, apply(abs(Tm), 2L, max)), p, p, byrow = TRUE)
    if (any(abs(diag(Tm) - 1) > 1e-8 * diag(ts)) || any(abs(Tm[lower.tri(Tm)]) > 1e-8 * ts[lower.tri(ts)]))
        return(NULL)
    Tm[lower.tri(Tm)] <- 0
    Tm[abs(Tm) < 1e-10 * ts] <- 0
    diag(Tm) <- 1
    list(shift = shift, T = xtx_center_exact(Tm, shift, proto, build, q, attr(X0, "assign"), tt))
}

## T of xtx_center_shift() to rounding in every entry. The least squares
## solve above finds each column of T to about 1e-16 of its largest entry
## only: with x and z near 1e7, the column of x:z has an entry near 1e14
## and entries near 1e7, which come out with nine digits, and coefficients
## mapped back with them (b = T b~) shift the fitted values. But T is
## known in closed form up to constants: the column of a term with shifted
## variables A is the product of (v_a - s_a) over A, so its entry on a
## column of a term with variables B is a constant times the product of
## the s_a, a in A but not in B. The constants are T at unit shifts (a
## well scaled solve, with the factorization `q` of the prototype design),
## the products come from the terms (`assign`: the term of each column).
## Returns `Tm` itself when the terms do not describe the columns, or when
## the two disagree beyond the accuracy of `Tm`.
xtx_center_exact <- function(Tm, shift, proto, build, q, assign, tt) {
    p <- ncol(Tm)
    fac <- if (!is.null(tt)) attr(tt, "factors")
    if (is.null(assign) || length(assign) != p || !is.matrix(fac) || max(assign) > ncol(fac)) return(Tm)
    unit <- proto
    for (nm in names(shift)) unit[[nm]] <- proto[[nm]] - 1
    X1 <- tryCatch(as.matrix(build(unit)), error = function(e) NULL)
    if (is.null(X1) || ncol(X1) != p) return(Tm)
    K <- unname(qr.coef(q, X1))
    if (!all(is.finite(K))) return(Tm)
    ## the shifted variables of each term (0: the intercept, none), and the
    ## product of the shifts of those a term has and another lacks
    vars <- lapply(c(0L, seq_len(ncol(fac))), function(t) {
        if (t == 0L) character(0) else intersect(rownames(fac)[fac[, t] > 0], names(shift))
    })
    nt <- length(vars)
    St <- matrix(1, nt, nt)
    for (j in seq_len(nt)) for (i in seq_len(nt)) St[i, j] <- prod(shift[setdiff(vars[[j]], vars[[i]])])
    K[abs(K) < 1e-10] <- 0
    Tx <- K * St[assign + 1L, assign + 1L, drop = FALSE]
    Tx[lower.tri(Tx)] <- 0
    diag(Tx) <- 1
    ts <- matrix(pmax(1, apply(abs(Tm), 2L, max)), p, p, byrow = TRUE)
    if (!all(is.finite(Tx)) || any(abs(Tx - Tm) > 1e-6 * ts)) return(Tm)
    Tx
}

## Centering of the numeric variables of an in-memory design (see
## xtx_p_center()), from its terms mt and model frame mf: list(shift, T), with
## the design of the shifted frame X~ = X T and T unit upper triangular, or
## NULL. As xtx_center() does for remote designs: only variables whose mean
## is large compared with their spread, and only when such a T exists
## (checked on prototype rows built, like the shifted design, from the
## frame's columns as stored).
xtx_center_memory <- function(mt, mf, contrasts) {
    if (!isTRUE(getOption("xtx.fe_center", TRUE))) return(NULL)
    if (attr(mt, "intercept") == 0L) return(NULL)
    ## the decision and the shift from a systematic sample of rows (any
    ## shift is exact; a shift near the mean keeps the sums accurate)
    tv <- xtx_terms_vars(mt)$names
    nr <- nrow(mf)
    idx <- if (nr > 1e5) unique(round(seq(1, nr, length.out = 1e5))) else seq_len(nr)
    st <- mf[idx, intersect(tv, names(mf)), drop = FALSE]
    attr(st, "terms") <- mt
    m0 <- mf[seq_len(min(nrow(mf), 2000L)), , drop = FALSE]
    attr(m0, "terms") <- mt
    p <- tryCatch(ncol(design_matrix(mt, m0, contrasts.arg = contrasts)), error = function(e) NA)
    if (is.na(p) || p > 2000L) return(NULL)
    xtx_center_shift(st, m0, function(fr) design_matrix(mt, fr, contrasts.arg = contrasts), p)
}

## The in-memory design of a model frame of the plan's terms.
xtx_frame_design <- function(D, mf) {
    plan <- D$plan
    d <- design_matrix(plan$terms, mf, contrasts.arg = plan$contrasts)
    if (isTRUE(plan$drop_intercept)) d <- xtx_drop_intercept(d)
    d
}


## A design whose X'X loses its scale: the dense QR of qr(), kept in the
## object, so that every method gives the dense results.
xtx_qr_dense <- function(x, tol) {
    q <- base::qr(as.matrix(x), tol = tol)
    p <- ncol(x)
    R <- matrix(0, p, p)
    Rd <- base::qr.R(q)
    R[seq_len(nrow(Rd)), ] <- Rd
    cn <- colnames(x)
    if (!is.null(cn)) colnames(R) <- cn[q$pivot]
    structure(list(design = x, R = R, rank = q$rank, pivot = q$pivot, tol = tol, dense = q),
              class = "xtx_qr")
}

#' Pivoted Cholesky factor with limited pivoting
#'
#' `pivoted_chol()` factors a symmetric positive semi-definite matrix `G`
#' (such as `X'X`) with the limited pivoting of [stats::lm()]: columns are
#' taken in order, and a column whose squared residual norm `d` after the
#' columns taken before it is not above `tol * scale[j]` is moved to the
#' end. This is the kernel of [qr()][design_qr] of a design, and of fits
#' with fixed effects.
#'
#' From `G` alone, `d` is accurate only to rounding, so columns near the
#' threshold cannot be judged. With `band`, columns with `lo <= d <= hi`
#' are *suspects*: `recheck(ch, cols)` is called with the current
#' factorization and the suspect columns, and returns their squared
#' residual norms computed from the data. The factorization is repeated
#' until every suspect was judged against the columns actually taken before
#' it.
#'
#' @param G A symmetric numeric matrix.
#' @param scale One positive number per column: the scale of its squared
#'   norm (such as `diag(G)`).
#' @param tol The tolerance, relative to `scale`.
#' @param ties `"keep"`: a column is taken when `d >= tol * scale[j]` (the
#'   rule of `lm()`); `"drop"`: when `d > tol * scale[j]`.
#' @param band `NULL`, or a list with `lo` and `hi`, one value per column:
#'   the range of `d` in which a column is a suspect.
#' @param recheck With `band`, a function of the factorization and the
#'   suspect columns that returns a list with `dres` (their squared residual
#'   norms) and `taken` (which of them pass the rule), or `NULL` to give up.
#' @param lup Columns are tested until this many are taken.
#' @param dres Squared residual norms known in advance (`NA`: not known).
#' @return A list with `R` (the upper triangular factor, columns in pivoted
#'   order), `rank`, `pivot` and `suspects`; `NULL` when `recheck()` gave
#'   up.
#' @examples
#' X <- cbind(1, x = 1:5, x2 = 2 * (1:5), z = c(1, 0, 2, 1, 3))
#' G <- crossprod(X)
#' ch <- pivoted_chol(G, diag(G), 1e-14)
#' ch$rank
#' ch$pivot
#' @keywords internal
#' @export
pivoted_chol <- function(G, scale, tol, ties = c("keep", "drop"), band = NULL, recheck = NULL,
                       lup = ncol(G), dres = rep(NA_real_, ncol(G))) {
    ties <- match.arg(ties)
    p <- ncol(G)
    lo <- if (is.null(band)) double(0) else as.double(band$lo)
    hi <- if (is.null(band)) double(0) else as.double(band$hi)
    repeat {
        ch <- .Call(C_xtx_chol_pivot, G, as.double(scale), as.double(tol), ties == "keep", lo, hi,
                    as.integer(lup), as.double(dres))
        names(ch) <- c("R", "rank", "pivot", "suspects")
        new <- ch$suspects[is.na(dres[ch$suspects])]
        if (length(new) == 0L) return(ch)
        v <- recheck(ch, new)
        if (is.null(v)) return(NULL)
        dres[new] <- v$dres
        taken <- new[v$taken]
        if (length(taken) == 0L) return(ch)
        ## a suspect turned out to be independent: every later column,
        ## suspects included, was judged without it
        dres[seq_len(p) > min(taken)] <- NA_real_
    }
}

## Pivoted Cholesky factor of G = X'X with lm()'s rule (pivoted_chol()).
## Columns whose residual norm is too small to judge from G, or whose norm
## in G is zero, are suspects: their residual norms are computed from the
## data (xtx_qr_resid2()). Returns NULL when a column with a zero norm in G is
## not zero in the data (its square underflowed).
xtx_qr_chol <- function(x, G, tol) {
    p <- ncol(G)
    ## columns taken from G alone are clearly above tol
    susp <- max(1e-8, 100 * tol^2)
    dg <- diag(G)
    dres <- rep(NA_real_, p)
    ## A zero norm in G is a zero column unless the column's non-zero
    ## entries can be so small that their squares underflow: when the
    ## design's smallest possible non-zero entry rules that out, zero
    ## columns need no check against the data.
    if (any(dg == 0) && xtx_p_min_nonzero(x)^2 > 1e-300) dres[dg == 0] <- 0
    recheck <- function(ch, new) {
        v <- xtx_qr_resid2(x, G, ch, new)
        if (any(v$nonzero & !(dg[new] > 0))) return(NULL)
        ## the kernel takes a column when dres >= tol^2 * G[j, j], that is,
        ## when its residual norm is at least tol times its norm (as
        ## dqrdc2, a tie with tol keeps the column)
        list(dres = ifelse(v$nonzero, v$ratio^2 * dg[new], 0),
             taken = v$nonzero & v$ratio > 0 & v$ratio >= tol)
    }
    ## (every column with d <= susp * G[j, j], or a zero norm in G, is a
    ## suspect)
    pivoted_chol(G, dg, tol^2, "keep", band = list(lo = rep(-Inf, p), hi = susp * dg), recheck = recheck,
               lup = min(nrow(x), p), dres = dres)
}

## For columns `cols` of x: the norm of the residual after projection on
## the columns taken before each of them, relative to the column's norm,
## both computed from the data (one refinement step) and scaled so that
## neither overflows nor underflows; and whether the column has a non-zero
## entry. In blocks of columns.
xtx_qr_resid2 <- function(x, G, ch, cols) {
    acc <- ch$pivot[seq_len(ch$rank)]
    R <- ch$R
    p <- ncol(G)
    ratio <- numeric(length(cols))
    nonzero <- logical(length(cols))
    for (blk in split(seq_along(cols), (seq_along(cols) - 1L) %/% 32L)) {
        k <- length(blk)
        M <- matrix(0, p, k)
        E <- matrix(0, p, k)
        A <- vector("list", k)
        for (u in seq_len(k)) {
            j <- cols[blk[u]]
            a <- sum(acc < j)
            A[[u]] <- acc[seq_len(a)]
            Ra <- R[seq_len(a), seq_len(a), drop = FALSE]
            cf <- xtx_qr_bs(Ra, xtx_qr_bs(Ra, G[A[[u]], j, drop = FALSE], transpose = TRUE))
            M[j, u] <- 1
            M[A[[u]], u] <- -cf
            E[j, u] <- 1
        }
        res <- xtx_p_mult(x, M)
        h <- crossprod(x, res)
        for (u in seq_len(k)) {
            a <- length(A[[u]])
            Ra <- R[seq_len(a), seq_len(a), drop = FALSE]
            M[A[[u]], u] <- M[A[[u]], u] -
                xtx_qr_bs(Ra, xtx_qr_bs(Ra, h[A[[u]], u, drop = FALSE], transpose = TRUE))
        }
        ## a column's norm is sqrt(G[j, j]), accurate to rounding, unless
        ## G[j, j] is zero (possibly underflow): then it comes from the data
        dz <- diag(G)[cols[blk]] == 0
        nrm <- xtx_p_colnorms(xtx_p_mult(x, cbind(M, E[, dz, drop = FALSE])))
        res_n <- nrm[seq_len(k)]
        col_n <- sqrt(diag(G)[cols[blk]])
        col_n[dz] <- nrm[k + seq_len(sum(dz))]
        nonzero[blk] <- col_n > 0
        ratio[blk] <- ifelse(col_n > 0, res_n / col_n, 0)
    }
    list(ratio = ratio, nonzero = nonzero)
}

#' @export
print.xtx_qr <- function(x, ...) {
    cat(sprintf("<design_qr: %d x %d design, rank %d>\n",
                nrow(x$design), ncol(x$design), x$rank))
    invisible(x)
}

## y as a finite double matrix with `rows` rows; remembers whether it was a
## vector.
xtx_qr_y <- function(y, rows, what = "y") {
    vec <- !is.matrix(y)
    Y <- if (vec) matrix(y, ncol = 1L) else y
    if (!is.numeric(Y) && !is.logical(Y))
        stop(sprintf("`%s` must be numeric.", what), call. = FALSE)
    if (nrow(Y) != rows)
        stop(sprintf("`%s` must have %d rows.", what, rows), call. = FALSE)
    if (!all(is.finite(Y)))
        stop(sprintf("`%s` must be finite.", what), call. = FALSE)
    storage.mode(Y) <- "double"
    list(Y = unname(Y), vec = vec,
         rows = if (vec) names(y) else rownames(y),
         cols = if (vec) NULL else colnames(y))
}

## R11 (the factor of the first `rank` pivoted columns) and triangular
## solves with it, which also work at rank 0.
xtx_qr_R11 <- function(qr) {
    k <- qr$rank
    unname(qr$R[seq_len(k), seq_len(k), drop = FALSE])
}

xtx_qr_bs <- function(R11, rhs, transpose = FALSE) {
    if (nrow(R11) == 0L) return(matrix(0, 0L, ncol(rhs)))
    backsolve(R11, rhs, transpose = transpose)
}

## Solve R11' R11 b = rhs.
xtx_qr_solve <- function(qr, rhs) {
    R11 <- xtx_qr_R11(qr)
    xtx_qr_bs(R11, xtx_qr_bs(R11, rhs, transpose = TRUE))
}

## Coefficients of the full design (zero for aliased columns, original
## column order) from those of the first `rank` pivoted columns.
xtx_qr_full <- function(qr, b, fill = 0) {
    p <- ncol(qr$design)
    out <- matrix(fill, p, ncol(b))
    out[qr$pivot[seq_len(qr$rank)], ] <- b
    out
}

## Least squares coefficients of the first `rank` pivoted columns, with
## iterative refinement on the residuals. Each step costs one product and
## one crossproduct with the design. Refinement stops once the correction
## is at the level of rounding, or stops shrinking (the attainable accuracy
## is reached), and after at most `steps` steps.
xtx_qr_coef <- function(qr, Y, steps = 3L) {
    k <- qr$rank
    D <- qr$design
    keep <- qr$pivot[seq_len(k)]
    if (!is.null(qr$dense))
        return(base::qr.coef(qr$dense, as.matrix(Y))[keep, , drop = FALSE])
    if (!is.null(qr$centered)) {
        ## a remote design with numeric variables shifted (xtx_qr_centered()):
        ## X~'Y on the shifted design, solved with its well-conditioned
        ## factor, and mapped back (beta = T beta~)
        cc <- qr$centered
        op <- options(xtx.remote_compress = FALSE)
        on.exit(options(op))
        ht <- xtx_qr_xty(cc, cc$design, Y)[keep, , drop = FALSE]
        bt <- backsolve(cc$Rt, backsolve(cc$Rt, ht, transpose = TRUE))
        if (is.null(qr$G)) {
            ## in memory: refined from the residuals of the shifted design,
            ## as below
            if (xtx_qr_kappa2(cc$Rt) > 64) bt <- xtx_qr_refine(cc$design, cc$Rt, keep, bt, Y, steps = steps)
            return(cc$T %*% bt)
        }
        b <- cc$T %*% bt
        ## the moments of the solution, for sums of squares without a pass
        ## over the rows (see `qr.coef()` in ?design_qr)
        attr(b, "moments") <- list(b = bt, h = ht, G = crossprod(cc$Rt))
        return(b)
    }
    h <- xtx_qr_xty(qr, D, Y)
    b <- xtx_qr_solve(qr, h[keep, , drop = FALSE])
    if (k == 0L) return(b)
    ## The normal equations lose about kappa^2 * eps of accuracy, where
    ## kappa is the condition number of X with columns scaled to unit norm
    ## (Cholesky is insensitive to column scaling). When that is already at
    ## the level of rounding, refinement cannot help.
    R11 <- xtx_qr_R11(qr)
    if (xtx_qr_kappa2(R11) > 64) b <- xtx_qr_refine(D, R11, keep, b, Y, qr$G, h, steps)
    if (!is.null(qr$G)) attr(b, "moments") <- list(b = b, h = h[keep, , drop = FALSE], G = qr$G[keep, keep, drop = FALSE])
    b
}

## Iterative refinement of b, the solution of R11'R11 b = X'Y for the
## columns `keep` of the design D (R11 their factor). Stops once the
## correction is at the level of rounding, or stops shrinking.
xtx_qr_refine <- function(D, R11, keep, b, Y, G = NULL, h = NULL, steps = 3L) {
    p <- ncol(D)
    full <- function(b) {
        out <- matrix(0, p, ncol(b))
        out[keep, ] <- b
        out
    }
    ## X'(Y - X b): from the data in memory; for a remote design, from the
    ## moments already computed (X'Y - X'X b), which refines the solution
    ## of the normal equations without another pass
    residual <- function(b) {
        if (!is.null(G)) h - G %*% full(b) else xtx_p_crossprod(D, Y - xtx_p_mult(D, full(b)))
    }
    last <- Inf
    for (it in seq_len(steps)) {
        db <- xtx_qr_bs(R11, xtx_qr_bs(R11, residual(b)[keep, , drop = FALSE], transpose = TRUE))
        b <- b + db
        size <- max(abs(db)) / max(abs(b), .Machine$double.xmin)
        if (!(size > 64 * .Machine$double.eps) || size > 0.25 * last) break
        last <- size
    }
    b
}

## The result as a vector or a matrix like `y`, with row names `rows` and
## (for a matrix) the column names of `y`.
xtx_qr_out <- function(out, yy, rows = yy$rows) {
    if (yy$vec) {
        out <- out[, 1L]
        names(out) <- rows
    } else {
        dimnames(out) <- if (is.null(rows) && is.null(yy$cols)) NULL else list(rows, yy$cols)
    }
    out
}

## A row-wise result of the methods (the columns `cols`, one or several):
## a matrix or vector like `y` in memory (xtx_qr_out()), a table with columns
## `prefix` (or prefix1, prefix2, ...) remotely (xtx_p_output()).
xtx_qr_result <- function(qr, cols, prefix, yy, rows = yy$rows) {
    k <- xtx_ncols(cols)
    nm <- if (k == 1L) prefix else paste0(prefix, seq_len(k))
    xtx_p_output(qr$design, stats::setNames(lapply(seq_len(k), function(j) xtx_c_col(cols, j)), nm),
               memory = function(d) xtx_qr_out(unname(as.matrix(d)), yy, rows))
}

## The fitted values of y (a column, one or several).
xtx_qr_fitted <- function(qr, Y) {
    if (!is.null(qr$dense)) return(unname(base::qr.fitted(qr$dense, Y)))
    xtx_p_mult(qr$design, xtx_qr_full(qr, xtx_qr_coef(qr, Y)))
}

#' @importFrom methods setOldClass
NULL
setOldClass("xtx_qr")

#' @rdname design_qr
#' @param qr A `"xtx_qr"` object, from `qr()` on a `design_matrix`.
#' @param y A numeric vector or matrix: the response for `qr.coef()`,
#'   `qr.fitted()`, `qr.resid()` and `qr.qty()` (`nrow(x)` rows), or the
#'   values to multiply by `Q` for `qr.qy()` (`rank` rows).
#' @param k The rank for `qr.fitted()`; only `qr$rank` is supported.
#' @param complete,Dvec Not supported; only the defaults are accepted.
#' @param a,b For `solve()`: a `"xtx_qr"` object and the right-hand side,
#'   as in [solve.qr()].
#' @exportMethod qr.coef
setMethod("qr.coef", "xtx_qr", function(qr, y) {
    yy <- xtx_p_y(qr$design, y)
    b <- xtx_qr_coef(qr, yy$Y)
    out <- xtx_qr_full(qr, b, fill = NA_real_)
    out <- xtx_qr_out(out, yy, colnames(qr$design))
    ## (the moments of the solution, when the factorization has them)
    if (!is.null(attr(b, "moments"))) attr(out, "moments") <- attr(b, "moments")
    out
})

#' @rdname design_qr
#' @exportMethod qr.fitted
setMethod("qr.fitted", "xtx_qr", function(qr, y, k = qr$rank) {
    if (!identical(as.integer(k), as.integer(qr$rank)))
        stop("Only `k = qr$rank` is supported for design matrices.", call. = FALSE)
    yy <- xtx_p_y(qr$design, y)
    xtx_qr_result(qr, xtx_qr_fitted(qr, yy$Y), ".fitted", yy)
})

#' @rdname design_qr
#' @exportMethod qr.resid
setMethod("qr.resid", "xtx_qr", function(qr, y) {
    yy <- xtx_p_y(qr$design, y)
    xtx_qr_result(qr, yy$Y - xtx_qr_fitted(qr, yy$Y), ".resid", yy)
})

#' @rdname design_qr
#' @exportMethod qr.qty
setMethod("qr.qty", "xtx_qr", function(qr, y) {
    k <- qr$rank
    keep <- qr$pivot[seq_len(k)]
    yy <- xtx_p_y(qr$design, y)
    if (!is.null(qr$dense))
        return(xtx_qr_out(base::qr.qty(qr$dense, yy$Y)[seq_len(k), , drop = FALSE], yy,
                        yy$rows[seq_len(k)]))
    out <- xtx_qr_bs(xtx_qr_R11(qr), xtx_p_crossprod(qr$design, yy$Y)[keep, , drop = FALSE],
                   transpose = TRUE)
    xtx_qr_out(out, yy, yy$rows[seq_len(k)])
})

#' @rdname design_qr
#' @exportMethod qr.qy
setMethod("qr.qy", "xtx_qr", function(qr, y) {
    k <- qr$rank
    yy <- xtx_qr_y(y, k)
    if (!is.null(qr$dense)) {
        Z <- matrix(0, nrow(qr$design), ncol(yy$Y))
        Z[seq_len(k), ] <- yy$Y
        return(xtx_qr_out(unname(base::qr.qy(qr$dense, Z)), yy, NULL))
    }
    b <- xtx_qr_bs(xtx_qr_R11(qr), yy$Y)
    xtx_qr_result(qr, xtx_p_mult(qr$design, xtx_qr_full(qr, b)), ".qy", yy, NULL)
})

#' @rdname design_qr
#' @exportMethod qr.Q
setMethod("qr.Q", "xtx_qr", function(qr, complete = FALSE, Dvec) {
    if (!isFALSE(complete))
        stop("`complete = TRUE` is not supported for design matrices.", call. = FALSE)
    if (!missing(Dvec))
        stop("`Dvec` is not supported for design matrices.", call. = FALSE)
    k <- qr$rank
    if (!is.null(qr$dense)) return(unname(base::qr.Q(qr$dense)[, seq_len(k), drop = FALSE]))
    xtx_qr_result(qr, xtx_p_mult(qr$design, xtx_qr_full(qr, xtx_qr_bs(xtx_qr_R11(qr), diag(1, k)))), ".Q",
                list(vec = FALSE, cols = NULL), NULL)
})

#' @rdname design_qr
#' @exportMethod qr.R
setMethod("qr.R", "xtx_qr", function(qr, complete = FALSE) {
    if (!isFALSE(complete))
        stop("`complete = TRUE` is not supported for design matrices.", call. = FALSE)
    qr$R
})

#' @rdname design_qr
#' @method solve xtx_qr
#' @export
solve.xtx_qr <- function(a, b, ...) {
    n <- nrow(a$design)
    p <- ncol(a$design)
    if (a$rank != min(n, p)) stop("singular matrix 'a' in solve", call. = FALSE)
    if (missing(b)) {
        if (n != p) stop("only square matrices can be inverted", call. = FALSE)
        b <- diag(1, n)
    }
    ## as solve.qr(): coefficients outside the basis are zero, so that the
    ## result solves the system
    res <- qr.coef(a, b)
    res[is.na(res)] <- 0
    res
}
