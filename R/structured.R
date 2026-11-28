## Column-structured algebra on a plain design (see src/structured.c).
##
## Term t of the design has X_t[i, colmap_t[u, k]] = N_t(i)[u] C_t[F_t(i), k],
## with F_t(i) the combined level of its factors, N_t(i) the products of its
## numerics and C_t the Kronecker product of its contrasts. Products with the
## design are sums of N's within levels (computed in C), followed by small
## dense products with the C_t (here). Used only when every input is finite
## and no product can overflow; the results then agree with the dense
## computation up to rounding. Option `xtx.structured = FALSE`
## turns these kernels off (the row kernels are used instead).

## Is the structured path available for this plain design? Cached.
xtx_struct_ok <- function(d) {
    if (!isTRUE(getOption("xtx.structured", TRUE))) return(FALSE)
    if (xtx_is_comb(d) || !kernel_ok(d)) return(FALSE)
    cache <- d$cache
    if (!is.null(cache) && !is.null(cache$struct_ok)) return(cache$struct_ok)
    ok <- .Call(C_xtx_struct_ok, d$setup, d$vars, 2^23)
    if (!is.null(cache)) cache$struct_ok <- ok
    ok
}

## Term descriptions: off, colmap (0-based), fvar (the variable of each
## factor), fcm (their contrast matrices; C_t is their Kronecker product, the
## first factor varying fastest), cols (1-based global columns, ns x kc) and
## L (the combined levels). Cached.
xtx_struct_info <- function(d) {
    cache <- d$cache
    if (!is.null(cache) && !is.null(cache$struct_info)) return(cache$struct_info)
    info <- lapply(.Call(C_xtx_struct_info, d$setup, d$vars), function(t) {
        names(t) <- c("off", "colmap", "fvar", "fcm")
        t$cols <- t$off + t$colmap + 1L
        t$L <- prod(vapply(t$fcm, nrow, 1L))
        t
    })
    if (!is.null(cache)) cache$struct_info <- info
    info
}

## Products with Kronecker products of contrasts, one factor (mode) at a
## time, instead of with the dense L x kc matrix C_t: X holds an array whose
## first dims are the factors' levels; mats[[f]] is list(A) (contract the
## mode with A: sum_a X[a, ...] A[a, i]) or list(A, B) (a factor two terms
## share: sum_a X[a, ...] A[a, i] B[a, j]). Returns the result with the
## remaining dims first, then the new dims in order. The cost is the size
## of the moments times the non-zero entries per level, not L x kc.
xtx_kron_modes <- function(X, mats) {
    for (m in mats) X <- .Call(C_xtx_kron_mode, X, m[[1L]], if (length(m) > 1L) m[[2L]])
    X
}

## t(X) %*% diag(w) %*% X, w NULL or finite. The moments are by cell of
## the union of each pair of terms' factors and are contracted with the
## factors' contrasts one at a time, in C (gram_block() in src/structured.c).
xtx_struct_gram <- function(d, w = NULL) {
    .Call(C_xtx_struct_gram, d$setup, d$vars, if (is.null(w)) NULL else as.double(w), design_threads())
}

## t(X) %*% Y, Y an n x K finite matrix.
xtx_struct_xty <- function(d, Y) {
    info <- xtx_struct_info(d)
    mom <- .Call(C_xtx_struct_xty, d$setup, d$vars, Y, design_threads())
    K <- ncol(Y)
    out <- matrix(0, ncol(d), K)
    for (t in seq_along(info)) {
        Tt <- info[[t]]
        out[Tt$cols, ] <- xtx_struct_contract(mom[[t]], Tt, K)
    }
    out
}

## Moments of one term, [L, ns, K], contracted with its contrasts: the
## (ns kc) x K matrix whose rows follow the term's columns (Tt$cols).
xtx_struct_contract <- function(M, Tt, K) {
    nf <- length(Tt$fcm)
    A <- xtx_kron_modes(M, lapply(Tt$fcm, list))
    if (nf == 0L) return(A)
    dim(A) <- c(nrow(Tt$cols), K, vapply(Tt$fcm, ncol, 1L))
    A <- aperm(A, c(1L, 2L + seq_len(nf), 2L))
    dim(A) <- c(length(Tt$cols), K)
    A
}

## X %*% B, B a p x K finite matrix.
xtx_struct_mult <- function(d, B) {
    info <- xtx_struct_info(d)
    K <- ncol(B)
    P <- lapply(info, function(Tt) {
        ## P_t = C_t B_t as [L, ns, K]: B_t is [ns, k_1, .., k_m, K]
        ns <- nrow(Tt$cols)
        nf <- length(Tt$fcm)
        Bt <- B[Tt$cols, , drop = FALSE]
        if (nf == 0L) return(Bt)
        dim(Bt) <- c(ns, vapply(Tt$fcm, ncol, 1L), K)
        Bt <- aperm(Bt, c(1L + seq_len(nf), 1L, nf + 2L))
        A <- xtx_kron_modes(Bt, lapply(Tt$fcm, function(cm) list(t(cm))))
        dim(A) <- c(ns * K, Tt$L)
        t(A)
    })
    .Call(C_xtx_struct_mult, d$setup, d$vars, P, K, design_threads())
}

## diag(X V t(X)) = rowSums((X %*% V) * X), V a p x p finite matrix, or NULL
## when the tables below would hold more than `cap` numbers (the caller
## then uses the row kernel). The columns of the design split into its
## dense part (the terms with one level: the intercept and the products of
## numeric variables, which are non-zero in every row) and its terms with
## factors, which have a few non-zero entries per row. The dense part goes
## through a matrix product, block of rows by block of rows; the terms
## with factors read small tables: the symmetric part S of V contracted
## with their contrasts (see xtx_struct_quad() in src/structured.c). The
## split follows from the terms of the formula, not from a threshold on
## the density of the rows.
xtx_struct_quad <- function(d, V, cap = 2^23) {
    info <- xtx_struct_info(d)
    L <- vapply(info, function(Tt) as.double(Tt$L), 0)
    ns <- vapply(info, function(Tt) nrow(Tt$cols), 1L)
    kc <- vapply(info, function(Tt) ncol(Tt$cols), 1L)
    fac <- which(L > 1)
    dcols <- unlist(lapply(info[L == 1], function(Tt) Tt$cols[, 1L]))
    qd <- length(dcols)
    pairs <- if (length(fac) > 1L) utils::combn(fac, 2L) else matrix(0L, 2L, 0L)
    size <- sum(L[fac] * (kc[fac] + ns[fac] * (qd + ns[fac]))) +
        sum(L[pairs[1L, ]] * L[pairs[2L, ]] * ns[pairs[1L, ]] * ns[pairs[2L, ]])
    if (size > cap) return(NULL)
    S <- (V + t(V)) / 2
    storage.mode(S) <- "double"
    ## C_t, L x kc: the Kronecker product of the term's contrasts, the
    ## first factor varying fastest
    C <- vector("list", length(info))
    for (t in fac) C[[t]] <- Reduce(function(acc, cm) kronecker(cm, acc), info[[t]]$fcm)
    A <- Mself <- vector("list", length(info))
    for (t in fac) {
        cols <- info[[t]]$cols
        At <- array(0, c(qd, L[t], ns[t]))
        Mt <- array(0, c(L[t], ns[t], ns[t]))
        for (u in seq_len(ns[t])) {
            if (qd > 0L) At[, , u] <- S[dcols, cols[u, ], drop = FALSE] %*% t(C[[t]])
            for (v in seq_len(ns[t])) {
                Mt[, u, v] <- rowSums((C[[t]] %*% S[cols[u, ], cols[v, ], drop = FALSE]) * C[[t]])
            }
        }
        A[[t]] <- At
        Mself[[t]] <- Mt
    }
    Mpair <- lapply(seq_len(ncol(pairs)), function(j) {
        s <- pairs[1L, j]
        t <- pairs[2L, j]
        cs <- info[[s]]$cols
        ct <- info[[t]]$cols
        M <- array(0, c(L[s], L[t], ns[s], ns[t]))
        for (u in seq_len(ns[s])) for (v in seq_len(ns[t])) {
            M[, , u, v] <- C[[s]] %*% S[cs[u, ], ct[v, ], drop = FALSE] %*% t(C[[t]])
        }
        M
    })
    pairs0 <- pairs - 1L
    storage.mode(pairs0) <- "integer"
    .Call(C_xtx_struct_quad, d$setup, d$vars, S[dcols, dcols, drop = FALSE], A, Mself,
          pairs0, Mpair, design_threads())
}
