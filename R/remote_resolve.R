## The remote frame of a model: a formula and remote data (see remote.R)
## with weights, an offset and a subset, turned into the remote design of
## the rows used, the response and the arguments as columns over them (in
## SQL, or evaluated in R, batch by batch, into a temporary table):
## xtx_frame_remote(), the remote path of design_frame(), on which the
## input stage of the fits is built (xtx.fit). Also design_resolve()
## (arguments, and a response, as columns of an existing design).

## The arguments of remote data that design_matrix() and design_frame() take.
xtx_remote_keys <- c("table", "sql", "engine", "chunk_size")

## The remote path of design_frame() (its arguments, documented there).
## `weights`, `offset` and `subset` are taken unevaluated; the other
## arguments in `...` are bare names of columns that must be present on the
## rows used. `.caller`: the frame design_frame() was called from, where R
## expressions are evaluated when the formula has no environment.
## The steps, in order (the errors of each are those of lm() and glm() on a
## data frame where they have one): `na.action`; unknown arguments; the
## source (one query for its columns); the response of the formula, by the
## type of its columns (xtx_response_spec()); the types of the columns
## given as weights and offset (xtx_frame_arg_types()); the design and the
## columns (xtx_remote_frame()); the response over the rows used
## (xtx_response_result()).
xtx_frame_remote <- function(formula, data, ..., subset, weights, na.action, offset, contrasts = NULL,
                             table = NULL, sql = NULL, engine = NULL, chunk_size = NULL, .caller = parent.frame()) {
    na <- if (missing(na.action)) "omit" else xtx_frame_na_kind(na.action)
    extra <- eval(substitute(alist(...)))
    ## (other arguments: unnamed bare names of columns that must be present;
    ## a named one is a misspelt argument, such as `weight = w`)
    nm <- names(extra) %||% rep("", length(extra))
    bad <- !vapply(extra, is.name, NA) | nzchar(nm)
    if (any(bad)) {
        nm <- ifelse(nzchar(nm), nm, vapply(extra, deparse2, ""))
        stop("Unknown arguments: ", paste(nm[bad], collapse = ", "), call. = FALSE)
    }
    xtx_need_duckdb()
    src <- xtx_source(data, table, sql)
    ## (R expressions are evaluated in the formula's environment, as
    ## model.frame() evaluates them in memory)
    env <- environment(formula) %||% .caller
    tt <- stats::terms(formula, data = src$proto)
    resp <- attr(tt, "response")
    spec <- if (resp > 0L) xtx_response_spec(attr(tt, "variables")[[resp + 1L]], src$types, src$holder$con)
            else list(kind = "none", ncol = 0L, args = list(), present = character())
    args <- list(weights = if (!missing(weights)) substitute(weights),
                 offset = if (!missing(offset)) substitute(offset),
                 subset = if (!missing(subset)) substitute(subset))
    tv <- xtx_terms_vars(tt)
    xtx_frame_arg_types(src$types, args$weights, args$offset,
                        lapply(attr(tt, "offset"), function(k) tv$exprs[[k]][[2L]]))
    dots <- Filter(Negate(is.null), list(table = table, sql = sql, engine = engine, chunk_size = chunk_size))
    f <- xtx_remote_frame(formula, src, c(spec$args, args), env, dots, contrasts, logical = "subset",
                          present = unique(c(vapply(extra, as.character, ""), spec$present)), na = na)
    D <- f$x
    structure(list(x = D, response = xtx_response_result(D, spec, f$cols, f$yk),
                   weights = f$cols$weights, offset = f$off, ncol = c(weights = f$wk, offset = f$offk),
                   types = src$types, terms = tt, contrasts = D$plan$contrasts, xlevels = D$plan$xlevels,
                   na.action = na),
              class = "xtx_remote_frame")
}

#' @export
print.xtx_remote_frame <- function(x, ...) {
    d <- dim(x$x)
    cat(sprintf("<remote model frame: %s rows, design of %d columns>\n", xtx_fmt_count(d[1L]), d[2L]))
    r <- x$response
    cat("response:", switch(r$kind, none = "none", matrix = sprintf("matrix of %d columns", r$ncol),
                            factor = sprintf("factor with %d levels", length(r$levels)), r$kind), "\n")
    given <- c(weights = !is.null(x$weights), offset = !is.null(x$offset))
    if (any(given)) cat("with:", paste(names(given)[given], collapse = ", "), "\n")
    invisible(x)
}

## na.action of a remote frame as one of "omit", "exclude", "fail", "pass":
## the functions of stats, or their names.
xtx_frame_na_kind <- function(na.action) {
    if (is.null(na.action)) return("omit")
    for (nm in c("omit", "exclude", "fail", "pass")) {
        f <- get(paste0("na.", nm), envir = asNamespace("stats"))
        if (identical(na.action, f) || identical(na.action, paste0("na.", nm))) return(nm)
    }
    stop("`na.action` must be na.omit, na.exclude, na.fail, na.pass.", call. = FALSE)
}

## The left-hand side `lhs` of a formula on remote data, by the types of
## the columns it names, as model.frame() returns a response of any type:
## list(kind, ncol, args, present).
##   kind "numeric": an expression (a numeric column, or any expression
##     over the columns), resolved as argument `y`;
##   kind "logical": a logical column, as 0/1 (`y`);
##   kind "factor": a character or factor column, which no SQL makes
##     numeric: it only selects the rows where it is present (`present`);
##   kind "other": a column of another type (a date, say), likewise;
##   kind "matrix": cbind() of `ncol` expressions; the first two are
##     resolved as arguments `a` and `b`, and the columns named beyond
##     them only select rows. cbind() of one expression is that expression,
##     as a one-column matrix is its column in memory (but cbind() of a
##     factor gives its codes there: it stays a matrix of one column, the
##     call itself resolved as `a`, which R evaluates).
xtx_response_spec <- function(lhs, types, con) {
    col <- function(e) if (is.name(e) && as.character(e) %in% names(types)) as.character(e)
    ty <- function(e) if (is.null(col(e))) "" else types[[col(e)]]
    if (is.call(lhs) && identical(lhs[[1L]], as.name("cbind"))) {
        es <- unname(as.list(lhs)[-1L])
        k <- length(es)
        if (k == 1L && ty(es[[1L]]) != "factor") return(xtx_response_spec(es[[1L]], types, con))
        if (k == 1L) return(list(kind = "matrix", ncol = 1L, args = list(a = lhs), present = col(es[[1L]])))
        return(list(kind = "matrix", ncol = k, args = if (k >= 2L) list(a = es[[1L]], b = es[[2L]]) else list(),
                    present = as.character(unlist(lapply(if (k >= 2L) es[-(1:2)] else es, col)))))
    }
    ## (factor(x), as.factor(x), as.character(x) of a text, ENUM or integer
    ## column: the factor of the column, xtx_factor_call())
    fc <- xtx_factor_call(lhs, types)
    if (!is.null(fc)) return(list(kind = "factor", ncol = 1L, args = list(), present = fc$col, text = fc$text))
    t <- ty(lhs)
    if (t %in% c("character", "factor", "other"))
        return(list(kind = if (t == "other") "other" else "factor", ncol = 1L, args = list(), present = col(lhs)))
    if (t == "logical") {
        y <- DBI::SQL(sprintf("CAST(s.%s AS DOUBLE)", DBI::dbQuoteIdentifier(con, col(lhs))))
        return(list(kind = "logical", ncol = 1L, args = list(y = y), present = character()))
    }
    list(kind = "numeric", ncol = 1L, args = list(y = lhs), present = character())
}

## The response of a remote frame over the rows of design D, from its
## description (xtx_response_spec()) and the columns resolved for it
## (`cols`: y, or a and b): list(kind, ncol) with
##   numeric, logical: `y`, a column; `ncol` is more than 1 when the
##     expression, evaluated in R, gave a matrix (`yk` columns; `y` is then
##     its first column);
##   factor: `y`, the values as a text column (which compares with a
##     string: y == "a"), and `levels`, those of the rows of D (in the
##     order of the levels of a factor column; sorted for a character
##     column): one query;
##   matrix: `cols`, the columns resolved (the first two, when `ncol` is at
##     least 2; the codes of cbind() of one factor column);
##   other, none: nothing.
xtx_response_result <- function(D, spec, cols, yk = 1L) {
    out <- list(kind = spec$kind, ncol = spec$ncol)
    if (spec$kind %in% c("numeric", "logical")) {
        out$y <- cols$y
        out$ncol <- yk
    } else if (spec$kind == "matrix") {
        out$cols <- unname(Filter(Negate(is.null), cols[c("a", "b")]))
    } else if (spec$kind == "factor") {
        col <- spec$present
        q <- xtx_qi(D, col)
        vals <- xtx_query(D, sprintf("SELECT DISTINCT s.%s AS v FROM (%s) AS s WHERE %s", q, D$src$select, D$filter))$v
        proto <- D$src$proto[[col]]
        f <- if (is.factor(proto) && !isTRUE(spec$text)) factor(as.character(vals), levels = levels(proto))
             else factor(if (isTRUE(spec$text)) as.character(vals) else vals)
        out$levels <- levels(droplevels(f))
        out$y <- xtx_column(D, sprintf("CAST(s.%s AS VARCHAR)", q), type = "character")
    }
    out
}

## The types of the columns that a model on remote data takes as its
## `weights` and `offset` (and the expressions of its offset() terms,
## `offt`), checked before their SQL translation, with the errors of lm()
## and glm() on a data frame, in their order: the weights, then the offset
## (model.offset()). Only bare names of columns are checked here (an
## expression over a column of the wrong type fails when evaluated, as in
## memory).
xtx_frame_arg_types <- function(types, weights, offset, offt) {
    ty <- function(e) if (is.name(e) && as.character(e) %in% names(types)) types[[as.character(e)]] else NA_character_
    wt <- ty(weights)
    if (!is.na(wt) && !wt %in% c("double", "integer")) stop("'weights' must be a numeric vector", call. = FALSE)
    omsg <- xtx_offset_type_msg(ty(offset), vapply(offt, ty, ""))
    if (!is.null(omsg)) stop(omsg, call. = FALSE)
    invisible()
}

## The error of model.offset() for an `offset` argument and offset() terms
## of these column types (NA: not a column, taken as numeric), or NULL: the
## argument alone must be numeric; with offset() terms, model.offset()
## adds them up (`+` fails on text, gives NA on a factor) and the sum must
## be numeric.
xtx_offset_type_msg <- function(arg, terms) {
    num <- c("double", "integer")
    if (!length(terms)) return(if (!is.na(arg) && !arg %in% num) "'offset' must be numeric")
    s <- if (is.na(arg)) "double" else arg
    for (t in terms) {
        if (is.na(t)) t <- "double"
        if ("character" %in% c(s, t)) return("non-numeric argument to binary operator")
        s <- if ("factor" %in% c(s, t)) "logical" else if ("other" %in% c(s, t)) "other"
             else if (s == "logical" && t == "logical") "integer" else "double"
    }
    if (!s %in% num) "'offset' must be numeric"
}

## The remote design of a model and the columns of its arguments (one code
## path for linear models, GLMs and fixed effects).
##   formula: the model formula; the design is that of its right-hand side
##     (offset() terms included; their values are columns too), and its
##     response is given by the caller in `args` (`y`, or `a` and `b`).
##   src: an opened source (xtx_source()).
##   args: named list of unevaluated arguments (NULL: absent), resolved in
##     env: column names, DBI::SQL() or R expressions over the columns.
##     Numeric, except those named in `logical` (`subset`, which restricts
##     the rows).
##   present: columns that must not be missing (FE and cluster variables, a
##     factor response).
##   na: "omit" (or "exclude": the same rows), "fail" or "pass" (lm()'s
##     errors when a row has a missing value).
## R expressions with no SQL translation are evaluated in R into a
## temporary table (`engine = "sql"` refuses). The
## rows used are those of the subset where no model variable, argument or
## `present` column is missing, and factor levels are those of these rows,
## as model.frame(drop.unused.levels = TRUE) gives them in memory.
## Returns list(x = the design, cols = named list of columns (NULL for
## absent arguments), off = the offset (offset() terms plus the `offset`
## argument) or NULL, offk = the largest number of columns of an offset
## evaluated in R (a matrix: the caller reports its length, as in memory),
## wk = the number of columns of the weights evaluated in R (likewise),
## yk = that of the response `y` (likewise)).
xtx_remote_frame <- function(formula, src, args, env, dots = list(), contrasts = NULL, logical = character(),
                             present = character(), na = "omit") {
    miss <- setdiff(present, names(src$types))
    if (length(miss)) stop(sprintf("No column `%s` in the remote data.", miss[1L]), call. = FALSE)
    tt <- stats::terms(formula, data = src$proto)
    tv <- xtx_terms_vars(tt)
    rhs <- if (attr(tt, "response") > 0L) stats::delete.response(tt) else tt
    attr(rhs, "predvars") <- NULL
    rhs <- stats::formula(rhs)
    environment(rhs) <- environment(formula)
    ## offset() terms are arguments too: their SQL, or the column computed
    ## in R, replaces the variable of the design (`override`)
    offt <- list()
    for (k in attr(tt, "offset")) offt[[tv$names[[k]]]] <- tv$exprs[[k]][[2L]]
    all <- c(Filter(Negate(is.null), args), offt)
    label <- function(nm) if (nm %in% names(offt)) "offset term" else if (nm %in% c("y", "a", "b")) "response" else nm
    con <- src$holder$con
    q <- function(x) as.character(DBI::dbQuoteIdentifier(con, x))
    res <- list()
    items <- list()
    exprs <- list()
    mat <- list()
    for (nm in names(all)) {
        a <- xtx_resolve_or_eval(all[[nm]], env, src$types, con, label(nm), bool = nm %in% logical)
        if (is.null(a$expr)) {
            res[[nm]] <- a
            next
        }
        col <- paste0("__xtx_eval", length(exprs) + 1L)
        items[[label(nm)]] <- a
        exprs[[col]] <- a$expr
        ## (`a` alone: cbind() of one factor column, a one-column matrix)
        if (nm %in% c("y", "offset", "weights") || nm %in% names(offt) || (nm == "a" && is.null(all[["b"]])))
            mat[[col]] <- "count"
        res[[nm]] <- list(sql = if (nm %in% logical) sprintf("s.%s", q(col)) else sprintf("CAST(s.%s AS DOUBLE)", q(col)),
                          inputs = col, logical = nm %in% logical)
    }
    if (length(exprs)) {
        dots <- xtx_untranslated_engine(items, dots)
        lg <- unlist(lapply(res, function(r) if (isTRUE(r$logical)) r$inputs))
        src <- xtx_materialize(src, exprs, env, logical = lg,
                              chunk_size = dots$chunk_size %||% getOption("xtx.chunk_size", 250000),
                              matrix = mat)
    }
    ## (the number of columns of the offsets and weights evaluated in R)
    nc <- attr(src, "ncol")
    ncols <- function(nms) max(c(1L, unlist(nc[unlist(lapply(res[intersect(nms, names(res))], `[[`, "inputs"))])))
    offk <- ncols(c("offset", names(offt)))
    wk <- ncols("weights")
    yk <- ncols("y")
    dots$table <- dots$sql <- NULL
    nonlog <- setdiff(names(res), c(logical, names(offt)))
    extra <- stats::setNames(lapply(res[nonlog], function(r) list(sql = r$sql, inputs = r$inputs, kind = "numeric")),
                             sprintf("(%s)", nonlog))
    for (col in present) extra[[paste0("(present:", col, ")")]] <- list(inputs = col, kind = "present")
    override <- lapply(res[names(offt)], function(r) list(sql = r$sql, inputs = r$inputs))
    if (offk > 1L && length(offt)) {
        ## (an offset() term gave a matrix, which the caller reports with
        ## the number of rows used: the offset() terms are then columns
        ## that only select rows, and the design that of the other terms)
        for (nm in names(offt)) extra[[sprintf("(%s)", nm)]] <- c(override[[nm]], kind = "numeric")
        override <- list()
        tr <- stats::terms(rhs, data = src$proto)
        tl <- attr(tr, "term.labels")
        rhs <- stats::reformulate(if (length(tl)) tl else "1", intercept = attr(tr, "intercept") > 0L,
                                  env = environment(rhs))
    }
    sub <- if (any(logical %in% names(res))) paste(sprintf("(%s)", vapply(res[intersect(logical, names(res))], `[[`, "", "sql")), collapse = " AND ")
    ## (no rows left: the error of lm())
    D <- tryCatch(xtx_remote(rhs, src, contrasts, NULL, dots, drop_unused = TRUE, subset = sub, extra = extra,
                            override = override),
                  xtx_no_rows = function(e) stop("0 (non-NA) cases", call. = FALSE))
    if (na %in% c("fail", "pass")) xtx_remote_fit_na(D, res, na, nonlog, names(offt))
    cols <- lapply(names(args), function(nm) if (!is.null(res[[nm]])) xtx_column(D, res[[nm]]$sql))
    names(cols) <- names(args)
    off <- c(res[["offset"]]$sql, vapply(res[names(offt)], `[[`, "", "sql"))
    list(x = D, cols = cols, off = if (length(off)) xtx_column(D, paste(sprintf("(%s)", off), collapse = " + ")),
         offk = offk, wk = wk, yk = yk)
}

## na.fail and na.pass: the errors of lm() when rows have missing values.
xtx_remote_fit_na <- function(D, res, na, args, offt) {
    total <- as.numeric(xtx_query(D, sprintf("SELECT COUNT(*) AS n FROM (%s) AS s", D$src$select))$n)
    if (total == D$n) return(invisible())
    if (na == "fail") stop("missing values in object", call. = FALSE)
    bad <- function(nm) {
        r <- res[[nm]]
        if (is.null(r)) return(FALSE)
        as.numeric(xtx_query(D, sprintf("SELECT COUNT(*) AS n FROM (%s) AS s WHERE (%s) IS NULL OR isnan(%s)",
                                       D$src$select, r$sql, r$sql))$n) > 0
    }
    if (any(vapply(c("y", "a", "b", "offset", offt), bad, NA))) stop("NA/NaN/Inf in 'y'", call. = FALSE)
    if (bad("weights")) stop("missing or negative weights not allowed", call. = FALSE)
    stop("NA/NaN/Inf in 'x'", call. = FALSE)
}

## R expressions over the remote data with no SQL translation (a response,
## an offset() term, the `weights`, `offset` or `subset` of a fit): they are
## evaluated in R, batch by batch, and written with all the columns of the
## source to a temporary table, whose new columns the fit then reads like
## any other. One pass over the data; the rows, levels, missing values and
## compressed sums are then exact, as for translated expressions.
##   exprs: named list of R expressions (the names of the new columns);
##   logical: names of the expressions to keep as logical (subset);
##   model: for the expressions that are model variables (named as in
##     exprs), list(inputs, factor, label), see xtx_eval_model(). Their
##     columns and kinds are returned in the attribute "model" of the
##     result: list(kind, cols) per expression.
##   matrix: for expressions that may give a matrix with one row per row
##     (named as in exprs): "count" (an offset or weights: its first column
##     is written, and its number of columns returned in the attribute
##     "ncol" of the result, for the caller's error), or the message of the
##     error for a matrix with several columns. A one-column matrix is its
##     column.
## Returns the new source (xtx_source() of the temporary table),
## which keeps the original one alive.
xtx_materialize <- function(src, exprs, env, logical = character(), chunk_size = 250000, model = list(),
                           matrix = list()) {
    info <- list()
    ncols <- list()
    con <- src$holder$con
    q <- function(x) as.character(DBI::dbQuoteIdentifier(con, x))
    from <- src$base_table
    raw <- NULL
    if (is.null(from)) {
        raw <- q(xtx_tmp_name("xtx_raw"))
        DBI::dbExecute(con, sprintf("CREATE TEMPORARY TABLE %s AS %s", raw, src$select))
        on.exit(DBI::dbExecute(con, sprintf("DROP TABLE IF EXISTS %s", raw)), add = TRUE)
        from <- raw
    }
    r <- DBI::dbGetQuery(con, sprintf("SELECT MIN(rowid) AS lo, MAX(rowid) AS hi FROM %s", from))
    out <- xtx_tmp_name("xtx_eval")
    created <- FALSE
    write <- function(df) {
        if (!created) {
            DBI::dbWriteTable(con, out, df, temporary = TRUE)
            created <<- TRUE
        } else {
            DBI::dbAppendTable(con, out, df)
        }
    }
    evaluate <- function(chunk) {
        for (nm in names(exprs)) {
            if (!is.null(model[[nm]])) {
                r <- xtx_eval_model(exprs[[nm]], chunk, env, model[[nm]], nm, info[[nm]])
                info[[nm]] <<- r$info
                for (cn in names(r$cols)) chunk[[cn]] <- r$cols[[cn]]
                next
            }
            v <- eval(exprs[[nm]], chunk, env)
            if (length(v) == 1L && nrow(chunk) != 1L) v <- rep(v, nrow(chunk))
            if (!is.null(matrix[[nm]]) && length(dim(v)) == 2L && nrow(v) == nrow(chunk)) {
                if (ncol(v) > 1L) {
                    if (!identical(matrix[[nm]], "count")) stop(matrix[[nm]], call. = FALSE)
                    ncols[[nm]] <<- max(ncols[[nm]] %||% 1L, ncol(v))
                }
                v <- v[, 1L]
            }
            if (length(v) != nrow(chunk) || !is.null(dim(v)))
                stop(sprintf("`%s` gives %d values for %d rows.", deparse2(exprs[[nm]]), length(v), nrow(chunk)),
                     call. = FALSE)
            if (nm %in% logical) {
                if (!is.logical(v)) stop(sprintf("`%s` is not logical.", deparse2(exprs[[nm]])), call. = FALSE)
            } else {
                if (!is.numeric(v) && !is.logical(v))
                    stop(sprintf("`%s` is not numeric.", deparse2(exprs[[nm]])), call. = FALSE)
                v <- as.double(v)
            }
            chunk[[nm]] <- v
        }
        chunk
    }
    lo <- as.numeric(r$lo)
    if (!is.na(lo)) {
        hi <- as.numeric(r$hi)
        a <- lo
        while (a <= hi) {
            chunk <- DBI::dbGetQuery(con, sprintf("SELECT * FROM %s WHERE rowid >= %.0f AND rowid < %.0f",
                                                  from, a, a + chunk_size))
            ## (evaluated before the write, so that its errors are not
            ## reported as errors of dbWriteTable())
            if (nrow(chunk)) {
                chunk <- evaluate(chunk)
                write(chunk)
            }
            a <- a + chunk_size
        }
    }
    if (!created) {
        chunk <- evaluate(src$proto)
        write(chunk)
    }
    new <- xtx_source(con, out, NULL)
    new$parent <- src$holder
    new$label <- paste(src$label, "(with columns computed in R)")
    attr(new, "model") <- info
    attr(new, "ncol") <- ncols
    new
}

## The value of a model variable (the R expression `expr`) on a batch of
## rows, as columns named after `nm`: one numeric column, or one per column
## of a numeric matrix (nm_1, nm_2, ...: poly(), ns()); a text column for a
## factor or character value, or when the specification holds levels for
## the variable (m$factor); a logical column. The expression sees only the
## rows where its input columns (m$inputs) are present; the others get
## missing values. m$class: the class of the variable in training (one
## named element of the plan's data classes), checked on each batch.
## `prev`: the kind and columns found on earlier batches,
## which every batch must give. Returns list(cols, info = list(kind, cols)).
xtx_eval_model <- function(expr, chunk, env, m, nm, prev = NULL) {
    n <- nrow(chunk)
    ok <- if (length(m$inputs)) stats::complete.cases(chunk[m$inputs]) else rep(TRUE, n)
    bad <- function(msg) stop(sprintf("`%s` %s.", m$label, msg), call. = FALSE)
    if (!any(ok) && !is.null(prev)) {
        v <- NULL
        info <- prev
    } else {
        v <- eval(expr, chunk[ok, , drop = FALSE], env)
        ## a variable of another type than in training: the error of
        ## predict.lm() (a character variable with levels is a factor, as
        ## model.frame() makes it)
        if (!is.null(m$class))
            stats::.checkMFClasses(m$class, stats::setNames(list(if (is.character(v) && isTRUE(m$factor)) factor(v) else v),
                                                            names(m$class)))
        k <- if (is.matrix(v)) ncol(v) else 1L
        if (!is.matrix(v) && length(v) == 1L && sum(ok) != 1L) v <- rep(v, sum(ok))
        if (NROW(v) != sum(ok)) bad(sprintf("gives %d values for %d rows", NROW(v), sum(ok)))
        kind <- if (is.matrix(v)) {
            if (!is.numeric(v) && !is.logical(v)) bad("gives a matrix that is not numeric")
            "numeric"
        } else if (isTRUE(m$factor) || is.factor(v) || is.character(v)) {
            "factor"
        } else if (is.logical(v)) {
            "logical"
        } else if (is.numeric(v)) {
            "numeric"
        } else {
            bad(sprintf("gives values of class %s, which a remote design cannot use", class(v)[1L]))
        }
        info <- list(kind = kind, cols = if (k == 1L) nm else paste0(nm, "_", seq_len(k)))
        if (!is.null(prev) && !identical(prev, info)) bad("gives values of different kinds on different rows")
    }
    fill <- function(x, na) {
        out <- rep(na, n)
        if (!is.null(v)) out[ok] <- x
        out
    }
    cols <- lapply(seq_along(info$cols), function(u) {
        switch(info$kind,
               numeric = fill(if (is.matrix(v)) as.double(v[, u]) else as.double(v), NA_real_),
               factor = fill(as.character(v), NA_character_),
               logical = fill(v, NA))
    })
    list(cols = stats::setNames(cols, info$cols), info = info)
}

## Resolve an argument (xtx_resolve_arg()), or, when its R expression has no
## SQL translation, return it for evaluation in R: list(sql, inputs) or
## list(expr, reason).
xtx_resolve_or_eval <- function(x, env, types, con, what, bool = FALSE) {
    if (is.null(x)) return(NULL)
    if (is.call(x) && !deparse2(x[[1L]]) %in% c("SQL", "DBI::SQL")) {
        tr <- if (bool) xtx_sql_bool(x, types, con) else xtx_sql_num(x, types, con)
        if (is.null(tr$sql)) return(list(expr = x, reason = tr$reason))
    }
    xtx_resolve_arg(x, env, types, con, what, bool)
}

## The expressions of a fit that SQL cannot compute, as one message:
## "weights `myfun(x)` cannot be translated to SQL (myfun() has no SQL
## translation)".
xtx_untranslated_msg <- function(items) {
    nm <- names(items)[1L]
    sprintf("%s `%s` cannot be translated to SQL (%s)", nm, deparse2(items[[1L]]$expr), items[[1L]]$reason)
}

## A fit with expressions evaluated in R: a message names the first one;
## `engine = "sql"` refuses.
xtx_untranslated_engine <- function(items, dots) {
    if (!length(items)) return(dots)
    why <- xtx_untranslated_msg(items)
    if (identical(dots$engine, "sql")) stop("The SQL engine cannot compute this fit: ", why, ".", call. = FALSE)
    message("Evaluated in R: ", why, ".")
    dots
}

#' Arguments of a model as columns of a design
#'
#' `design_resolve()` turns arguments of a model, given unevaluated, into
#' columns over the rows of a design: for a remote design, column names,
#' `DBI::SQL()` expressions, or R expressions over the columns (those with a
#' SQL translation run in the database; the others are evaluated in R,
#' batch by batch, into a temporary table, which becomes the source of the
#' design returned); in memory, the arguments are evaluated in `data`.
#' The remote path of [design_frame()] does the same for a formula and
#' remote data, and builds the design too.
#'
#' Names may change before the first release.
#'
#' @param x A design: a `design_matrix` or a remote design.
#' @param exprs A named list of unevaluated arguments (a column name, a
#'   `DBI::SQL()` expression, an R expression), `NULL` for an absent one.
#' @param env The environment in which R expressions are evaluated.
#' @param logical The names of the arguments that are logical (such as a
#'   subset); the others are numeric.
#' @param data In memory, the data frame of the rows of `x` in which the
#'   arguments are evaluated.
#' @param response For a remote design, `NULL` or the unevaluated left-hand
#'   side of a formula (a column name, `cbind()` of two, an expression):
#'   it is resolved by the types of its columns and returned as `response`,
#'   as the remote path of [design_frame()] returns it. No rows are dropped:
#'   the rows are those of `x`.
#' @return A list with `x` (the design, on a new source when an argument
#'   was evaluated in R), `cols` (a named list of columns: vectors in
#'   memory, [design_column()]s remotely; `NULL` for absent arguments),
#'   `evaluated` (the arguments evaluated in R, with the reason) and, with
#'   `response`, `response` (see [design_frame()]).
#' @examples
#' d <- data.frame(y = rnorm(5), x = 1:5)
#' D <- design_matrix(y ~ x, d)
#' design_resolve(D, list(w = quote(x / 2)), environment(), data = d)$cols$w
#' @keywords internal
#' @export
design_resolve <- function(x, exprs, env = parent.frame(), logical = character(), data = NULL, response = NULL) {
    checkmate::assert_multi_class(x, c("xtx_matrix", "xtx_remote"))
    checkmate::assert_list(exprs, names = "unique")
    if (is.null(response)) return(xtx_p_resolve(x, exprs, env, logical = logical, data = data))
    if (!inherits(x, "xtx_remote")) stop("`response` applies only to a remote design.", call. = FALSE)
    spec <- xtx_response_spec(response, x$src$types, xtx_con(x))
    if (any(names(spec$args) %in% names(exprs))) stop("With `response`, `exprs` cannot be named `y`, `a` or `b`.", call. = FALSE)
    r <- xtx_p_resolve(x, c(spec$args, exprs), env, logical = logical)
    r$response <- xtx_response_result(r$x, spec, r$cols)
    r$cols <- r$cols[names(exprs)]
    r
}
