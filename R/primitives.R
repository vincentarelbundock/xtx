## The primitive layer of the fitting core (Kata design.matrix#vavc).
##
## Fitting code (linear models, and on this layer GLMs and fixed effects) is
## written once against the operations below, which take a *design* and
## *columns*, and dispatch to two backends:
##
##   design                     column
##   ------------------------   -----------------------------------------
##   in memory: design_matrix   numeric vector (or matrix: several columns)
##   remote:    xtx_remote       design_column: SQL expressions over the rows
##                              of a frame (see below; made by xtx_column())
##
## The remote column class, design_column, is also the public column type
## (R/column.R: design_column(), as.vector(), as.matrix(), col_compute(),
## as.data.frame(), print/format, and the methods below for arithmetic). Everything else
## here is internal; names start with xtx_p_ (operations on designs), xtx_c_
## (operations on columns) or xtx_.
##
## ---------------------------------------------------------------------------
## Columns
## ---------------------------------------------------------------------------
##
## A column holds one value per row of a design (in the order of its rows in
## memory; remotely, rows have no order). In memory it is a double vector of
## length nrow(x), or an n x k matrix for k columns. Remotely it is an
## object of class design_column, made by
##
##   xtx_column(frame, sql)
##
## with `frame` a xtx_remote (its source `s` holds what the expressions read,
## its filter selects its rows) and `sql` a character vector, one SQL
## expression over `s` per column (k = length(sql)). Columns are made by
##
##   xtx_p_column(x, sql)    a column of design x from SQL over its source
##                         (in memory: `sql` is returned, i.e. pass a vector)
##   xtx_p_mult(x, B)        products with the design (below)
##   arithmetic on columns
##
## Operations on columns (same R code for both backends):
##   +  -  *  /  ^         (numbers recycle; columns must have the same k or 1)
##   ==  !=  <  >  <=  >=  &  |  !
##                         (NA where either side is NA or NaN, as in R)
##   is.na is.nan is.finite is.infinite (R's semantics; never NA)
##   exp log sqrt abs log1p expm1 log2 log10 cos sin tan floor ceiling sign
##                         (Math group; log(x, base) with a number or a
##                         column as base, as R: xtx_column_log())
##   xtx_c_pmin xtx_c_pmax xtx_c_ifelse xtx_c_rowsums xtx_c_cbind xtx_c_col(x, j) xtx_ncols
##   xtx_c_const(like, value) xtx_c_lgamma xtx_c_round (half to even)
##   xtx_c_is_finite (NaN-safe) xtx_c_present (not NA or NaN: !is.na)
##   xtx_c_on(col, x) (the column carried over to design x, made from its own)
##   xtx_c_count(cond, negate) (1 where cond is TRUE, to sum)
##   xtx_c_freeze(col, table) (computed once: remotely a projection node;
##                         with table = TRUE, a temporary table of the rows
##                         used with the columns of the source and the values)
##   xtx_c_pnorm xtx_c_qnorm xtx_c_dnorm xtx_c_pcauchy
##                         (remotely the DuckDB macros of xtx_sql_macros(), as
##                         R computes them: pnorm/dnorm/pcauchy to rounding,
##                         qnorm to one ulp; xtx_c_mexpm1m(u) = -expm1(-u))
##
## NA semantics: SQL NULL is R's NA, IEEE NaN is NaN; arithmetic propagates
## both as R does (NULL op x is NULL). Sums over rows (xtx_colsum) return NA
## when a value is NULL, as sum() does, and NaN when a value is NaN.
## Comparisons with NULL are NULL (R: NA).
##
## Products (xtx_p_mult, xtx_p_quad), lookups (xtx_c_lookup, xtx_c_dense_codes)
## and projections (xtx_c_freeze, col_compute(table = FALSE)) are *nodes*
## of the frame of their column: immutable values (the SQL of a product
## over the rows of its design; the join of a lookup table, already
## written; the expressions of a projection) that every column made from
## them carries. Nothing is created in the database for a node. An
## operation that reads rows (sums, crossprod, groups, xtx_p_collect,
## as.vector, as.matrix, as.data.frame, output) or stores them (xtx_c_freeze
## with a table, col_compute, materialize) takes all its columns at once
## and renders, in its one statement, the source relation with one derived
## projection or join per node that the columns name, in dependency order
## (see "frames, source relations and nodes" below). A read changes no
## column, so columns of one design combine in the same way whatever was
## read before: `D %*% b1 + D %*% b2` reads both products over the rows of
## the design. Columns combine when they are columns of the same design
## and their source relations are the data of the design or tables stored
## from it one on the other (the later one holds the columns of the
## earlier ones); two tables stored separately cannot be combined (no key
## aligns their rows). In memory a product is computed at once.
##
## Getting values into R: xtx_p_collect(col) returns the values (in memory the
## vector itself; remotely a SELECT of all rows: only for small data). The
## public as.vector() and as.matrix() call it, and col_compute() calls
## xtx_c_freeze(table = TRUE) (with a name, a permanent table of the same
## rows and columns).
##
## ---------------------------------------------------------------------------
## Operations on designs
## ---------------------------------------------------------------------------
##
##   xtx_p_scale(x, s)          the design with rows scaled by column s (s * x);
##                            scales multiply. Same rows and columns as x.
##   xtx_p_crossprod(x)         X'X, p x p matrix
##   xtx_p_crossprod(x, y)      y a design (x itself, possibly scaled: X'SY);
##                            or a column (k columns): X'c, p x k matrix
##   xtx_p_colsums(x, group)    column sums of x (possibly scaled): vector of
##                            length p; with group (a column, values used as
##                            group labels; NA a group of its own, with
##                            rowsum()'s warning): G x p matrix, groups sorted
##                            as rowsum() sorts them, row names the labels
##   xtx_p_groupsum(col, group) sums of a column (k columns) within groups:
##                            G x k matrix, as rowsum(col, group) (NA for a
##                            group with a missing value)
##   xtx_p_y(x, y)              the response of the qr() methods as a column
##   xtx_p_min_nonzero(x)       a lower bound on the non-zero entries (0 if none)
##   xtx_c_codes(codes, G)      a grouping as integer codes 1..G: grouped sums
##                            then return all G groups in code order, without
##                            labels, from one-pass C kernels in memory
##   xtx_p_colsums_groups(x, groups, extra)
##                            several coded groupings at once: per grouping,
##                            list(S = G x p, n = G sums of the scale, extra =
##                            G x k sums of scale * extra); one pass of the
##                            structured kernel in memory
##   xtx_p_fe_sums(x, w, groups, extra, pairs)
##                            xtx_p_colsums_groups() of w * x and the
##                            crosstabs of pairs of coded groupings at once
##                            (fewer passes in memory)
##   xtx_c_groupsums(col, groups, total)
##                            sums of col within each coded grouping (a G x k
##                            matrix each, as xtx_p_groupsum()) and of `total`
##                            over all rows (as xtx_colsum()), in one query
##                            remotely
##   xtx_p_crosstab(a, b, w)    for two grouping columns, the sums of w (counts
##                            without w) for each pair of groups: a
##                            Matrix::dgCMatrix, rows the sorted groups of a,
##                            columns those of b (labels as.character() of
##                            the values, as in memory); with xtx_c_codes(), the
##                            Ga x Gb matrix in code order
##   xtx_p_mult(x, B)           X B (B a p-vector or p x k matrix) as a column
##                            (in memory a vector for a vector B; remotely
##                            a node, see above)
##   xtx_p_quad(x, V)           rowSums((X V) * X), V p x p symmetric positive
##                            semi-definite, as a column (remotely nodes)
##   xtx_design_with(x, cols)   the design with the nodes of the columns of a
##                            list in its source, so that their SQL reads
##                            it; in memory x
##   xtx_colsum(col)           sums over rows: a number per column
##   xtx_c_sums(...)            sums over rows of several columns (one query)
##   xtx_p_colnorms(col)        Euclidean norms of the columns
##   xtx_nrows(x)              number of rows of a design or column
##   xtx_engine(x)             "sql" for a remote design, ""
##                            in memory (for the description of a fit)
##   xtx_p_patterns(x, z, w, compress, keys = NULL, fit = "lm")
##                            compressed least squares: when rows with the
##                            same design inputs are worth grouping, a list
##                            with x (the design of the patterns, in memory),
##                            z (the weighted mean of the column z per
##                            pattern), w (the sum of the weights per pattern,
##                            NULL when all 1), within (the weighted sum of
##                            squares of z around the pattern means) and rows
##                            (the number of rows); otherwise NULL. compress:
##                            "auto", TRUE or FALSE. A weighted least squares
##                            fit of z on x with weights w gives the fit on
##                            every row; its residual sum of squares plus
##                            `within` is that of the fit on every row. With
##                            `keys` (a column, k columns), patterns are also
##                            split by their values, returned as `keys`
##                            (one row per pattern). fit: "lm" or "glm",
##                            for the remote rule of "auto"
##                            (xtx_compress_rule(); a weighted "lm" is
##                            "wlm").
##   xtx_p_materialize(x, cols = NULL)
##                            the design over a copy of its rows: remotely a
##                            temporary table of the rows used with every
##                            column of its source (later queries scan it
##                            instead of recomputing joins), or, given
##                            `cols` (a list of the columns to carry over
##                            to it), the columns of the data and the
##                            internal ones (stored values, nodes) that the
##                            design and `cols` read; a new design (carry
##                            columns over with xtx_c_on()); in memory x
##   xtx_p_subset(x, keep)      the design of the rows where the logical column
##                            keep is TRUE (remotely also a DBI::SQL()
##                            predicate); a new remote design
##   xtx_p_resolve(x, exprs, env, logical, data)
##                            columns for the arguments of a fit (response,
##                            weights, offset, subset): column names, SQL(),
##                            R expressions (in SQL when translatable, else
##                            evaluated in R into a temporary table); returns
##                            list(x, cols, evaluated)
##   xtx_p_output(x, cols, memory, na, design, attrs)
##                            row-wise results of an algorithm: in memory
##                            memory(d), d the data frame of the named
##                            columns `cols` (one row per row of x; default
##                            d); remotely a temporary view with the columns
##                            of the remote data and `cols`, a string of
##                            class xtx_remote_table (padded by napredict(na),
##                            with attributes attrs), or with design = TRUE
##                            the remote design with `cols` added to its
##                            data. One code path computes the columns; only
##                            the container differs.
##
## Remote methods: crossprod, colSums, products and patterns run in the
## database; xtx_p_groupsum and
## xtx_p_crosstab are single GROUP BY queries.
##
## ---------------------------------------------------------------------------
## User-facing operations of a remote design (operator goal: every matrix
## operation of an in-memory design works on a remote one, through this
## layer; status as of this contract)
## ---------------------------------------------------------------------------
##
##   operation                    remote meaning                    status
##   dim, dimnames                rows, columns                     done
##   as.matrix(D, force)          collected; refused above          done
##                                remote_max_rows unless force
##   D %*% b, D %*% B             column (design_column, k = ncol(B)) done
##   crossprod(D), (D, s * D),    p x p / p x k matrices            done
##   (D, "col" / SQL() / column)
##   crossprod(D, numeric w)      no row alignment remotely: error; done (error)
##                                use a column or s * D
##   colSums/colMeans(s * D,      vectors / G x p matrices          done
##     group = col)
##   rowSums/rowMeans(D)          column (labels .rowsum/.rowmean)  done
##   quadratic_diagonal(D, V)     column (label .quad)              done
##   design_column(D, "y"/SQL())  column of the remote data         done
##   as.vector, as.matrix,        values in R / data frame (with    done
##     as.data.frame,             `keep` columns) / temporary or
##     col_compute (of a column)  permanent table
##   qr(D)                        xtx_qr_chol on X'X with lm()'s      done
##                                rule, suspects from products
##   qr.coef, qr.qty, qr.R, solve numbers (as in memory)            done
##   qr.fitted, qr.resid, qr.qy,  tables (.fitted, .resid, .qy,     done
##     qr.Q                       .Q1...)
##   as(D, "CsparseMatrix")       collected, same guard as          done
##                                as.matrix (option to raise)
##   D[i, ]                       i a logical remote column or a    done
##                                DBI::SQL() predicate; numeric or
##                                logical vectors refused
##   s * D, D / s, -D             row scaling (number, column       done
##                                name, SQL(), column)
##   D1 + D2, D1 - D2             refused (operator decision)       done (error)
##   finite_rows(D) (internal)    a logical column                  done
##
## Tests: test_primitives.R (the layer), test_parity_ops.R (every operation
## above on the same data, in memory and remote: a table, Parquet).

## ---- remote columns ----

## A remote column: `frame` (the row domain, source relation and nodes, see
## below) and, per output column, `sql` (its expression over the relation
## `s` of the frame), `notnull` (known never to be NULL: its sums need no
## count of missing values), `lgl` (logical-valued: BOOLEAN in SQL, cast to
## DOUBLE where a number is needed, xtx_column_arith()) and `type` (the
## type of a value read as it is stored, "character", "integer", "logical",
## "factor" or "other"; "" for a double or a computed expression, which can
## be NaN: xtx_c_present()). Optional: `label` (set by xtx_column_label() in
## the public calls, names the values for as.data.frame()) and `neg` (the
## logical column that `!` negated). The flags are per output column and
## survive slicing, cbind(), storage and translation to another design.
xtx_column <- function(frame, sql, notnull = FALSE, lgl = FALSE, type = "") {
    ## (parenthesized, so that expressions compose with any operator)
    sql <- as.character(sql)
    sql <- sprintf("(%s)", sql)
    k <- length(sql)
    structure(list(frame = frame, sql = sql, notnull = rep_len(notnull, k), lgl = rep_len(lgl, k),
                   type = rep_len(type, k)), class = "xtx_column")
}

## The per-column flags of a column.
xtx_column_lgl <- function(x) {
    u <- unclass(x)
    rep_len(u$lgl %||% FALSE, length(u$sql))
}
xtx_column_type <- function(x) {
    u <- unclass(x)
    rep_len(u$type %||% "", length(u$sql))
}

## The output columns `idx` of a column, with their flags and labels (the
## expressions as they are: not parenthesized again).
xtx_column_pick <- function(x, idx) {
    u <- unclass(x)
    out <- list(frame = u$frame, sql = u$sql[idx], notnull = (u$notnull %||% rep(FALSE, length(u$sql)))[idx],
                lgl = xtx_column_lgl(x)[idx], type = xtx_column_type(x)[idx])
    if (!is.null(u$label)) out$label <- u$label[idx]
    structure(out, class = "xtx_column")
}

xtx_is_column <- function(x) inherits(x, "xtx_column")

## The number of columns k of a column: of a remote column its number of
## SQL expressions, in memory NCOL(). (dim() of a remote column is c(n, k)
## and length() is n * k, as for a matrix in memory: the fitting core counts
## columns with xtx_ncols(), never with length().)
xtx_ncols <- function(x) if (xtx_is_column(x)) length(unclass(x)$sql) else NCOL(x)

## The number of rows of a remote column: the row count cached with its
## design (every frame keeps it), else one COUNT query, cached thereafter.
xtx_column_nrows <- function(x) {
    fr <- unclass(x)$frame
    n <- xtx_nrows(fr)
    if (!is.null(n)) return(n)
    xtx_column_valid(x)
    ## (the rows of the design: nodes add columns, not rows)
    fr <- xtx_frame_ready(fr)
    n <- as.numeric(xtx_query(fr, sprintf("SELECT COUNT(*) AS n FROM (%s) AS s WHERE %s",
                                         fr$src$select, fr$filter))$n)
    if (is.environment(fr$cache)) fr$cache$n <- n
    n
}

#' @export
dim.xtx_column <- function(x) c(xtx_column_nrows(x), xtx_ncols(x))

#' @export
length.xtx_column <- function(x) prod(dim.xtx_column(x))

## ---- frames, source relations and nodes ----
##
## The frame of a column is a remote design (its row domain `id`, filter,
## connection and schema) with
##
##   a source relation (the *root*, xtx_frame_root()): the data of the
##       design, or a relation stored from it (xtx_root_child()). A root is
##       the `src` record of a design (holder, select, proto, types) with
##       `path`, the tokens of the stores it descends from (none for the
##       data), and `carried`, the ids of the nodes whose columns it holds;
##   nodes (xtx_frame_nodes()): the dependencies that the expressions of
##       the column name and the root does not hold, in dependency order.
##
## A node (xtx_node()) is an immutable value: its id (the stem of the names
## of its columns, unique in the session), the SQL of its columns over the
## relation `s` below it (`items`), a join (`join`) and the ids of the nodes
## it reads (`deps`). The kinds are "mult" (a product X B: its SQL is that
## of xtx_sql_predict(), written when the product is made), "lookup" (a
## table of values joined on a key) and "proj" (expressions computed once:
## col_compute(table = FALSE)). Nothing is created in the database for a
## node, except the table of values of a lookup, written when it is made.
##
## A read renders the relation of its frame (xtx_frame_ready()): the root,
## and on top of it one derived projection (or join) per node that the SQL
## of the read names, each once, in dependency order. Rendering changes
## neither the frame nor the columns: the same columns give the same SQL at
## their next read, whatever was read in between.
##
## A design keeps the same fields, with its nodes (those its scale and
## filter name) already rendered in `src` (and the root in `base`), so that
## code that reads `D$src$select` with `D$filter` needs no change
## (xtx_design_with()).

xtx_frame_root <- function(fr) fr$base %||% fr$src
xtx_frame_nodes <- function(fr) fr$nodes %||% list()

## The stem of the column names of a new node.
xtx_node_stem <- function(prefix) paste0("__xtx_", prefix, substr(xtx_tmp_name("n"), 3L, 10L))

## A node; `on`: the nodes of the frame it is made on (those its SQL names
## become its dependencies).
xtx_node <- function(kind, id, cols, on, items, join = "", ref = NULL) {
    list(id = id, kind = kind, cols = cols, deps = names(xtx_nodes_used(on, c(items, join))),
         items = paste(items, collapse = ", "), join = join, ref = ref)
}

## Frame `fr` with one more node.
xtx_frame_add <- function(fr, node) {
    xtx_frame_set(fr, xtx_frame_root(fr), c(xtx_frame_nodes(fr), stats::setNames(list(node), node$id)))
}

## Frame `fr` on root `root` with the nodes `nodes`; with `render`, as a
## design: the nodes rendered in its `src`.
xtx_frame_set <- function(fr, root, nodes = list(), render = FALSE) {
    fr$base <- NULL
    fr$nodes <- NULL
    fr$src <- root
    if (!length(nodes)) return(fr)
    fr$nodes <- nodes
    if (render) {
        fr$src <- xtx_root_render(root, nodes)
        fr$base <- root
    }
    fr
}

## The relation of a root with nodes on top: one derived projection per
## node that the root does not carry.
xtx_root_render <- function(root, nodes) {
    nodes <- nodes[!names(nodes) %in% root$carried]
    if (!length(nodes)) return(root)
    sel <- root$select
    for (nd in nodes) sel <- sprintf("SELECT s.*, %s FROM (%s) AS s%s", nd$items, sel, nd$join)
    root$select <- sel
    ## (not a table that batches read by rowid)
    root$base_table <- NULL
    root
}

## The nodes that SQL fragments name (by the stems of their columns), with
## the nodes those read.
xtx_nodes_used <- function(nodes, sql) {
    if (!length(nodes)) return(nodes)
    txt <- paste(sql, collapse = " ")
    keep <- vapply(nodes, function(nd) grepl(nd$id, txt, fixed = TRUE), NA)
    if (all(keep)) return(nodes)
    ids <- names(nodes)
    ## (dependencies come first: one pass from the end)
    for (i in rev(seq_along(nodes))) if (keep[i]) keep[ids %in% nodes[[i]]$deps] <- TRUE
    nodes[keep]
}

## A relation stored from root `parent` (a table, or the data of a design
## with columns added): `root` is its source record; it descends from the
## parent and carries the nodes `nodes` that were rendered into it.
xtx_root_child <- function(root, parent, nodes = list()) {
    xtx_tmp$roots <- (xtx_tmp$roots %||% 0) + 1
    root$path <- c(parent$path, sprintf("t%x", as.integer(xtx_tmp$roots)))
    root$carried <- unique(c(parent$carried, names(nodes)))
    root
}

## Of two roots of one design, the one that holds what both hold: the
## descendant. Two relations stored separately (neither descends from the
## other) have no rows in common order and no key: an error.
xtx_root_later <- function(a, b) {
    pa <- a$path %||% character()
    pb <- b$path %||% character()
    la <- length(pa)
    lb <- length(pb)
    if (la == lb && (la == 0L || identical(pa, pb))) return(a)
    if (la < lb && identical(pa, pb[seq_len(la)])) return(b)
    if (lb < la && identical(pb, pa[seq_len(lb)])) return(a)
    stop("Columns built on two different branches of the same remote design cannot be combined: ",
         "each reads a table stored separately from the design (col_compute()), and ",
         "remote rows have no key to align them. Combine the columns before storing them (store one ",
         "column made from both, or one product `D %*% cbind(b1, b2)` and then `[, 1]` and `[, 2]`), or ",
         "read them with as.data.frame(keep = ) and match them in R.",
         call. = FALSE)
}

## The frame on which the columns of frames a and b both evaluate: the
## later of their roots, with the nodes of both (those of a first, each
## once). They must be frames of the same design. The result has the
## fields of `a`.
xtx_frame_merge <- function(a, b) {
    if (is.null(a)) return(b)
    if (is.null(b)) return(a)
    if (!identical(a$id, b$id))
        stop("Columns of different remote designs cannot be combined (a subset or a materialized design ",
             "is another design: make the columns again on it).", call. = FALSE)
    ra <- xtx_frame_root(a)
    root <- xtx_root_later(ra, xtx_frame_root(b))
    na <- xtx_frame_nodes(a)
    nb <- xtx_frame_nodes(b)
    new <- if (length(nb)) nb[!names(nb) %in% names(na)]
    if (!length(new) && length(root$path) == length(ra$path)) return(a)
    xtx_frame_set(a, root, c(na, new))
}

## The frame of a read: the frame with the nodes that the SQL of the read
## (and the filter of the design) name rendered in its source, so that the
## query is `SELECT ... FROM (fr$src$select) AS s WHERE fr$filter`.
xtx_frame_ready <- function(fr, sql = character()) {
    if (is.null(fr)) return(fr)
    nodes <- xtx_frame_nodes(fr)
    if (!length(nodes)) return(fr)
    xtx_frame_set(fr, xtx_frame_root(fr), xtx_nodes_used(nodes, c(sql, fr$filter)), render = TRUE)
}

## The frame of a column with only the nodes its SQL (and filter) name.
xtx_column_frame1 <- function(x) {
    u <- unclass(x)
    fr <- u$frame
    nodes <- xtx_frame_nodes(fr)
    if (!length(nodes)) return(fr)
    used <- xtx_nodes_used(nodes, c(u$sql, fr$filter))
    if (length(used) == length(nodes)) fr else xtx_frame_set(fr, xtx_frame_root(fr), used)
}

## The frame of the remote columns in `cols`: a column, a xtx_c_codes()
## grouping, NULL, or a list of those (nested); NULL when none is remote.
xtx_cols_frame <- function(cols) {
    f <- NULL
    walk <- function(v) {
        if (xtx_is_column(v)) f <<- xtx_frame_merge(f, xtx_column_frame1(v))
        else if (xtx_is_codes(v)) walk(v$codes)
        else if (is.list(v) && is.null(oldClass(v))) for (e in v) walk(e)
        invisible()
    }
    walk(cols)
    f
}

## `cols` (as xtx_cols_frame() takes it) with fun() applied to each column.
xtx_cols_map <- function(v, fun) {
    if (xtx_is_column(v)) return(fun(v))
    if (xtx_is_codes(v)) {
        v$codes <- xtx_cols_map(v$codes, fun)
        return(v)
    }
    if (is.list(v) && is.null(oldClass(v))) return(lapply(v, xtx_cols_map, fun))
    v
}

xtx_column_frame <- function(...) {
    f <- NULL
    for (v in list(...)) if (xtx_is_column(v)) f <- xtx_frame_merge(f, unclass(v)$frame)
    f
}

## The identity of the values of a column: its design, its source relation
## and its expressions, which name their nodes by ids unique in the session
## (the same SQL text over another design or another stored relation is
## another column).
xtx_column_key <- function(x) {
    u <- unclass(x)
    list(id = u$frame$id, root = xtx_frame_root(u$frame)$path %||% character(), sql = u$sql)
}

## SQL of an operand (a remote column, or a number)
xtx_column_operand <- function(x) {
    if (xtx_is_column(x)) return(unclass(x)$sql)
    ## (NA is NULL; NaN stays NaN, as in R arithmetic)
    if (length(x) == 1L && is.na(x) && !is.nan(x)) return("CAST(NULL AS DOUBLE)")
    if (is.logical(x) && length(x) == 1L) return(if (x) "TRUE" else "FALSE")
    if (is.numeric(x) && length(x) == 1L) return(xtx_sql_double(x))
    stop("Remote columns combine only with numbers and other remote columns.", call. = FALSE)
}

#' @export
Ops.xtx_column <- function(e1, e2) {
    op <- .Generic
    if (nargs() == 1L) {
        u <- unclass(e1)
        if (op == "!") {
            ## (the negation of a negation is the column negated: !is.na(x)
            ## is the short SQL of "present", not NOT (NOT ...))
            if (!is.null(u$neg)) return(u$neg)
            lg <- xtx_column_lgl(e1)
            ## (a number is FALSE where it is 0, as in R: !x is x == 0, NA
            ## where x is NA or NaN)
            sql <- ifelse(lg, sprintf("(NOT %s)", u$sql), sprintf("(NULLIF(%s, CAST('NaN' AS DOUBLE)) = 0)", u$sql))
            out <- unclass(xtx_column(u$frame, sql, notnull = xtx_operand_notnull(e1) & lg, lgl = TRUE))
            ## (only the negation of a logical column is that column again)
            if (all(lg)) out$neg <- e1
            return(structure(out, class = "xtx_column"))
        }
        a <- xtx_column_arith(e1)
        sql <- switch(op, "-" = sprintf("(-%s)", a), "+" = a,
                      stop(sprintf("`%s` is not supported for remote columns.", op), call. = FALSE))
        return(xtx_column(u$frame, sql))
    }
    frame <- xtx_column_frame(e1, e2)
    ## a column of text against a string (the response of a remote frame
    ## that is a factor: design_frame()): equal or not, NULL where the
    ## column is
    if (op %in% c("==", "!=") && (is.character(e1) || is.character(e2))) {
        str <- if (is.character(e1)) e1 else e2
        if (length(str) != 1L || is.na(str)) stop("A remote column compares with one string.", call. = FALSE)
        col <- unclass(if (is.character(e1)) e2 else e1)$sql
        return(xtx_column(frame, sprintf("(%s %s %s)", col, if (op == "==") "=" else "<>",
                                         DBI::dbQuoteString(xtx_con(frame), str)), lgl = TRUE))
    }
    a <- xtx_column_operand(e1)
    b <- xtx_column_operand(e2)
    k <- max(length(a), length(b))
    if (!(length(a) %in% c(1L, k)) || !(length(b) %in% c(1L, k)))
        stop("Remote columns with different numbers of columns.", call. = FALSE)
    a <- rep_len(a, k)
    b <- rep_len(b, k)
    ## Comparisons with NaN are NA in R, but true or false in DuckDB (which
    ## orders NaN above every number and has NaN = NaN): compare with NaN
    ## made NULL (NULLIF(v, NaN), as NaN = NaN only for NaN, keeps v once in
    ## the expression), so that NaN or NULL on either side gives NULL.
    if (op %in% c("==", "!=", "<", ">", "<=", ">=")) {
        ## (two logical operands compare as they are; a logical column
        ## against a number is 0 or 1, as in R)
        la <- rep_len(xtx_operand_lgl(e1), k)
        lb <- rep_len(xtx_operand_lgl(e2), k)
        nn <- function(v, e, cast) {
            if (xtx_is_column(e))
                ## (a value stored as text or integer is never NaN)
                return(ifelse(cast, sprintf("CAST(%s AS DOUBLE)", v),
                              ifelse(rep_len(!xtx_column_type(e) %in% c("", "double"), k), v,
                                     sprintf("NULLIF(%s, CAST('NaN' AS DOUBLE))", v))))
            if (is.nan(e)) return("CAST(NULL AS DOUBLE)")
            if (is.logical(e) && !is.na(e)) return(ifelse(la & lb, v, xtx_sql_double(as.numeric(e))))
            v
        }
        a <- nn(a, e1, la & !lb)
        b <- nn(b, e2, lb & !la)
        ## (two logical columns: no NaN to guard)
        both <- la & lb
        if (any(both)) {
            if (xtx_is_column(e1)) a[both] <- rep_len(xtx_column_operand(e1), k)[both]
            if (xtx_is_column(e2)) b[both] <- rep_len(xtx_column_operand(e2), k)[both]
        }
    }
    if (op %in% c("&", "|")) {
        ## (a number is TRUE where it is not 0, as in R)
        tr <- function(v, e) {
            if (!xtx_is_column(e)) return(if (is.numeric(e) && !is.na(e)) (if (e != 0) "TRUE" else "FALSE") else v)
            ifelse(rep_len(xtx_column_lgl(e), k), v, sprintf("(NULLIF(%s, CAST('NaN' AS DOUBLE)) <> 0)", v))
        }
        a <- tr(a, e1)
        b <- tr(b, e2)
    }
    if (op %in% c("+", "-", "*", "/", "^")) {
        ## (a logical column is 0/1 in arithmetic, as in R)
        a <- rep_len(xtx_column_arith(e1), k)
        b <- rep_len(xtx_column_arith(e2), k)
    }
    sql <- switch(op,
        "+" = , "-" = , "*" = , "/" = sprintf("(%s %s %s)", a, op, b),
        "^" = xtx_sql_pow(a, b, e2),
        "==" = sprintf("(%s = %s)", a, b), "!=" = sprintf("(%s <> %s)", a, b),
        "<" = , ">" = , "<=" = , ">=" = sprintf("(%s %s %s)", a, op, b),
        "&" = sprintf("(%s AND %s)", a, b), "|" = sprintf("(%s OR %s)", a, b),
        stop(sprintf("`%s` is not supported for remote columns.", op), call. = FALSE))
    ## (AND, OR of columns that are never NULL are never NULL; the sums of a
    ## column known never NULL need no count of NULL values)
    nn <- op %in% c("&", "|") & rep_len(xtx_operand_notnull(e1), k) & rep_len(xtx_operand_notnull(e2), k)
    xtx_column(frame, sql, notnull = nn, lgl = op %in% c("==", "!=", "<", ">", "<=", ">=", "&", "|"))
}

## Is an operand logical-valued (per column)? A logical column, TRUE, FALSE.
xtx_operand_lgl <- function(x) if (xtx_is_column(x)) xtx_column_lgl(x) else is.logical(x) && length(x) == 1L && !is.na(x)

## A column marked as logical (comparisons, &, |, !, is.na() and the like):
## SQL has no arithmetic on BOOLEAN, so that such a column is cast to DOUBLE
## when it enters arithmetic (xtx_column_arith()), as R takes TRUE as 1.
xtx_column_logical <- function(x) {
    u <- unclass(x)
    u$lgl <- rep(TRUE, length(u$sql))
    structure(u, class = "xtx_column")
}

## The operand of an arithmetic operator: xtx_column_operand(), with logical
## columns (and TRUE, FALSE) as numbers.
xtx_column_arith <- function(x) {
    if (!xtx_is_column(x)) return(xtx_column_operand(if (is.logical(x) && length(x) == 1L && !is.na(x)) as.numeric(x) else x))
    u <- unclass(x)
    ifelse(xtx_column_lgl(x), sprintf("CAST(%s AS DOUBLE)", u$sql), u$sql)
}

## Is an operand never NULL (per column)? A remote column: its flags; a
## number: when it is not NA.
xtx_operand_notnull <- function(x) {
    if (xtx_is_column(x)) return(unclass(x)$notnull %||% rep(FALSE, length(unclass(x)$sql)))
    !is.na(x)
}

## is.na() and friends, with R's semantics (NA is NULL remotely): is.na
## TRUE for NA or NaN, is.nan for NaN only, is.finite FALSE for NA, NaN and
## +-Inf, is.infinite TRUE for +-Inf only. Never NA.
#' @export
is.na.xtx_column <- function(x) !xtx_c_present(x)

#' @export
is.nan.xtx_column <- function(x)
    xtx_column(unclass(x)$frame, sprintf("COALESCE(isnan(CAST(%s AS DOUBLE)), FALSE)", unclass(x)$sql),
               notnull = TRUE, lgl = TRUE)

#' @export
is.finite.xtx_column <- function(x) xtx_c_is_finite(x)

#' @export
is.infinite.xtx_column <- function(x)
    xtx_column(unclass(x)$frame, sprintf("COALESCE(isinf(CAST(%s AS DOUBLE)), FALSE)", unclass(x)$sql),
               notnull = TRUE, lgl = TRUE)

#' @export
Math.xtx_column <- function(x, ...) {
    if (.Generic == "log") return(xtx_column_log(x, ...))
    ## (log2 and log10 as R: -Inf at 0, NaN below, where DuckDB refuses)
    if (.Generic == "log2") return(xtx_column_log(x, 2))
    if (.Generic == "log10") return(xtx_column_log(x, 10))
    ## (a logical column is 0/1)
    a <- xtx_column_arith(x)
    sql <- switch(.Generic,
        exp = sprintf("exp(%s)", a), sqrt = sprintf("sqrt(%s)", a),
        abs = sprintf("abs(%s)", a), log1p = xtx_sql_log1p(sprintf("CAST(%s AS DOUBLE)", a)), expm1 = {
            ## (accurate near zero: the macro xtx_expm1)
            xtx_sql_macro_expm1(unclass(x)$frame)
            sprintf("xtx_expm1(CAST(%s AS DOUBLE))", a)
        },
        cos = sprintf("cos(%s)", a), sin = sprintf("sin(%s)", a), tan = sprintf("tan(%s)", a),
        floor = sprintf("floor(%s)", a), ceiling = sprintf("ceil(%s)", a),
        ## (DuckDB's sign(NaN) is 0, R's NaN)
        sign = sprintf("CASE WHEN isnan(CAST(%s AS DOUBLE)) THEN CAST('NaN' AS DOUBLE) ELSE sign(%s) END", a, a),
        ## (round, signif, trunc, gamma, lgamma and the errors: R/columns_public.R)
        xtx_math_sql(x, .Generic, ...))
    xtx_column(unclass(x)$frame, sql)
}

## log(x) and log(x, base) of a remote column, as R computes them (nm8t):
## R_log(v) is ln(v) for v > 0, -Inf for 0 and NaN for v < 0 (DuckDB's ln()
## refuses 0 and negative numbers: it is reached only for v > 0, NaN or
## NULL, which it propagates); log(x, base) is log10(x) where base is 10,
## log2(x) where it is 2, else R_log(x) / R_log(base) (R's logbase()), so
## that base 1 gives x / 0 (Inf, -Inf or NaN), base 0 R_log(x) / -Inf and a
## negative base NaN; NA (NULL) in x or base gives NA, else NaN gives NaN.
## `base`: a number (its R_log() computed in R) or a column of the design
## (one column, or as many as x).
xtx_column_log <- function(x, base) {
    xv <- sprintf("CAST(%s AS DOUBLE)", unclass(x)$sql)
    if (missing(base)) return(xtx_column(unclass(x)$frame, xtx_sql_rlog(xv)))
    if (xtx_is_column(base)) {
        f <- xtx_column_frame(x, base)
        bv <- sprintf("CAST(%s AS DOUBLE)", unclass(base)$sql)
        k <- max(length(xv), length(bv))
        if (!(length(xv) %in% c(1L, k)) || !(length(bv) %in% c(1L, k)))
            stop("log() of remote columns with different numbers of columns.", call. = FALSE)
        xv <- rep_len(xv, k)
        bv <- rep_len(bv, k)
        sql <- sprintf("(CASE WHEN %1$s = 10 THEN %2$s WHEN %1$s = 2 THEN %3$s ELSE %4$s / %5$s END)",
                       bv, xtx_sql_rlog(xv, "log10"), xtx_sql_rlog(xv, "log2"), xtx_sql_rlog(xv), xtx_sql_rlog(bv))
        return(xtx_column(f, sql))
    }
    if (!(is.numeric(base) || is.logical(base)) || length(base) != 1L || !is.null(dim(base)))
        stop("log() of a remote column: `base` must be one number or a column of the same design.", call. = FALSE)
    xtx_column(unclass(x)$frame, xtx_sql_log(xv, base))
}

## R_log(v) in SQL, `v` a DOUBLE operand (named three times: a short one);
## `fun` ln, log2 or log10.
xtx_sql_rlog <- function(v, fun = "ln")
    sprintf("(CASE WHEN %1$s = 0 THEN CAST('-Infinity' AS DOUBLE) WHEN %1$s < 0 THEN CAST('NaN' AS DOUBLE) ELSE %2$s(%1$s) END)",
            v, fun)

## log(v, base) in SQL as R computes it (see xtx_column_log()), `v` a DOUBLE
## operand, `base` one number or logical.
xtx_sql_log <- function(v, base) {
    base <- as.double(base)
    if (!is.na(base) && base == 10) return(xtx_sql_rlog(v, "log10"))
    if (!is.na(base) && base == 2) return(xtx_sql_rlog(v, "log2"))
    ## R_log(base), in R (NA is NULL: the result is NA)
    den <- if (is.na(base) && !is.nan(base)) "CAST(NULL AS DOUBLE)"
           else xtx_sql_double(if (is.nan(base)) NaN else if (base > 0) log(base) else if (base == 0) -Inf else NaN)
    sprintf("(%s / %s)", xtx_sql_rlog(v), den)
}

## log1p(v) in SQL, `v` a DOUBLE operand (named several times: a short
## one). DuckDB has no log1p(), and ln(1 + v) loses the digits of v near 0:
## with u = 1 + v rounded, ln(u) * v / (u - 1) is accurate to a few ulps
## (Goldberg 1991, Theorem 4), and v itself where u rounds to 1. As R: -Inf
## at -1, NaN below (and at -Inf), Inf at Inf, NaN for NaN, NA for NA.
xtx_sql_log1p <- function(v)
    sprintf(paste("(CASE WHEN %1$s = -1 THEN CAST('-Infinity' AS DOUBLE)",
                  "WHEN %1$s < -1 THEN CAST('NaN' AS DOUBLE)",
                  "WHEN 1 + %1$s = 1 OR isinf(%1$s) THEN %1$s",
                  "ELSE ln(1 + %1$s) * %1$s / ((1 + %1$s) - 1) END)"), v)

## x^y in SQL (a, b the SQL of the operands, e2 the R operand y), with R's
## R_pow() where DuckDB's pow() differs: x^0 and 1^y are 1 for any x and y
## (NA and NaN included; pow(NULL, 0) is NULL), and a negative base with
## an infinite power, or -Inf with a power that is not an integer, gives
## NaN. A power that is a nonzero integer (as y^2) needs none of this.
xtx_sql_pow <- function(a, b, e2) {
    if (!xtx_is_column(e2) && is.numeric(e2) && length(e2) == 1L && is.finite(e2) && e2 != 0 && e2 == round(e2))
        return(sprintf("pow(%s, %s)", a, b))
    sprintf(paste("(CASE WHEN %2$s = 0 OR %1$s = 1 THEN CAST(1 AS DOUBLE)",
                  "WHEN %1$s < 0 AND (isinf(CAST(%2$s AS DOUBLE)) OR (isinf(CAST(%1$s AS DOUBLE)) AND %2$s <> floor(%2$s)))",
                  "THEN CAST('NaN' AS DOUBLE) ELSE pow(%1$s, %2$s) END)"), a, b)
}

## pmax/pmin with R's missing values: DuckDB's greatest() and least() skip
## NULL and order NaN above every number, while in R an NA or NaN gives NA or
## NaN (that of y where both are missing, as R's pmax takes the later
## missing argument). NULLIF(v, NaN) is NULL where v is NA or NaN. With a
## finite number on one side, the column is kept where it is missing or
## beyond the number. (The fits use xtx_fast_pmax/xtx_fast_pmin, below.)
xtx_c_pmax <- function(x, y) {
    f <- xtx_column_frame(x, y)
    if (is.null(f)) return(pmax(x, y))
    xtx_sql_pext(f, x, y, "greatest", "<=")
}

xtx_c_pmin <- function(x, y) {
    f <- xtx_column_frame(x, y)
    if (is.null(f)) return(pmin(x, y))
    xtx_sql_pext(f, x, y, "least", ">=")
}

xtx_sql_pext <- function(f, x, y, fun, cmp) {
    num <- function(v) !xtx_is_column(v) && is.numeric(v) && length(v) == 1L && is.finite(v)
    a <- xtx_column_arith(x)
    b <- xtx_column_arith(y)
    nan <- "CAST('NaN' AS DOUBLE)"
    if (num(x) || num(y)) {
        if (num(x)) { cc <- a; v <- b } else { cc <- b; v <- a }
        return(xtx_column(f, sprintf("CASE WHEN NULLIF(%1$s, %3$s) %4$s %2$s THEN %2$s ELSE %1$s END",
                                   v, cc, nan, cmp)))
    }
    xtx_column(f, sprintf(paste("CASE WHEN NULLIF(%2$s, %3$s) IS NULL THEN %2$s",
                              "WHEN NULLIF(%1$s, %3$s) IS NULL THEN %1$s ELSE %4$s(%1$s, %2$s) END"),
                        a, b, nan, fun))
}

## Internal pmax/pmin for the fits (qpr8), only where the inputs are never
## missing: remotely DuckDB's greatest() and least(), which name each
## operand once. The CASE of xtx_c_pmax/xtx_c_pmin repeats its operand, and the
## links nest them around eta inside the pnorm/dnorm macros, which repeat
## their argument again (a small probit fit: 1.3 s -> 3.4 s). They give
## pmax/pmin wherever no input is NA or NaN (greatest() skips NULL; least()
## skips NULL and NaN); each use states why its inputs are not. In memory:
## pmax/pmin.
xtx_fast_pmax <- function(x, y) {
    f <- xtx_column_frame(x, y)
    if (is.null(f)) return(pmax(x, y))
    xtx_column(f, sprintf("greatest(%s, %s)", xtx_column_arith(x), xtx_column_arith(y)))
}

xtx_fast_pmin <- function(x, y) {
    f <- xtx_column_frame(x, y)
    if (is.null(f)) return(pmin(x, y))
    xtx_column(f, sprintf("least(%s, %s)", xtx_column_arith(x), xtx_column_arith(y)))
}

## Internal log for the fits (qpr8), only where x is positive by
## construction (or NULL): remotely DuckDB's ln(), which names x once (the
## CASE of log(), R's log of 0 and of negative numbers, names it three
## times; at 0 and below, ln() gives -Inf and NaN in DuckDB 2.0 but older
## versions refuse them). In memory: log().
xtx_fast_log <- function(x) {
    if (!xtx_is_column(x)) return(log(x))
    xtx_column(unclass(x)$frame, sprintf("ln(CAST(%s AS DOUBLE))", unclass(x)$sql))
}

## (in memory a numeric vector: `no` with `yes` where test is TRUE and NA
## where it is NA, as ifelse() but without its copies and attributes)
xtx_c_ifelse <- function(test, yes, no) {
    f <- xtx_column_frame(test, yes, no)
    if (is.null(f)) {
        n <- length(test)
        ## (one pass in C, which allocates only the result)
        if (is.logical(test) && n > 0L && length(yes) %in% c(1L, n) && length(no) >= 1L &&
            is.numeric(yes) && is.numeric(no)) {
            if (!is.double(yes)) yes <- as.double(yes)
            if (!is.double(no)) no <- as.double(no)
            return(.Call(C_xtx_ifelse, test, yes, no))
        }
        out <- rep_len(as.double(no), n)
        i <- which(test)
        out[i] <- if (length(yes) == 1L) yes else yes[i]
        if (anyNA(test)) out[is.na(test)] <- NA
        return(out)
    }
    ## (NULL where test is NULL, as NA in memory; the simple CASE names the
    ## test once; never NULL when no operand is)
    nn <- all(xtx_operand_notnull(test)) && all(xtx_operand_notnull(yes)) && all(xtx_operand_notnull(no))
    br <- xtx_branches(yes, no)
    xtx_column(f, sprintf("(CASE %s WHEN TRUE THEN %s WHEN FALSE THEN %s END)", xtx_column_test(test),
                        br$yes, br$no), notnull = nn, lgl = br$lgl)
}

## The branches of an ifelse (columns or numbers) as SQL of one type:
## logical when both are, else numbers (a logical branch as 0/1, as R
## coerces). list(yes, no, lgl).
xtx_branches <- function(yes, no) {
    lgl <- all(xtx_operand_lgl(yes)) && all(xtx_operand_lgl(no))
    if (lgl) list(yes = xtx_column_operand(yes), no = xtx_column_operand(no), lgl = TRUE)
    else list(yes = xtx_column_arith(yes), no = xtx_column_arith(no), lgl = FALSE)
}

## The test of an ifelse as a BOOLEAN: a logical column as it is, a number
## TRUE where it is not 0 (NULL where it is NA or NaN).
xtx_column_test <- function(test) {
    if (!xtx_is_column(test)) return(xtx_column_operand(test))
    ifelse(xtx_column_lgl(test), unclass(test)$sql,
           sprintf("(NULLIF(%s, CAST('NaN' AS DOUBLE)) <> 0)", unclass(test)$sql))
}

## Columns side by side (cbind); each keeps its flags.
xtx_c_cbind <- function(...) {
    args <- list(...)
    f <- do.call(xtx_column_frame, args)
    if (is.null(f)) return(do.call(cbind, args))
    ## (the expressions of columns as they are: already parenthesized)
    out <- xtx_column(f, unlist(lapply(args, function(a) if (xtx_is_column(a)) "" else xtx_column_operand(a))),
             notnull = unlist(lapply(args, function(a) if (xtx_is_column(a)) unclass(a)$notnull else !is.na(a))),
             lgl = unlist(lapply(args, xtx_operand_lgl)),
             type = unlist(lapply(args, function(a) if (xtx_is_column(a)) xtx_column_type(a) else "")))
    out <- unclass(out)
    out$sql <- unlist(lapply(args, function(a) if (xtx_is_column(a)) unclass(a)$sql else sprintf("(%s)", xtx_column_operand(a))))
    structure(out, class = "xtx_column")
}

xtx_c_rowsums <- function(x) {
    if (!xtx_is_column(x)) return(if (is.matrix(x)) rowSums(x) else x)
    xtx_column(unclass(x)$frame, sprintf("(%s)", paste(xtx_column_arith(x), collapse = " + ")))
}

## distribution functions: stats in memory, DuckDB macros remotely
xtx_c_macro <- function(x, fun, macro) {
    if (!xtx_is_column(x)) return(fun(x))
    xtx_sql_macros(unclass(x)$frame)
    xtx_column(unclass(x)$frame, sprintf("%s(%s)", macro, xtx_column_arith(x)))
}
xtx_c_pnorm <- function(x) xtx_c_macro(x, stats::pnorm, "xtx_pnorm")
xtx_c_qnorm <- function(x) xtx_c_macro(x, stats::qnorm, "xtx_qnorm")
xtx_c_dnorm <- function(x) xtx_c_macro(x, stats::dnorm, "xtx_dnorm")
xtx_c_pcauchy <- function(x) xtx_c_macro(x, stats::pcauchy, "xtx_pcauchy")
xtx_c_mexpm1m <- function(x) xtx_c_macro(x, function(u) -expm1(-u), "xtx_mexpm1m")

## The values of a column in R (remotely: all its rows; for small data).
xtx_p_collect <- function(x) {
    if (!xtx_is_column(x)) return(x)
    sql <- unclass(x)$sql
    fr <- xtx_frame_ready(unclass(x)$frame, sql)
    res <- xtx_query(fr, sprintf("SELECT %s FROM (%s) AS s WHERE %s",
                                paste(sprintf("CAST(%s AS DOUBLE) AS c%d", sql, seq_along(sql)), collapse = ", "),
                                fr$src$select, fr$filter))
    if (length(sql) == 1L) res[[1L]] else as.matrix(res)
}

## ---- sums over rows ----

xtx_colsum <- function(x) {
    if (!xtx_is_column(x)) {
        if (is.matrix(x)) return(colSums(x))
        return(sum(x))
    }
    sql <- unclass(x)$sql
    fr <- xtx_frame_ready(unclass(x)$frame, sql)
    k <- length(sql)
    nn <- unclass(x)$notnull %||% rep(FALSE, k)
    items <- c(sprintf("SUM(CAST(%s AS DOUBLE)) AS s%d", sql, seq_len(k)),
               sprintf("SUM(CASE WHEN (%s) IS NULL THEN 1 ELSE 0 END) AS n%d", sql, seq_len(k))[!nn])
    res <- xtx_query(fr, sprintf("SELECT %s FROM (%s) AS s WHERE %s", paste(items, collapse = ", "),
                                fr$src$select, fr$filter))
    vapply(seq_len(k), function(j) {
        v <- as.double(res[[paste0("s", j)]])
        if (!nn[j] && as.numeric(res[[paste0("n", j)]]) > 0) return(NA_real_)
        if (is.na(v)) 0 else v
    }, 0)
}

## Sums over rows of several columns (vectors or numbers in memory; one
## query remotely): a number per argument.
xtx_c_sums <- function(...) {
    a <- list(...)
    if (any(vapply(a, xtx_is_column, NA))) return(xtx_colsum(do.call(xtx_c_cbind, a)))
    vapply(a, function(v) sum(v), 0)
}

## "sql" for a remote design, "" in memory.
xtx_engine <- function(x) if (inherits(x, "xtx_remote")) "sql" else ""

xtx_nrows <- function(x) {
    if (inherits(x, "xtx_remote")) return(x$cache$n %||% x$n)
    if (xtx_is_column(x)) return(xtx_nrows(unclass(x)$frame))
    NROW(x)
}

xtx_p_colnorms <- function(x) {
    if (!xtx_is_column(x)) return(.Call(C_xtx_col_norms, if (is.matrix(x)) x else matrix(x, ncol = 1L)))
    sqrt(xtx_colsum(x * x))
}

## ---- designs ----

xtx_p_column <- function(x, sql) {
    if (is.null(sql)) return(NULL)
    if (inherits(x, "xtx_remote")) return(xtx_column(x, sql))
    sql
}

## The remote design with the nodes of the columns `cols` (a design, or
## columns as xtx_cols_frame() takes them: those their SQL names) rendered
## in its source, on the later of its root and theirs, so that the SQL of
## the columns, the scale and the filter of the design all read
## `x$src$select` (same plan, rows and scale). The columns must be columns
## of the design. In memory x.
xtx_design_with <- function(x, cols, short = TRUE) {
    if (!inherits(x, "xtx_remote")) return(x)
    fr <- if (inherits(cols, "xtx_remote")) cols else xtx_cols_frame(cols)
    if (is.null(fr)) return(x)
    if (short && !identical(x$id, fr$id)) stop("Columns of another remote design.", call. = FALSE)
    m <- xtx_frame_merge(x, fr)
    xtx_frame_set(x, xtx_frame_root(m), xtx_frame_nodes(m), render = TRUE)
}

xtx_p_scale <- function(x, s) {
    if (is.null(s)) return(x)
    if (!xtx_is_column(s)) return(s * x)
    x <- xtx_design_with(x, s)
    ## (a logical column scales by 0 or 1)
    x$scale <- xtx_remote_scale_sql(x, DBI::SQL(xtx_column_arith(s)), "*")
    x
}

## A remote design with the columns `col` added to its source (named nm):
## its source is a relation derived from the design's (a new root, which
## carries the nodes of the columns), not a table.
xtx_with_columns <- function(x, col, nm) {
    x <- xtx_design_with(x, col, short = FALSE)
    ext <- xtx_remote_extend(xtx_remote_unscaled(x), stats::setNames(unclass(col)$sql, nm), "column")
    xtx_frame_set(x, xtx_root_child(ext$src, xtx_frame_root(x), xtx_frame_nodes(x)))
}

xtx_p_crossprod <- function(x, y = NULL, compress = TRUE) {
    ## A refinement step reads the residual column once. Skip the pattern
    ## estimate and aggregate its products directly.
    if (!compress && inherits(x, "xtx_remote") && xtx_is_column(y)) {
        yc <- xtx_remote_cols(x, y, "y")
        return(xtx_remote_moments(xtx_remote_unscaled(yc$D), yc$sql, yc$D$scale,
                                 gram = FALSE, compress = FALSE)$xty)
    }
    if (is.null(y)) return(crossprod(x))
    if (inherits(x, "xtx_remote") && !inherits(y, "xtx_remote") && !xtx_is_column(y))
        stop("`y` must be a design or a column.", call. = FALSE)
    crossprod(x, y)
}

xtx_p_colsums <- function(x, group = NULL) {
    if (is.null(group)) colSums(x) else colSums(x, group = group)
}

## SQL GROUP BY over the frame of columns: list(groups (a data frame of the
## group values, sorted as rowsum() sorts them), sums (matrix)).
xtx_group_query <- function(groups, values, nulls = FALSE) {
    gs <- vapply(groups, function(g) unclass(g)$sql, "")
    vs <- unlist(lapply(values, function(v) unclass(v)$sql))
    fr <- xtx_frame_ready(do.call(xtx_column_frame, c(groups, values)), c(gs, vs))
    items <- c(sprintf("%s AS g%d", gs, seq_along(gs)),
               if (length(vs)) sprintf("SUM(CAST(%s AS DOUBLE)) AS v%d", vs, seq_along(vs)) else "COUNT(*) AS v1",
               if (nulls && length(vs)) sprintf("SUM(CASE WHEN (%s) IS NULL THEN 1 ELSE 0 END) AS n%d", vs, seq_along(vs)))
    res <- xtx_query(fr, sprintf("SELECT %s FROM (%s) AS s WHERE %s GROUP BY ALL",
                                paste(items, collapse = ", "), fr$src$select, fr$filter))
    res
}

xtx_p_groupsum <- function(col, group) {
    if (xtx_is_codes(group)) return(xtx_p_groupsum_codes(col, group))
    if (!xtx_is_column(col) && !xtx_is_column(group)) {
        return(rowsum(if (is.matrix(col)) col else matrix(col, ncol = 1L), group, reorder = TRUE))
    }
    res <- xtx_group_query(list(group), list(col), nulls = TRUE)
    ug <- sort(unique(res$g1), na.last = TRUE, method = "quick")
    if (anyNA(ug)) warning("missing values for 'group'", call. = FALSE)
    i <- match(res$g1, ug)
    k <- xtx_ncols(col)
    out <- matrix(0, length(ug), k, dimnames = list(as.character(ug), NULL))
    for (j in seq_len(k)) {
        v <- as.double(res[[paste0("v", j)]])
        ## a group with a missing value sums to NA, as rowsum() gives
        v[as.numeric(res[[paste0("n", j)]]) > 0] <- NA_real_
        out[i, j] <- v
    }
    out
}

xtx_p_crosstab <- function(a, b, w = NULL) {
    if (xtx_is_codes(a) && xtx_is_codes(b)) return(xtx_p_crosstab_codes(a, b, w))
    if (!xtx_is_column(a) && !xtx_is_column(b) && !xtx_is_column(w)) {
        ua <- sort(unique(a), na.last = TRUE, method = "quick")
        ub <- sort(unique(b), na.last = TRUE, method = "quick")
        return(Matrix::sparseMatrix(i = match(a, ua), j = match(b, ub),
                                    x = if (is.null(w)) rep(1, length(a)) else as.double(w),
                                    dims = c(length(ua), length(ub)),
                                    dimnames = list(as.character(ua), as.character(ub))))
    }
    res <- xtx_group_query(list(a, b), if (is.null(w)) list() else list(w))
    ## the group values as in memory: sorted in R, labels as.character()
    ua <- sort(unique(res$g1), na.last = TRUE, method = "quick")
    ub <- sort(unique(res$g2), na.last = TRUE, method = "quick")
    Matrix::sparseMatrix(i = match(res$g1, ua), j = match(res$g2, ub), x = as.double(res$v1),
                         dims = c(length(ua), length(ub)),
                         dimnames = list(as.character(ua), as.character(ub)))
}

## (remotely a node: the SQL of the product over the rows of the design, its
## scale included, is written now and rendered by each read that names the
## column; nothing is created in the database)
xtx_p_mult <- function(x, B) {
    if (!inherits(x, "xtx_remote")) return(if (is.matrix(B)) x %*% B else drop(x %*% B))
    B <- if (is.matrix(B)) B else matrix(B, ncol = 1L)
    checkmate::assert_matrix(B, mode = "numeric", nrows = x$p, any.missing = FALSE, .var.name = "beta")
    checkmate::assert_true(all(is.finite(B)), .var.name = "all(is.finite(beta))")
    storage.mode(B) <- "double"
    stem <- xtx_node_stem("p")
    nm <- paste0(stem, "_", seq_len(ncol(B)))
    pr <- xtx_sql_predict(x, B)
    preds <- if (is.null(x$scale)) pr$preds else sprintf("(%s) * (%s)", pr$preds, x$scale)
    node <- xtx_node("mult", stem, nm, xtx_frame_nodes(x), sprintf("%s AS %s", preds, xtx_qi(x, nm)), pr$joins,
                     ref = list(design = x$id, B = B, scale = x$scale))
    xtx_column(xtx_frame_add(xtx_remote_unscaled(x), node), sprintf("s.%s", xtx_qi(x, nm)))
}

xtx_p_quad <- function(x, V) {
    if (!inherits(x, "xtx_remote")) {
        J <- as_jac(x)
        storage.mode(V) <- "double"
        if (!jac_kernel_ok(J)) {
            X <- as.matrix(J)
            return(unname(rowSums((X %*% V) * X)))
        }
        g <- xtx_dedup_groups(J, cost = xtx_row_cost(J)^2)
        ## A plain design within the bounds of the structured kernels: its
        ## dense columns through a matrix product, its factors through
        ## tables (xtx_struct_quad()). The row kernel below costs the
        ## square of the entries of each row, which is fine for a few
        ## factor levels and several times slower than the dense
        ## computation for numeric variables.
        ## With groups of identical rows, the same kernel on the distinct
        ## rows: each row's result is the same bit for bit either way.
        st <- xtx_struct_single(J)
        if (!is.null(st) && xtx_in_range(V, 1e-100, 1e100)) {
            Dq <- if (is.null(g)) st$design else st$design[g$first, ]
            q <- if (xtx_struct_ok(Dq)) xtx_struct_quad(Dq, V)
            if (!is.null(q)) {
                if (!is.null(g)) q <- q[g$map]
                if (!is.null(st$scale)) q <- st$scale^2 * q
                if (xtx_all_finite(q)) return(q)
            }
        }
        q <- do.call(.Call, c(list(C_xtx_jac_quad), jac_args(J),
                              list(V, if (is.null(g)) NULL else g$first, design_threads())))
        return(if (is.null(g)) q else q[g$map])
    }
    zero <- xtx_column(x, "CAST(0 AS DOUBLE)")
    if (ncol(V) == 0L) return(zero)
    ## x'Vx = x'Sx with S = (V + V') / 2 the symmetric part of V. With S =
    ## sum_j l_j u_j u_j', it is sum_j l_j (x'u_j)^2: the squares of one
    ## product X L, added for the positive eigenvalues and subtracted for
    ## the negative ones (V need not be positive semi-definite).
    e <- eigen((V + t(V)) / 2, symmetric = TRUE)
    pos <- which(e$values > 0)
    neg <- which(e$values < 0)
    keep <- c(pos, neg)
    if (!length(keep)) return(zero)
    L <- e$vectors[, keep, drop = FALSE] %*% diag(sqrt(abs(e$values[keep])), length(keep))
    P <- xtx_p_mult(x, L)
    P2 <- P * P
    if (!length(neg)) return(xtx_c_rowsums(P2))
    qn <- xtx_c_rowsums(P2[, length(pos) + seq_along(neg)])
    if (!length(pos)) return(-qn)
    xtx_c_rowsums(P2[, seq_along(pos)]) - qn
}

## Compressed least squares (see the header); the in-memory method groups
## identical rows (exports.R), the remote one is the scan of the SQL engine
## grouped by the cells of the design (remote_compress.R).
xtx_p_patterns <- function(x, z, w = NULL, compress = "auto", keys = NULL, fit = "lm") {
    if (inherits(x, "xtx_remote")) {
        mode <- if (isTRUE(compress)) TRUE else if (isFALSE(compress)) FALSE else xtx_compress_mode()
        ## (with the nodes of z and w)
        x <- xtx_design_with(x, list(z, w))
        ## (extra key columns: more keys of the GROUP BY)
        if (!is.null(keys)) x <- xtx_remote_cols(x, keys, "keys")$D
        pt <- xtx_sql_patterns(x, unclass(z)$sql, if (!is.null(w)) unclass(w)$sql,
                              keys = if (!is.null(keys)) unclass(keys)$sql, mode = mode, fit = fit)
        if (is.null(pt)) return(NULL)
        ## (the rows of the design, which later queries reuse: for a linear
        ## model, as its fit kept them)
        if (identical(fit, "lm")) assign("n", sum(pt$n), envir = x$cache)
        pos <- which(pt$sw > 0)
        ws <- pt$sw[pos]
        if (is.null(w) && all(ws == 1)) ws <- NULL
        return(list(x = xtx_dense_design(pt$x[pos, , drop = FALSE]), z = pt$sy[pos] / pt$sw[pos], w = ws,
                    within = sum(pt$within), rows = sum(pt$n),
                    keys = if (!is.null(keys)) pt$keys[pos, , drop = FALSE]))
    }
    g <- xtx_fit_groups(x, compress)
    if (is.null(g)) return(NULL)
    if (!is.null(keys)) {
        ## patterns within the values of the key columns
        K <- if (is.matrix(keys)) keys else matrix(keys, ncol = 1L)
        key <- do.call(paste, c(list(g$map), lapply(seq_len(ncol(K)), function(j) sprintf("%a", K[, j])), sep = "\r"))
        map <- match(key, unique(key))
        g <- list(map = map, first = match(seq_len(max(map)), map))
        if (length(g$first) == nrow(x) || (!isTRUE(compress) && length(g$first) > 0.5 * nrow(x))) return(NULL)
    }
    ss <- xtx_fit_suffstats(z, w, g$map, length(g$first))
    ## sums that overflow (w * z beyond the range of doubles): no compression
    if (!all(is.finite(ss$w), is.finite(ss$mean), is.finite(ss$within))) return(NULL)
    keep <- ss$w > 0
    ws <- ss$w[keep]
    if (is.null(w) && all(ws == 1)) ws <- NULL
    list(x = x[g$first[keep], ], z = ss$mean[keep], w = ws, within = ss$within, rows = nrow(x),
         keys = if (!is.null(keys)) (if (is.matrix(keys)) keys else matrix(keys, ncol = 1L))[g$first[keep], , drop = FALSE])
}

## Row-wise results: in memory memory(d), d the data frame of the columns
## (default: d); remotely a table (a temporary view with the columns of the
## remote data and `cols`), padded by napredict(na) for a fit with
## na.exclude, with attributes `attrs`; with design = TRUE, the remote
## design with the columns added to its data (so that they can scale its
## rows).
xtx_p_output <- function(x, cols, memory = identity, na = NULL, design = FALSE, attrs = list()) {
    if (!inherits(x, "xtx_remote")) return(memory(as.data.frame(lapply(cols, function(v) as.vector(v)))))
    taken <- intersect(names(cols), names(x$src$types))
    if (length(taken)) stop(sprintf("The remote data already has a column `%s`.", taken[1L]), call. = FALSE)
    if (design) return(xtx_with_columns(x, do.call(xtx_c_cbind, unname(cols)), names(cols)))
    fr <- xtx_design_with(xtx_remote_unscaled(x), unname(cols))
    srccols <- names(x$src$types)
    items <- vapply(names(cols), function(nm) {
        v <- cols[[nm]]
        sprintf("%s AS %s", if (xtx_is_column(v)) unclass(v)$sql else xtx_column_operand(v), xtx_qi(fr, nm))
    }, "")
    sel <- sprintf("SELECT %s FROM (%s) AS s WHERE %s",
                   paste(c(paste0("s.", xtx_qi(fr, srccols)), items), collapse = ", "), fr$src$select, fr$filter)
    name <- xtx_tmp_name("xtx_out")
    xtx_execute(fr, sprintf("CREATE TEMPORARY VIEW %s AS %s", xtx_qi(fr, name), sel))
    out <- stats::napredict(na, structure(name, class = "xtx_remote_table", connection = xtx_con(fr),
                                          holder = fr$src$holder))
    for (a in names(attrs)) attr(out, a) <- attrs[[a]]
    out
}

#' Row-wise results of a design
#'
#' `design_output()` returns results with one value per row of a design (such
#' as fitted values and residuals), on both backends: in memory
#' `memory(d)`, with `d` the data frame of the columns; for a remote design, a
#' table in the database (the name of a temporary view, of class
#' `xtx_remote_table`, with every column of the remote data and the columns
#' `cols`), padded by [napredict()] with `na` (for instance the `na.action`
#' of a remote fit of 'xtx.fit' with `na.exclude`), with the attributes `attrs`; with `design = TRUE`, the
#' remote design with the columns added to its data (so that they can scale
#' its rows: `crossprod(R, DBI::SQL('".resid"^2') * R)`).
#'
#' Names may change before the first release.
#'
#' @param x A design (`design_matrix` or remote design).
#' @param cols A named list of columns of `x` (vectors in memory,
#'   [design_column()]s remotely) or numbers.
#' @param memory In memory, a function of the data frame of `cols` that
#'   returns the result (default: the data frame).
#' @param na Remotely, the `na.action` whose [napredict()] method pads the
#'   table (`NULL`: none).
#' @param design Remotely, return the design with the columns added instead
#'   of a table.
#' @param attrs Remotely, attributes of the table returned.
#' @return In memory, `memory(d)`; remotely a `xtx_remote_table` (or, with
#'   `design = TRUE`, a remote design).
#' @examples
#' d <- data.frame(y = rnorm(10), x = rnorm(10))
#' D <- design_matrix(y ~ x, d)
#' f <- drop(D %*% c(0, 1))
#' head(design_output(D, list(.fitted = f, .resid = d$y - f)))
#' @keywords internal
#' @export
design_output <- function(x, cols, memory = identity, na = NULL, design = FALSE, attrs = list()) {
    checkmate::assert_multi_class(x, c("xtx_matrix", "xtx_remote"))
    checkmate::assert_list(cols, names = "unique")
    checkmate::assert_function(memory)
    checkmate::assert_flag(design)
    xtx_p_output(x, cols, memory = memory, na = na, design = design, attrs = attrs)
}

## DuckDB macros for the distribution functions of remote columns, as R
## computes them (from the GLM work, kata/402b-glm): xtx_pnorm, xtx_qnorm
## (one ulp), xtx_dnorm, xtx_pcauchy and xtx_mexpm1m(u) = -expm1(-u). Created
## once per connection.
## The macro of expm1() of a remote column, created at its first use on a
## connection (apart from the macros of the GLM links, so that fits that
## never call expm1() send no statement for it): Kahan's trick for either
## sign, (w - 1) x / ln(w) with w = exp(x); x itself when w rounds to 1, -1
## when w - 1 rounds to -1, exp(x) when it is large (and for NaN, which
## DuckDB orders above every number).
xtx_sql_macro_expm1 <- function(D) {
    if (isTRUE(D$src$holder$xtx_macro_expm1)) return(invisible())
    xtx_execute(D, paste("CREATE OR REPLACE TEMP MACRO xtx_expm1(x) AS (CASE WHEN x > 40 THEN exp(x)",
                         "WHEN exp(x) = 1 THEN x WHEN exp(x) - 1 = -1 THEN -1.0",
                         "ELSE (exp(x) - 1) * x / ln(exp(x)) END)"))
    D$src$holder$xtx_macro_expm1 <- TRUE
    invisible()
}

xtx_sql_macros <- function(D) {
    if (isTRUE(D$src$holder$xtx_macros)) return(invisible())
    horner <- function(x, co) {
        s <- sprintf("%.17g", co[1L])
        for (k in co[-1L]) s <- sprintf("(%s * %s + %.17g)", s, x, k)
        s
    }
    a <- c(2.2352520354606839287, 161.02823106855587881, 1067.6894854603709582, 18154.981253343561249,
           0.065682337918207449113)
    b <- c(47.20258190468824187, 976.09855173777669322, 10260.932208618978205, 45507.789335026729956)
    cc <- c(0.39894151208813466764, 8.8831497943883759412, 93.506656132177855979, 597.27027639480026226,
            2494.5375852903726711, 6848.1904505362823326, 11602.651437647350124, 9842.7148383839780218,
            1.0765576773720192317e-8)
    d <- c(22.266688044328115691, 235.38790178262499861, 1519.377599407554805, 6485.558298266760755,
           18615.571640885098091, 34900.952721145977266, 38912.003286093271411, 19685.429676859990727)
    p <- c(0.21589853405795699, 0.1274011611602473639, 0.022235277870649807, 0.001421619193227893466,
           2.9112874951168792e-5, 0.02307344176494017303)
    q <- c(1.28426009614491121, 0.468238212480865118, 0.0659881378689285515, 0.00378239633202758244,
           7.29751555083966205e-5)
    num1 <- sprintf("(x * %s)", horner("x*x", c(a[5L], a[1:4])))
    den1 <- horner("x*x", c(1, b[1:4]))
    num2 <- horner("abs(x)", c(cc[9L], cc[1:8]))
    den2 <- horner("abs(x)", c(1, d[1:8]))
    ix <- "(1.0/(x*x))"
    num3 <- sprintf("(%s * %s)", ix, horner(ix, c(p[6L], p[1:5])))
    den3 <- horner(ix, c(1, q[1:5]))
    xs <- "(trunc(abs(x)*16)/16)"
    tail <- function(temp) sprintf("(exp(-%s*%s*0.5) * exp(-((abs(x) - %s)*(abs(x) + %s))*0.5) * %s)", xs, xs, xs, xs, temp)
    t2 <- tail(sprintf("(%s / %s)", num2, den2))
    t3 <- tail(sprintf("((0.398942280401432677939946059934 - %s / %s) / abs(x))", num3, den3))
    pn <- sprintf(paste("CASE WHEN abs(x) <= 0.67448975 THEN 0.5 + %s / %s",
                        "WHEN abs(x) <= 5.656854249492380195206754896838 THEN (CASE WHEN x > 0 THEN 1 - %s ELSE %s END)",
                        "WHEN abs(x) < 50 THEN (CASE WHEN x > 0 THEN 1 - %s ELSE %s END)",
                        "ELSE (CASE WHEN x > 0 THEN 1.0 ELSE 0.0 END) END"), num1, den1, t2, t2, t3, t3)
    qa <- c(2509.0809287301226727, 33430.575583588128105, 67265.770927008700853, 45921.953931549871457,
            13731.693765509461125, 1971.5909503065514427, 133.14166789178437745, 3.387132872796366608)
    qb <- c(5226.495278852545925, 28729.085735721942674, 39307.89580009271061, 21213.794301586595867,
            5394.1960214247511077, 687.1870074920579083, 42.313330701600911252, 1)
    qc <- c(7.7454501427834140764e-4, .0227238449892691845833, .24178072517745061177, 1.27045825245236838258,
            3.64784832476320460504, 5.7694972214606914055, 4.6303378461565452959, 1.42343711074968357734)
    qd <- c(1.05075007164441684324e-9, 5.475938084995344946e-4, .0151986665636164571966, .14810397642748007459,
            .68976733498510000455, 1.6763848301838038494, 2.05319162663775882187, 1)
    qe <- c(2.01033439929228813265e-7, 2.71155556874348757815e-5, .0012426609473880784386, .026532189526576123093,
            .29656057182850489123, 1.7848265399172913358, 5.4637849111641143699, 6.6579046435011037772)
    qf <- c(2.04426310338993978564e-15, 1.4215117583164458887e-7, 1.8463183175100546818e-5, 7.868691311456132591e-4,
            .0148753612908506148525, .13692988092273580531, .59983220655588793769, 1)
    r1 <- "(0.180625 - (p - 0.5)*(p - 0.5))"
    rr <- "sqrt(-ln(least(p, 1 - p)))"
    r2 <- sprintf("(%s - 1.6)", rr)
    r3 <- sprintf("(%s - 5)", rr)
    qn <- sprintf(paste("CASE WHEN abs(p - 0.5) <= 0.425 THEN (p - 0.5) * %s / %s",
                        "ELSE (CASE WHEN p < 0.5 THEN -1 ELSE 1 END) * (CASE WHEN %s <= 5 THEN %s / %s ELSE %s / %s END) END"),
                  horner(r1, qa), horner(r1, qb), rr, horner(r2, qc), horner(r2, qd), horner(r3, qe), horner(r3, qf))
    ## dnorm as dnorm.c: for |x| >= 5, exp() of a split square
    x1 <- "(round(abs(x) * 65536) / 65536)"
    dn <- sprintf(paste("CASE WHEN abs(x) < 5 THEN 0.398942280401432677939946059934 * exp(-0.5 * x * x)",
                        "ELSE 0.398942280401432677939946059934 * exp(-0.5 * %s * %s) * exp((-0.5 * (abs(x) - %s) - %s) * (abs(x) - %s)) END"),
                  x1, x1, x1, x1, x1)
    pc <- "CASE WHEN abs(x) > 1 THEN (CASE WHEN x > 0 THEN 1 - atan(1 / x) / pi() ELSE -atan(1 / x) / pi() END) ELSE 0.5 + atan(x) / pi() END"
    ## -expm1(-u) for u >= 0, Kahan's trick: (1 - w) u / -ln(w), w = exp(-u)
    em <- "CASE WHEN u > 40 THEN 1.0 WHEN exp(-u) = 1 THEN u ELSE (1 - exp(-u)) * u / (-ln(exp(-u))) END"

    for (m in list(c("xtx_pnorm", "x", pn), c("xtx_qnorm", "p", qn), c("xtx_dnorm", "x", dn),
                   c("xtx_pcauchy", "x", pc), c("xtx_mexpm1m", "u", em)))
        xtx_execute(D, sprintf("CREATE OR REPLACE TEMP MACRO %s(%s) AS (%s)", m[1L], m[2L], m[3L]))
    D$src$holder$xtx_macros <- TRUE
    invisible()
}


## ---- coded groups, multi-group sums (the kernels of the fixed effects) ----

## In memory, codes are integers 1..G. A remote FE grouping carries the
## source key and its level table; public coded groupings use integer codes.
## Grouped sums with coded groups return all G groups in code order
## (zeros for empty groups), without labels, and use one-pass C kernels in
## memory (src/fe.c).
xtx_c_codes <- function(codes, G) {
    if (!xtx_is_column(codes)) codes <- as.integer(codes)
    structure(list(codes = codes, G = as.integer(G)), class = "xtx_codes")
}

## Fixed effects can group directly on source keys. Translate the returned
## group labels to solver indices only after DuckDB has aggregated the rows.
xtx_group_index <- function(g, values) {
    if (is.null(g$levels)) as.integer(values) else match(values, g$levels[[1L]])
}

xtx_is_codes <- function(x) inherits(x, "xtx_codes")

## sums of columns within coded groups: G x k
xtx_p_groupsum_codes <- function(col, group) {
    G <- group$G
    if (!xtx_is_column(col) && !xtx_is_column(group$codes)) {
        ## (the kernel takes a vector as one column: no copy of a double)
        W <- col
        if (!is.double(W)) storage.mode(W) <- "double"
        return(.Call(C_xtx_fe_gsum, group$codes, G, W, design_threads()))
    }
    res <- xtx_group_query(list(group$codes), list(col), nulls = TRUE)
    k <- xtx_ncols(col)
    out <- matrix(0, G, k)
    i <- xtx_group_index(group, res$g1)
    for (j in seq_len(k)) {
        v <- as.double(res[[paste0("v", j)]])
        v[as.numeric(res[[paste0("n", j)]]) > 0] <- NA_real_
        out[i, j] <- v
    }
    out
}

## Sums of a column (k columns) within each of several coded groupings, and
## of another column `total` (or NULL) over all rows: list(groups = one G x k
## matrix per grouping, as xtx_p_groupsum(); total = xtx_colsum(total)). In
## memory those calls; remotely one query (GROUP BY GROUPING SETS), whose
## rows the frame computes once (the frame of residuals, with its joins).
xtx_c_groupsums <- function(col, groups, total = NULL) {
    remote <- xtx_is_column(col) && all(vapply(groups, function(g) xtx_is_column(g$codes), NA))
    if (!remote || length(groups) + (!is.null(total)) < 2L)
        return(list(groups = lapply(groups, function(g) xtx_p_groupsum(col, g)),
                    total = if (!is.null(total)) xtx_colsum(total)))
    gs <- vapply(groups, function(g) unclass(g$codes)$sql, "")
    vs <- unclass(col)$sql
    ts <- if (!is.null(total)) unclass(total)$sql else character()
    fr <- xtx_frame_ready(do.call(xtx_column_frame, c(list(col, total), lapply(groups, `[[`, "codes"))), c(gs, vs, ts))
    k <- length(vs)
    sums <- function(e, p) c(sprintf("SUM(CAST(%s AS DOUBLE)) AS %s%d", e, p, seq_along(e)),
                             sprintf("SUM(CASE WHEN (%s) IS NULL THEN 1 ELSE 0 END) AS %sn%d", e, p, seq_along(e)))
    sets <- c(sprintf("(%s)", gs), if (length(ts)) "()")
    res <- xtx_query(fr, sprintf("SELECT %s FROM (%s) AS s WHERE %s GROUP BY GROUPING SETS (%s)",
                                paste(c(sprintf("%s AS g%d", gs, seq_along(gs)), sums(vs, "v"), sums(ts, "t")),
                                      collapse = ", "),
                                fr$src$select, fr$filter, paste(sets, collapse = ", ")))
    ## (codes are never missing: the rows of grouping j are those with g_j)
    gk <- lapply(seq_along(gs), function(j) !is.na(res[[paste0("g", j)]]))
    out <- lapply(seq_along(groups), function(j) {
        r <- gk[[j]]
        for (l in seq_along(gs)[-j]) r <- r & !gk[[l]]
        m <- matrix(0, groups[[j]]$G, k)
        i <- xtx_group_index(groups[[j]], res[[paste0("g", j)]][r])
        for (u in seq_len(k)) {
            v <- as.double(res[[paste0("v", u)]][r])
            v[as.numeric(res[[paste0("vn", u)]][r]) > 0] <- NA_real_
            m[i, u] <- v
        }
        m
    })
    tot <- NULL
    if (length(ts)) {
        r <- Reduce(`&`, lapply(gk, `!`))
        tot <- vapply(seq_along(ts), function(u) {
            if (!any(r)) return(0)
            if (as.numeric(res[[paste0("tn", u)]][r][1L]) > 0) return(NA_real_)
            v <- as.double(res[[paste0("t", u)]][r][1L])
            if (is.na(v)) 0 else v
        }, 0)
    }
    list(groups = out, total = tot)
}

## weighted counts by pair of coded groups: Ga x Gb dgCMatrix
xtx_p_crosstab_codes <- function(a, b, w = NULL) {
    if (!xtx_is_column(a$codes) && !xtx_is_column(b$codes) && !xtx_is_column(w)) {
        ## (the kernel returns the transposes of the crosstabs of its first
        ## argument with the others, and is fastest when the first is
        ## sorted)
        ct <- .Call(C_xtx_fe_crosstabs, b$codes, list(a$codes), if (is.null(w)) NULL else as.double(w), b$G,
                    a$G, design_threads())[[1L]]
        return(xtx_ct_matrix(ct, a$G, b$G))
    }
    res <- xtx_group_query(list(a$codes, b$codes), if (is.null(w)) list() else list(w))
    ## (the pairs, one per row, through the kernel of the rows in memory)
    ct <- .Call(C_xtx_fe_crosstabs, xtx_group_index(b, res$g2), list(xtx_group_index(a, res$g1)), as.double(res$v1), b$G,
                a$G, design_threads())[[1L]]
    xtx_ct_matrix(ct, a$G, b$G)
}

## A crosstab from the kernels (compressed columns list(p, i, x)) as an
## nr x nc dgCMatrix.
xtx_ct_matrix <- function(ct, nr, nc) {
    methods::new("dgCMatrix", i = ct[[2L]], p = ct[[1L]], x = ct[[3L]], Dim = c(as.integer(nr), as.integer(nc)))
}

## Column sums of a design within several coded groups at once:
## xtx_p_colsums_groups(x, groups, extra = NULL), x possibly scaled (s * x),
## groups a list of xtx_c_codes(), extra a column (k columns) or NULL. Returns
## one list(S = G x p sums of the scaled design, n = G sums of the scale (the
## counts without one), extra = G x k sums of scale * extra) per grouping.
## In memory: one pass of the structured kernel for all groupings when it
## applies (xtx_fe_struct_gsum), else the row kernels; remotely, one GROUP BY
## per grouping combines the design, scale, and extra column sums.
xtx_p_colsums_groups <- function(x, groups, extra = NULL) {
    if (inherits(x, "xtx_remote")) {
        s <- x$scale
        one <- xtx_column(xtx_remote_unscaled(x), if (is.null(s)) "CAST(1 AS DOUBLE)" else s)
        totals <- if (is.null(extra)) one else xtx_c_cbind(one, one * extra)
        x <- xtx_design_with(x, list(totals, groups))
        return(lapply(groups, function(g) {
            cs <- xtx_sql_colsums(x, unclass(g$codes)$sql, s, FALSE,
                                 group_totals = unclass(totals)$sql)
            i <- xtx_group_index(g, cs$keys)
            S <- matrix(0, g$G, ncol(x))
            S[i, ] <- cs$sums
            ns <- matrix(0, g$G, ncol(cs$totals))
            ns[i, ] <- cs$totals
            list(S = S, n = ns[, 1L], extra = if (!is.null(extra)) ns[, -1L, drop = FALSE])
        }))
    }
    ## in memory: the plain design and its row scale
    w <- NULL
    D <- x
    if (xtx_is_comb(x)) {
        if (!is.null(x$design2)) stop("A combination of two designs has no grouped sums here.", call. = FALSE)
        D <- x$design
        s <- x$scale
        w <- if (is.null(s)) NULL
             else if (is.double(s) && length(s) == nrow(D) && is.null(attributes(s))) s
             else rep_len(as.double(s), nrow(D))
    }
    xtx_gsums_memory(D, w, groups, extra)$gs
}

## In memory: the grouped sums of xtx_p_colsums_groups() of design D scaled by
## w (NULL for 1), and the crosstabs (as xtx_p_crosstab(), sums of w) of the
## pairs of groupings `pairs` (a list of c(a, b), indices in `groups`: the
## G_b x G_a matrix of t(crosstab(a, b))) whose dense table is small enough,
## all in one pass of the structured kernel when it applies; list(gs, ct),
## ct with NULL for the pairs left out.
xtx_gsums_memory <- function(D, w, groups, extra = NULL, pairs = list()) {
    G <- vapply(groups, function(g) g$G, 1L)
    codes <- lapply(groups, function(g) g$codes)
    ct <- vector("list", length(pairs))
    ## (the kernels take a vector as one column: no copy of a double)
    Y <- extra
    if (!is.null(Y) && !is.double(Y)) storage.mode(Y) <- "double"
    if (xtx_colsums_struct_ok(D, G)) {
        ## dense crosstabs: up to 2^22 cells, and at most two per row
        dense <- vapply(pairs, function(pr) as.double(G[pr[1L]]) * G[pr[2L]] <= min(2^22, 2 * nrow(D)), NA)
        pm <- if (any(dense)) matrix(as.integer(unlist(pairs[dense])) - 1L, nrow = 2L)
        M <- .Call(C_xtx_fe_struct_gsum, D$setup, D$vars, w, codes, as.integer(G), Y, TRUE, pm, design_threads())
        for (q in seq_along(which(dense))) {
            pr <- pairs[[which(dense)[q]]]
            ct[[which(dense)[q]]] <- xtx_ct_matrix(M[[length(groups) + q]], G[pr[2L]], G[pr[1L]])
        }
        info <- xtx_struct_info(D)
        icpt <- isTRUE(D$setup$intercept > 0L)
        p <- ncol(D)
        gs <- lapply(seq_along(groups), function(k) {
            m <- M[[k]]
            S <- matrix(0, G[k], p)
            q <- 1L
            for (t in seq_along(info)) {
                Tt <- info[[t]]
                if (t == 1L && icpt) {
                    S[, Tt$cols[1L, 1L]] <- m[1L, ]
                    next
                }
                L <- Tt$L
                ns <- nrow(Tt$cols)
                rows <- q + seq_len(L * ns)
                S[, Tt$cols] <- t(xtx_struct_contract(m[rows, , drop = FALSE], Tt, G[k]))
                q <- q + L * ns
            }
            list(S = S, n = m[1L, ], extra = if (is.null(Y)) NULL else t(m[q + seq_len(NCOL(Y)), , drop = FALSE]))
        })
        return(list(gs = gs, ct = ct))
    }
    one <- if (is.null(w)) rep(1, nrow(D)) else w
    nt <- design_threads()
    gs <- lapply(seq_along(groups), function(k) {
        J <- as_jac(D)
        W <- matrix(one, ncol = 1L)
        S <- if (jac_kernel_ok(J)) t(xtx_crossprod_kernel(J, W, codes[[k]], G[k]))
             else .Call(C_xtx_fe_gsum, codes[[k]], G[k], one * as.matrix(D), nt)
        ns <- .Call(C_xtx_fe_gsum, codes[[k]], G[k], cbind(one, if (!is.null(Y)) one * Y), nt)
        list(S = S, n = ns[, 1L], extra = if (!is.null(Y)) ns[, -1L, drop = FALSE])
    })
    list(gs = gs, ct = ct)
}

## The sums of a fit with fixed effects (absorb.R of xtx.fit): the grouped sums of
## xtx_p_colsums_groups(w * x, groups, extra) (x unscaled, w NULL for 1) and
## the crosstabs t(xtx_p_crosstab(groups[[b]], groups[[a]], w)) (G_a x G_b
## sums of w, transposed) for each pair c(a, b) of `pairs`: list(gs, ct).
## In memory one pass of the structured kernel for the grouped sums and
## the small crosstabs (dense tables), and one pass for the crosstabs of
## each grouping a with all its b (xtx_fe_crosstabs(): the rows by level of
## a); remotely the queries of xtx_p_colsums_groups() and xtx_p_crosstab().
xtx_p_fe_sums <- function(x, w, groups, extra = NULL, pairs = list()) {
    if (inherits(x, "xtx_remote") || xtx_is_comb(x)) {
        gs <- xtx_p_colsums_groups(if (is.null(w)) x else xtx_p_scale(x, w), groups, extra)
        ct <- lapply(pairs, function(pr) xtx_p_crosstab(groups[[pr[2L]]], groups[[pr[1L]]], w))
        return(list(gs = gs, ct = ct))
    }
    if (!is.null(w)) w <- as.double(w)
    r <- xtx_gsums_memory(x, w, groups, extra, pairs)
    ## the other crosstabs: one pass per grouping a, for all its b
    left <- which(vapply(r$ct, is.null, NA))
    for (a in unique(vapply(pairs[left], `[`, 1L, 1L))) {
        qs <- left[vapply(pairs[left], `[`, 1L, 1L) == a]
        bs <- vapply(pairs[qs], `[`, 1L, 2L)
        cts <- .Call(C_xtx_fe_crosstabs, groups[[a]]$codes, lapply(groups[bs], `[[`, "codes"), w,
                     groups[[a]]$G, vapply(groups[bs], function(g) g$G, 1L), design_threads())
        for (j in seq_along(qs)) r$ct[[qs[j]]] <- xtx_ct_matrix(cts[[j]], groups[[bs[j]]]$G, groups[[a]]$G)
    }
    r
}

xtx_colsums_struct_ok <- function(D, G) {
    if (!isTRUE(getOption("xtx.structured", TRUE))) return(FALSE)
    if (!xtx_struct_ok(D)) return(FALSE)
    info <- xtx_struct_info(D)
    Q <- 2 + sum(vapply(info, function(Tt) Tt$L * nrow(Tt$cols), 0))
    all(Q * as.double(G) <= 2^25)
}

## ---- more column functions (requested by the GLM work) ----

## one column of a multi-column column
xtx_c_col <- function(x, j) {
    if (!xtx_is_column(x)) return(if (is.matrix(x)) x[, j] else x)
    xtx_column_pick(x, j)
}

## a constant column over the rows of a column's frame (or a design)
xtx_c_const <- function(like, value) {
    if (inherits(like, "xtx_remote")) return(xtx_column(like, xtx_sql_double(value)))
    if (!xtx_is_column(like)) return(rep(as.double(value), NROW(like)))
    xtx_column(unclass(like)$frame, xtx_sql_double(value))
}

xtx_c_lgamma <- function(x) {
    if (!xtx_is_column(x)) return(lgamma(x))
    xtx_column(unclass(x)$frame, sprintf("lgamma(%s)", xtx_column_arith(x)))
}

## round half to even, as R's round(x)
xtx_c_round <- function(x) {
    if (!xtx_is_column(x)) return(round(x))
    xtx_column(unclass(x)$frame, sprintf("round_even(%s, 0)", xtx_column_arith(x)))
}

## is.finite, NaN-safe (in DuckDB, NaN = NaN and NaN > 0 are true)
xtx_c_is_finite <- function(x) {
    if (!xtx_is_column(x)) return(is.finite(x))
    xtx_column(unclass(x)$frame, sprintf("COALESCE(isfinite(CAST(%s AS DOUBLE)), FALSE)", unclass(x)$sql),
               notnull = TRUE, lgl = TRUE)
}

## TRUE where x is not missing, FALSE where it is NA or NaN (NULL remotely);
## infinite values count as present (= !is.na(x); was xtx_c_notnan, renamed
## in tapx: it never tested NaN alone.) Never NA, so that !xtx_c_present(x)
## counts the missing values. (In memory one pass in C for a plain vector.)
## Remotely, a column of the data as it is stored (s."v", a column of a
## key: design_key()) that is not double cannot be NaN, and may not be
## numeric at all: only IS NOT NULL (isnan() would not take text). A double
## column, or any other expression, is NaN-checked as well.
xtx_c_present <- function(x) {
    if (!xtx_is_column(x)) {
        if (is.null(attributes(x)) && (is.double(x) || is.integer(x) || is.logical(x)))
            return(.Call(C_xtx_notnan, x))
        return(!is.na(x))
    }
    sql <- unclass(x)$sql
    ## (by the type of the column: a logical column or a value stored as
    ## text, integer or logical is never NaN)
    plain <- xtx_column_lgl(x) | !xtx_column_type(x) %in% c("", "double")
    one <- ifelse(plain, sprintf("(%s IS NOT NULL)", sql), sprintf("COALESCE(NOT isnan(CAST(%s AS DOUBLE)), FALSE)", sql))
    xtx_column(unclass(x)$frame, one, notnull = TRUE, lgl = TRUE)
}

## The rows where a logical column is TRUE (with negate, where it is not
## TRUE: FALSE or missing), as a column to sum: 1 or 0. In memory a logical
## vector, which sums as 0/1.
xtx_c_count <- function(cond, negate = FALSE) {
    if (!xtx_is_column(cond)) {
        if (anyNA(cond)) cond <- !is.na(cond) & cond
        return(if (negate) !cond else cond)
    }
    sql <- unclass(cond)$sql
    xtx_column(unclass(cond)$frame, sprintf(if (negate) "CASE WHEN %s THEN 0 ELSE 1 END" else "CASE WHEN %s THEN 1 ELSE 0 END", sql),
             notnull = TRUE)
}

## A column computed once. Remotely, without `table`: a node whose columns
## are the expressions, so that later expressions name them instead of
## repeating them (DuckDB expands macros such as xtx_qnorm() inline, and
## nested expressions would grow without bound); a projection in the
## queries that read it, nothing in the database. Columns computed
## separately this way combine (their nodes are rendered on the same
## source). In memory, the column itself. With table = TRUE, a temporary
## table of the rows used (computed once for all later queries, which then
## scan it: for columns that many queries read, such as the working weights
## and responses of a GLM with fixed effects) with the columns of the
## source and the values (xtx_store()); its rows keep the filter of the
## design.
xtx_c_freeze <- function(col, table = FALSE) {
    if (!xtx_is_column(col)) return(col)
    u <- unclass(col)
    sql <- u$sql
    stem <- xtx_node_stem("f")
    nm <- paste0(stem, "_", seq_along(sql))
    if (table) {
        fr <- xtx_store(u$frame, stats::setNames(sql, nm))
    } else {
        fr <- xtx_column_frame1(col)
        fr <- xtx_frame_add(fr, xtx_node("proj", stem, nm, xtx_frame_nodes(fr), sprintf("%s AS %s", sql, xtx_qi(fr, nm))))
    }
    xtx_column(fr, sprintf("s.%s", xtx_qi(fr, nm)), notnull = u$notnull %||% FALSE, lgl = xtx_column_lgl(col),
               type = xtx_column_type(col))
}

## The values `exprs` (named SQL expressions of a column of frame `frame`)
## stored with the rows of the frame: a table (temporary, or the permanent
## table `name`, replaced if `overwrite`) of the rows used, with the columns
## of the source of the frame, of the nodes that the expressions name, and
## the values. Returns the frame on that table: a root that descends from
## the frame's and carries those nodes, so that the columns of the design
## evaluate on it. The expressions are checked (and the columns of the
## table read) by a query without rows before the table is written; a
## failed statement creates nothing and leaves a table of that name as it
## was.
xtx_store <- function(frame, exprs, name = NULL, overwrite = FALSE) {
    fr <- xtx_frame_ready(xtx_remote_unscaled(frame), exprs)
    root <- xtx_frame_root(fr)
    nodes <- xtx_frame_nodes(fr)
    tab <- xtx_remote_to_table(xtx_remote_extend(fr, exprs, "column"), name = name, overwrite = overwrite)
    xtx_frame_set(fr, xtx_root_child(tab$src, root, nodes))
}

## The frame of a remote design as a temporary table of its rows used (with
## every column of the frame; the design's filter stays valid on them).
## (`keep`: the columns of the frame to copy, NULL for all; `name`: a
## permanent table of that name instead, replaced if `overwrite`, for the
## public col_compute())
xtx_remote_to_table <- function(x, keep = NULL, name = NULL, overwrite = FALSE) {
    tn <- xtx_qi(x, if (is.null(name)) xtx_tmp_name("xtx_mat") else name)
    create <- if (is.null(name)) "CREATE TEMPORARY TABLE" else if (overwrite) "CREATE OR REPLACE TABLE" else "CREATE TABLE"
    sel <- if (is.null(keep)) "*" else paste(xtx_qi(x, keep), collapse = ", ")
    xtx_execute(x, sprintf("%s %s AS SELECT %s FROM (%s) AS s WHERE %s", create, tn, sel, x$src$select, x$filter))
    x$src$select <- sprintf("SELECT * FROM %s", tn)
    x$src$base_table <- tn
    if (!is.null(keep) && !is.null(x$src$proto)) {
        proto <- x$src$proto[intersect(keep, names(x$src$proto))]
        ## (columns of nodes, which the prototype of the source does not
        ## list: doubles)
        for (cn in setdiff(keep, names(proto))) proto[[cn]] <- double(0)
        x$src$proto <- proto[keep]
        x$src$types <- vapply(x$src$proto, xtx_col_type, "")
    }
    x
}

## The columns of the frame of remote design x to copy for the design and
## the columns `cols` (a list of columns, xtx_c_codes() groupings or NULL):
## the columns of the data (whose names do not start with "__xtx_"; results
## on the rows keep them), and the internal columns (stored values, and the
## columns of the nodes rendered in the source of x) that the design (the
## SQL and inputs of its variables, its scale) or `cols` read: those whose
## names occur in their SQL.
xtx_remote_used <- function(x, cols) {
    sql <- function(v) {
        if (is.null(v)) return(character())
        if (xtx_is_codes(v)) return(sql(v$codes))
        if (xtx_is_column(v)) return(unclass(v)$sql)
        if (is.list(v)) return(unlist(lapply(v, sql), use.names = FALSE))
        character()
    }
    txt <- c(sql(cols), unlist(lapply(x$vars, function(v) c(v$sql, v$col, v$inputs))), x$scale %||% character())
    txt <- paste(as.character(txt), collapse = " ")
    carried <- xtx_frame_root(x)$carried
    nodes <- xtx_frame_nodes(x)
    nm <- c(names(x$src$proto), unlist(lapply(nodes[!names(nodes) %in% carried], `[[`, "cols"), use.names = FALSE))
    nm <- unique(nm)
    nm[!startsWith(nm, "__xtx_") | vapply(nm, grepl, NA, x = txt, fixed = TRUE)]
}

## The design over a copy of its rows (see the header): a new design, whose
## source holds the rows used, the columns of the data and (given `cols`)
## the internal columns that the design and `cols` read. Columns of x carry
## over to it with xtx_c_on().
xtx_p_materialize <- function(x, cols = NULL) {
    if (!inherits(x, "xtx_remote")) return(x)
    ## (with the nodes of the columns: the copy holds them)
    if (!is.null(cols)) x <- xtx_design_with(x, cols)
    sc <- x$scale
    root <- xtx_frame_root(x)
    ids <- unique(c(root$carried, names(xtx_frame_nodes(x))))
    keep <- if (!is.null(cols) && !is.null(x$src$proto)) xtx_remote_used(x, cols)
    if (!length(keep)) keep <- NULL
    tab <- xtx_remote_to_table(xtx_remote_unscaled(x), keep)$src
    ## (the nodes whose columns were copied)
    if (!is.null(keep)) ids <- ids[vapply(ids, function(id) any(startsWith(keep, id)), NA)]
    tab$path <- NULL
    tab$carried <- ids
    tab$origin <- root$path %||% character()
    out <- xtx_frame_set(x, tab)
    out$filter <- "TRUE"
    out$cache <- new.env(parent = emptyenv())
    out$from <- c(x$id, x$from)
    out$id <- xtx_tmp_name("xtx")
    out$scale <- sc
    out
}


## ---- inputs of the methods of a design ----

## Which rows of a design matrix are finite: rowSums(!is.finite(as.matrix(x)))
## == 0, computed without forming the dense matrix (TRUE for each row whose
## entries are all finite, FALSE for a row with any NA, NaN or infinite entry,
## including entries that overflow, such as the product of two large numeric
## variables, and, for a combination of designs, scaled sums that overflow).
finite_rows <- function(x) {
    ## a remote design: a logical column (every entry of the row finite)
    if (inherits(x, "xtx_remote")) {
        P <- xtx_p_mult(x, diag(1, ncol(x)))
        return(xtx_column(unclass(P)$frame,
                        paste(sprintf("COALESCE(isfinite(%s), FALSE)", unclass(P)$sql), collapse = " AND ")))
    }
    checkmate::assert_class(x, "xtx_matrix")
    J <- as_jac(x)
    if (!jac_kernel_ok(J)) {
        return(unname(rowSums(!is.finite(as.matrix(x))) == 0))
    }
    do.call(.Call, c(list(C_xtx_jac_finite), jac_args(J), list(design_threads())))
}

## A lower bound on the absolute value of the non-zero entries of x (0 when
## no bound is available: always for a remote design).
xtx_p_min_nonzero <- function(x) {
    if (inherits(x, "xtx_remote")) return(0)
    terms <- xtx_terms(x)
    out <- Inf
    for (t in terms) {
        d <- t$design
        b <- if (kernel_ok(d)) .Call(C_xtx_min_nonzero, d$setup, d$vars) else {
            X <- as.matrix(d)
            a <- abs(X[X != 0 & is.finite(X)])
            if (length(a)) min(a) else 1
        }
        s <- t$scale
        if (!is.null(s)) {
            a <- abs(s[s != 0 & is.finite(s)])
            b <- b * (if (length(a)) min(a) else 1)
        }
        out <- min(out, b)
    }
    ## a sum of two designs can cancel to small values
    if (length(terms) > 1L) 0 else out
}

## The `y` of the qr() methods as a column of design x: in memory a finite
## double matrix (from a vector or matrix with nrow(x) rows); remotely a
## remote column from column names, a DBI::SQL() expression or a remote
## column. Returns list(Y, vec (was y one column?), rows, cols (names)).
xtx_p_y <- function(x, y) {
    if (!inherits(x, "xtx_remote")) return(xtx_qr_y(y, nrow(x)))
    D <- xtx_remote_unscaled(x)
    if (xtx_is_column(y)) return(list(Y = y, vec = xtx_ncols(y) == 1L, cols = NULL))
    if (inherits(y, "SQL")) return(list(Y = xtx_column(D, as.character(y)), vec = TRUE, cols = NULL))
    if (is.character(y)) {
        xtx_remote_column(D, y, "y", several = TRUE)
        return(list(Y = xtx_column(D, sprintf("CAST(s.%s AS DOUBLE)", xtx_qi(D, y))), vec = length(y) == 1L,
                    cols = if (length(y) > 1L) y))
    }
    stop("For a remote design, `y` must name columns of its data, or be a DBI::SQL() expression.", call. = FALSE)
}

## ---- rows of a design ----

## The design of the rows where `keep` (a logical column, NA counting as
## FALSE) is TRUE: x[keep, ] in memory; remotely, the same design with its
## filter restricted: a new design, another row domain even when every row
## is kept (columns of x are carried over to it explicitly, xtx_c_on()).
xtx_p_subset <- function(x, keep) {
    if (!inherits(x, "xtx_remote")) {
        if (xtx_is_column(keep)) stop("A remote column cannot select rows of a design in memory.", call. = FALSE)
        return(x[if (is.logical(keep) && !anyNA(keep)) keep else !is.na(keep) & keep, ])
    }
    if (inherits(keep, "SQL")) keep <- xtx_column(x, as.character(keep))
    if (!xtx_is_column(keep))
        stop("Rows of a remote design are selected by a logical remote column or a DBI::SQL() predicate.",
             call. = FALSE)
    out <- xtx_design_with(x, keep, short = FALSE)
    out$filter <- sprintf("(%s) AND COALESCE(%s, FALSE)", out$filter, unclass(keep)$sql)
    out$n <- as.numeric(xtx_query(out, sprintf("SELECT COUNT(*) AS n FROM (%s) AS s WHERE %s",
                                              out$src$select, out$filter))$n)
    out$cache <- new.env(parent = emptyenv())
    out$from <- c(x$id, x$from)
    out$id <- xtx_tmp_name("xtx")
    out
}

## A column of a design translated onto design x, a design made from it (a
## subset, a copy of its rows, the centered design: `x$from` names the
## designs x was made from): the same expressions and nodes over the rows
## of x. A column of any other design is refused. In memory, the column
## itself.
xtx_c_on <- function(col, x) {
    if (!xtx_is_column(col)) return(col)
    u <- unclass(col)
    fr <- xtx_column_frame1(col)
    if (!identical(fr$id, x$id) && !fr$id %in% x$from)
        stop("Columns of another remote design: only columns of the design that this one was made from ",
             "(by a subset, a copy of its rows or centering) carry over to it.", call. = FALSE)
    rx <- xtx_frame_root(x)
    rc <- xtx_frame_root(fr)
    if (!identical(fr$id, x$id) && !is.null(rx$origin)) {
        ## (x reads a copy of the rows: it holds what the relation it was
        ## copied from held)
        pc <- rc$path %||% character()
        if (length(pc) > length(rx$origin) || !identical(pc, rx$origin[seq_along(pc)]))
            stop("This column reads a table stored apart from the rows that the design copied.", call. = FALSE)
        root <- rx
    } else {
        root <- xtx_root_later(rx, rc)
    }
    nx <- xtx_frame_nodes(x)
    nc <- xtx_frame_nodes(fr)
    u$frame <- xtx_frame_set(xtx_remote_unscaled(x), root, c(nx, nc[!names(nc) %in% names(nx)]))
    u$neg <- NULL
    structure(u, class = "xtx_column")
}

## ---- arguments of a fit (round 7 of 2hba) ----

## Columns for the arguments of a fit over the rows of a remote design:
## exprs, a named list of unevaluated arguments (each a column name, a
## DBI::SQL() expression, or an R expression over the columns, evaluated in
## env); logical: the names of logical arguments (a subset). R expressions
## with a SQL translation run in the database; others are evaluated in R,
## batch by batch, into a temporary table (xtx_materialize()), which then is
## the source of the design. Returns list(x = the design (possibly on that
## new source), cols = named list of columns (NULL for NULL arguments),
## evaluated = names evaluated in R, with their reasons).
## In memory, the arguments are evaluated in `data` (a data frame, the rows
## of x) and returned as vectors.
xtx_p_resolve <- function(x, exprs, env, logical = character(), data = NULL) {
    if (!inherits(x, "xtx_remote")) {
        cols <- lapply(exprs, function(e) {
            if (is.null(e)) return(NULL)
            if (is.character(e) && length(e) == 1L) return(data[[e]])
            eval(e, data, env)
        })
        return(list(x = x, cols = cols, evaluated = list()))
    }
    con <- xtx_con(x)
    items <- list()
    later <- list()
    sqls <- list()
    for (nm in names(exprs)) {
        e <- exprs[[nm]]
        if (is.null(e)) next
        a <- xtx_resolve_or_eval(e, env, x$src$types, con, nm, bool = nm %in% logical)
        if (!is.null(a$expr)) {
            items[[nm]] <- a
            later[[paste0("__xtx_eval_", nm)]] <- a$expr
        } else {
            sqls[[nm]] <- a$sql
        }
    }
    if (length(later)) {
        ## (a table of the data with the values: a relation stored from the
        ## source, which holds every column of it)
        x$src <- xtx_root_child(xtx_materialize(x$src, later, env, logical = paste0("__xtx_eval_", logical),
                                                chunk_size = x$chunk_size), x$src)
        for (nm in names(items)) {
            col <- DBI::dbQuoteIdentifier(con, paste0("__xtx_eval_", nm))
            sqls[[nm]] <- if (nm %in% logical) sprintf("s.%s", col) else sprintf("CAST(s.%s AS DOUBLE)", col)
        }
    }
    cols <- lapply(names(exprs), function(nm) if (is.null(sqls[[nm]])) NULL else xtx_column(x, sqls[[nm]]))
    names(cols) <- names(exprs)
    list(x = x, cols = cols, evaluated = items)
}
