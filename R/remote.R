## Out-of-core designs: the rows of the design live in a DuckDB table, a SQL
## query or Parquet files, and never all at once in R.
##
## A remote design holds the source (connection and SQL), a design_plan
## (terms with fixed transformations, levels, contrasts), a description of
## each model variable (how to compute it, or its missingness, in SQL).
## Products run in the database: each term t has X_t[i, ] = N_t(i)
## C_t[F_t(i), ], with F_t(i) the combined level of its factors and N_t(i)
## the product of its numeric variables (see structured.R). Sums over rows
## of products of such columns are GROUP BY aggregates over factor codes;
## the small dense products with the C_t are done in R (remote_sql.R).
## A variable that SQL cannot compute is evaluated in R once, batch by
## batch, into columns of a temporary table (xtx_remote_evaluate()), and is
## then a column like any other.

#' Out-of-core designs on DuckDB
#'
#' With a DBI connection to DuckDB, or the path of Parquet files, as `data`,
#' [design_matrix()] returns a *remote design*: the rows stay in the
#' database, and products with the design run there, so that data larger
#' than memory can be used.
#'
#' ```
#' design_matrix(formula, con, table = "flights")
#' design_matrix(formula, con, sql = "SELECT * FROM flights WHERE year > 2000")
#' design_matrix(formula, "path/to/*.parquet")
#' design_matrix(spec, con, table = "new_flights")   # fit in R, score in the database
#' ```
#'
#' @section Sources:
#' `con` must be a [DBI][DBI::dbConnect()] connection to DuckDB (packages
#' DBI and duckdb). Give either `table`, the name of a table or view (a
#' character string or a [DBI::Id()]), or `sql`, a `SELECT` query. A
#' character vector of Parquet paths (globs allowed) opens a private
#' in-memory DuckDB connection that reads them with `read_parquet()`; it is
#' closed when the design is garbage collected.
#'
#' @section Training and new data:
#' With a formula, the design is learned from the remote data, as
#' `design_matrix(formula, data)` learns it from a data frame:
#' * Levels: those of an `ENUM` column (a factor in R) are its declared
#'   levels; those of a character column are its distinct values in the rows
#'   used, sorted as [factor()] sorts them in R (not by SQL collation).
#'   Logical columns are factors with levels `FALSE` and `TRUE`.
#' * Rows: rows with a missing value (SQL `NULL`, or `NaN`) in any model
#'   variable, response included, are dropped, as [stats::na.omit()] drops
#'   them.
#' * Transformations whose parameters depend on the data, such as
#'   [stats::poly()], [splines::ns()] or [scale()], are fitted as
#'   [stats::model.frame()] fits them, on the whole of their input columns:
#'   those columns (only those) are read into R once. To avoid this, fit
#'   the specification in R, on a sample or on other data, and pass it
#'   instead of the formula.
#' * Variables computed by an R expression that yields a factor, such as
#'   `factor(cyl)` or `interaction(a, b)`, take the levels that the
#'   expression gives on the distinct values of its input columns.
#'
#' With a specification (a design, a plan, a fitted model or another remote
#' design), the levels, contrasts and transformation parameters come from
#' it, and the remote data is only scored. Levels absent from the
#' specification raise the error that [stats::model.frame()] raises, and
#' columns of another type than in training (a numeric variable stored as
#' text, say) the error of [stats::predict.lm()]. As in memory, the rows
#' with missing values are kept by default (`na.action = na.pass`): their
#' predictions are `NULL`, and products involving their columns are `NA`.
#' With `na.action = na.omit`, they are dropped.
#'
#' @section Engine:
#' Products with a remote design are computed in DuckDB. SQL computes
#' factors, character and logical columns, numeric columns, interactions,
#' `factor(x)`, `as.factor(x)` and `as.character(x)` of a text, factor or
#' integer column (with the levels R gives them: in numeric order for
#' `factor()` of an integer column, sorted as text for `as.character()`),
#' and numeric expressions built from `+`, `-`, `*`, `/`, `^`, [I()], [log()]
#' (with base 2, 10 or a number), [log2()], [log10()], [exp()], [sqrt()] and
#' [abs()]. Any other variable, such as [stats::poly()], [splines::ns()],
#' `factor()` of a numeric expression or with `levels`, or a function of
#' yours, is evaluated in R once, in batches
#' of `chunk_size` rows, and written with the columns of the remote data to
#' a temporary table (one pass over the data, and a copy of it in the
#' database), whose columns the design then reads like any others; a
#' message names the first such term. The expression is evaluated on the
#' rows where its input columns are present, and a row where it gives `NA`
#' is a row with a missing value. With `engine = "sql"`, such a formula is
#' refused instead. `summary(D)` shows, for each term, whether SQL computes
#' it, and `summary(D, sql = TRUE)` the generated SQL.
#'
#' A term whose factors have more than about eight million combined
#' level-column cells is refused; for a factor with many levels, use fixed
#' effects (`y ~ x | f` in `xtx.fit::lm_xtx()` and `xtx.fit::glm_xtx()`).
#'
#' @section Operations:
#' The operations of a design in memory (see [design_algebra] and
#' [design_arithmetic]), computed in the database:
#' * `crossprod(D)`: `X'X`; `crossprod(D, "y")`: `X'y` for one or more
#'   columns of the remote data (or a `DBI::SQL()` expression).
#' * Scaled rows: `"w" * D` (a column name), `DBI::SQL("1 / v") * D` (an
#'   expression over the columns), `2 * D`, `D / "w"` and `-D` scale each
#'   row, as `w * D` with a numeric vector in memory. So `crossprod(D, "w" *
#'   D)` is `X'WX`, `crossprod("w" * D, "y")` is `X'Wy` and `colSums("w" * D,
#'   group = "g")` the weighted sums within groups. Only one design is
#'   involved: `crossprod(x, y)` needs `y` to be `x` with other scales.
#' * `colSums(D)`, `colMeans(D)`, and with `group = "g"` (a column, or a
#'   `DBI::SQL()` expression), the sums or means within groups (sorted as
#'   [base::rowsum()] sorts them).
#' * `D %*% beta`, `rowSums(D)`, `rowMeans(D)` and
#'   `quadratic_diagonal(D, V)`: a remote column ([design_column]), a
#'   lazy SQL expression over the rows of the design, where a design in
#'   memory returns a vector. Arithmetic on it builds another column;
#'   `as.vector(p)` reads its values into R, `as.data.frame(p, keep = "id")`
#'   reads them with columns of the remote data (to join them back), and
#'   [col_compute()] stores it in a table (temporary, or permanent with a
#'   `name`). `design_column(D, "y")`
#'   makes a column from a column of the remote data. Columns of one design
#'   combine, separate products included (`D %*% b1 + D %*% b2`): a
#'   product is SQL that each query reading it computes over the rows of
#'   the design, with the other columns of that query; reading a column
#'   creates nothing in the database.
#' * `as.matrix(D)` collects the rows and builds the dense matrix, only when
#'   there are at most `getOption("xtx.remote_max_rows", 1e5)` rows
#'   or with `force = TRUE`.
#' * `design_matrix(D, newdata)` applies the design to a data frame, and
#'   `design_matrix(D, con, table = ...)` to other remote data.
#' * `summary(D)`: the engine and, for each term, whether SQL computes it.
#'
#' @section Compression:
#' When the rows have few distinct *design patterns* (the same factor levels
#' and the same values of the numeric variables of the design, as with
#' factors and discrete numeric variables), `crossprod()`, `colSums()`,
#' `colMeans()`, `xtx.fit::lm_xtx()` and `xtx.fit::glm_xtx()` can compress: the one scan
#' of the rows also groups by the numeric variables, which returns each
#' pattern with its number of rows and the sums of the weights and
#' responses, and the result is computed from the design of the patterns,
#' in memory. The option `xtx.remote_compress` controls it.
#' `TRUE` always compresses (unless every row is distinct); `FALSE` never
#' does. With `"auto"` (default), the choice depends on the operation:
#' * `colSums()` and `colMeans()` are never compressed (grouping by the
#'   numeric variables costs more than it saves);
#' * `crossprod()` is compressed only for designs with many numeric columns
#'   (currently at least 8 distinct products of numeric variables; the
#'   threshold is provisional), where compression saves many per-row
#'   products, and when the patterns are at most 1% of the rows. Other
#'   designs skip the decision entirely;
#' * `xtx.fit::glm_xtx()` is compressed when the patterns are at most 20% of the
#'   rows (each IRLS iteration then reads the patterns instead of the rows);
#' * `xtx.fit::lm_xtx()` without weights is compressed when the patterns are at
#'   most 0.1% of the rows;
#' * weighted `xtx.fit::lm_xtx()` is compressed only for designs with many numeric
#'   columns (as `crossprod()`), when the patterns are at most 0.1% of the
#'   rows (the sum of squares within the patterns takes a second pass over
#'   the rows).
#'
#' These thresholds are provisional. The option
#' `xtx.remote_compress_ratio` (unset by default), when set, caps
#' every one of these ratios: under `"auto"`, the patterns must then also be
#' at most that fraction of the rows.
#'
#' In every case the number of patterns must also be at most
#' `getOption("xtx.remote_compress_max", 1e6)`. To decide, a design
#' with a continuous numeric column (a `DOUBLE` column with many distinct
#' values in a sample of the rows) is not compressed, at almost no cost;
#' otherwise the number of patterns is estimated over all the rows
#' (`approx_count_distinct()`). Designs that keep rows with missing values
#' (`na.action = na.pass`) are not compressed.
#'
#' Remote designs keep only the *Algebra* promise of the package: results
#' agree with the in-memory computation up to floating point rounding. With
#' infinite values in the data, the placement of `NaN` and `Inf` in sums can
#' differ from the dense computation, and `NaN` may appear where the dense
#' result has `NA`.
#'
#' @section Fits:
#' Linear and generalized linear models on remote data, with or without
#' fixed effects, are in the package 'xtx.fit': `xtx.fit::lm_xtx()` and
#' `xtx.fit::glm_xtx()`; its help page `remote_fits` describes how they run
#' on the remote products and how they differ from fits in memory.
#'
#' @section Arguments of remote data:
#' With remote data, [design_matrix()], [design_frame()], `xtx.fit::lm_xtx()`, `xtx.fit::glm_xtx()` and
#' `predict()` take, through `...`: `table` or `sql` (the remote rows: a
#' table name, or a `SELECT` query; not with Parquet paths), `engine`
#' (`"auto"`, or `"sql"` to refuse variables that must be evaluated in R)
#' and `chunk_size` (rows per batch when variables are evaluated in R,
#' default `getOption("xtx.chunk_size", 250000)`).
#' @aliases coerce,xtx_remote,CsparseMatrix-method coerce,xtx_remote,dgCMatrix-method
#' @examples
#' if (requireNamespace("duckdb", quietly = TRUE)) {
#'     con <- DBI::dbConnect(duckdb::duckdb())
#'     d <- data.frame(y = rnorm(100), x = rnorm(100), f = gl(4, 25))
#'     DBI::dbWriteTable(con, "d", d)
#'     D <- design_matrix(y ~ x * f, con, table = "d")
#'     D
#'     beta <- solve(crossprod(D), crossprod(D, "y"))
#'     p <- D %*% beta
#'     p
#'     head(as.vector(p))
#'     r <- design_column(D, "y") - p
#'     head(as.data.frame(r, nm = "resid"))
#'     summary(D)
#'     DBI::dbDisconnect(con)
#' }
#' @name design_remote
NULL

## Is `x` remote data for design_matrix()? (A character vector is never
## valid data for model.frame(), so it can only mean Parquet paths.)
xtx_is_remote_source <- function(x) {
    inherits(x, "DBIConnection") || (is.character(x) && length(x) > 0L && is.null(dim(x)))
}

xtx_need_duckdb <- function() {
    for (p in c("DBI", "duckdb")) {
        if (!requireNamespace(p, quietly = TRUE))
            stop(sprintf("Remote designs need the '%s' package.", p), call. = FALSE)
    }
}


## Entry point from design_matrix(), and from the fits (xtx_remote_frame()).
## For the fits: `subset` (SQL over the source) restricts the rows; `extra`
## (named list of list(sql, inputs, kind), kind "numeric" for a numeric
## argument such as weights, "present" for columns that must only be
## present) adds variables in no term, which select rows; factor levels
## absent from the rows used are dropped with `drop_unused`, as lm() builds
## its model frame. `source` can also be an opened source
## (xtx_source()); `override` gives SQL for model variables whose
## values were computed in R (xtx_materialize()).
xtx_remote <- function(object, source, contrasts.arg, xlev, dots, na.action = stats::na.pass,
                      drop_unused = FALSE, subset = NULL, extra = list(), override = list()) {
    xtx_need_duckdb()
    nms <- names(dots)
    if (is.null(nms)) nms <- rep("", length(dots))
    ra <- dots[nms %in% xtx_remote_keys]
    rest <- dots[!nms %in% xtx_remote_keys]
    engine <- ra$engine %||% "auto"
    checkmate::assert_choice(engine, c("auto", "sql"), .var.name = "engine")
    chunk_size <- ra$chunk_size %||% getOption("xtx.chunk_size", 250000)
    checkmate::assert_number(chunk_size, lower = 1, finite = TRUE, .var.name = "chunk_size")
    src <- if (inherits(source, "xtx_src")) source else xtx_source(source, ra$table, ra$sql)
    if (!is.null(subset)) {
        src$select <- sprintf("SELECT s.* FROM (%s) AS s WHERE %s", src$select, subset)
        src$base_table <- NULL
        src$label <- paste(src$label, "(subset)")
    }

    if (inherits(object, "xtx_formula")) {
        contrasts.arg <- object$contrasts
        xlev <- object$xlevels
        object <- object$formula
    }
    if (inherits(object, "formula")) {
        if (length(rest)) stop("Unknown arguments: ", paste(names(rest), collapse = ", "), call. = FALSE)
        D <- xtx_remote_formula(object, src, contrasts.arg, xlev, extra, drop_unused, override)
    } else {
        if (!is.null(contrasts.arg) || !is.null(xlev))
            stop("`contrasts.arg` and `xlev` apply only to a formula.", call. = FALSE)
        plan <- if (inherits(object, "xtx_plan")) {
            object
        } else if (inherits(object, "xtx_matrix")) {
            xtx_plan_for_design(object)
        } else {
            do.call(xtx_plan_of, c(list(object), rest))
        }
        D <- xtx_remote_spec(plan, src, xtx_remote_na(na.action))
    }
    D$chunk_size <- as.integer(chunk_size)
    xtx_remote_finish(D, engine)
}

## na.action of a specification applied to remote data: "pass" (rows with
## missing values kept, as by default in memory) or "omit".
xtx_remote_na <- function(na.action) {
    if (identical(na.action, stats::na.pass) || identical(na.action, "na.pass")) return("pass")
    if (identical(na.action, stats::na.omit) || identical(na.action, "na.omit") ||
        identical(na.action, stats::na.exclude) || identical(na.action, "na.exclude")) return("omit")
    stop("For remote data, `na.action` must be na.pass or na.omit.", call. = FALSE)
}

## The source: remote data before a design is built on it. The connection
## (in a holder environment, so that a private connection is closed once
## nothing refers to it), the SELECT statement that reads the rows, and the
## columns with their types. `source`: a DuckDB connection (with `table`, the
## name of a table, or `sql`, a query), or Parquet paths. Returns a list of
## class `xtx_src`: `holder` (holder$con, the connection), `select`,
## `base_table` (a table that can be read in batches by rowid), `label`,
## `proto` (a data frame with no rows and the columns of the data) and
## `types` (the kind of each column: "double", "integer", "logical",
## "character", "factor" or "other").
xtx_source <- function(source, table = NULL, sql = NULL) {
    holder <- new.env(parent = emptyenv())
    base_table <- NULL
    if (is.character(source)) {
        if (!is.null(table) || !is.null(sql))
            stop("Give `table` or `sql` only with a connection.", call. = FALSE)
        holder$con <- DBI::dbConnect(duckdb::duckdb())
        holder$own <- TRUE
        reg.finalizer(holder, function(e) {
            try(DBI::dbDisconnect(e$con), silent = TRUE)
        }, onexit = TRUE)
        con <- holder$con
        paths <- paste(DBI::dbQuoteString(con, source), collapse = ", ")
        select <- sprintf("SELECT * FROM read_parquet([%s])", paths)
        label <- sprintf("parquet %s", paste(source, collapse = ", "))
    } else {
        if (!inherits(source, "duckdb_connection"))
            stop("Remote designs support DuckDB connections only.", call. = FALSE)
        holder$con <- source
        holder$own <- FALSE
        con <- source
        if (is.null(table) == is.null(sql))
            stop("Give exactly one of `table` and `sql`.", call. = FALSE)
        if (!is.null(table)) {
            qt <- as.character(DBI::dbQuoteIdentifier(con, table))
            select <- sprintf("SELECT * FROM %s", qt)
            label <- sprintf("table %s", qt)
            ## a base table can be read in batches by rowid, without a copy
            ok <- tryCatch({
                DBI::dbGetQuery(con, sprintf("SELECT rowid FROM %s LIMIT 0", qt))
                TRUE
            }, error = function(e) FALSE)
            if (ok) base_table <- qt
        } else {
            checkmate::assert_string(sql, .var.name = "sql")
            select <- sub(";[[:space:]]*$", "", sql)
            label <- "SQL query"
        }
    }
    proto <- DBI::dbGetQuery(con, sprintf("SELECT * FROM (%s) AS s LIMIT 0", select))
    types <- vapply(proto, xtx_col_type, "")
    structure(list(holder = holder, select = select, base_table = base_table, label = label,
         proto = proto, types = types), class = "xtx_src")
}

xtx_con <- function(D) D$src$holder$con

xtx_col_type <- function(x) {
    if (is.factor(x)) return("factor")
    if (is.character(x)) return("character")
    if (is.logical(x)) return("logical")
    if (is.integer(x) && !is.object(x)) return("integer")
    if (is.double(x) && !is.object(x)) return("double")
    "other"
}

xtx_qi <- function(D, x) as.character(DBI::dbQuoteIdentifier(xtx_con(D), x))

## Run a query (with the option xtx.remote_debug, print it).
xtx_query <- function(D, sql) {
    if (isTRUE(getOption("xtx.remote_debug"))) message(sql)
    DBI::dbGetQuery(xtx_con(D), sql)
}

xtx_execute <- function(D, sql) {
    if (isTRUE(getOption("xtx.remote_debug"))) message(sql)
    DBI::dbExecute(xtx_con(D), sql)
}

## Names of temporary tables, views and columns: a counter (8 hex digits,
## right after the prefix) and a token of the process, without using (or
## being reset by) R's random numbers.
xtx_tmp <- new.env(parent = emptyenv())
xtx_tmp_name <- function(prefix) {
    if (is.null(xtx_tmp$token))
        xtx_tmp$token <- sprintf("%x%x", Sys.getpid(), as.integer((as.numeric(Sys.time()) * 1000) %% 1e9))
    xtx_tmp$i <- (xtx_tmp$i %||% 0) + 1
    sprintf("%s_%08x%s", prefix, as.integer(xtx_tmp$i), xtx_tmp$token)
}

## Describe one model variable: how SQL can compute it, or test whether it
## is missing.
##   kind: "factor" (a character or ENUM column), "logical" (a logical
##         column), "numeric" (an expression with a SQL translation) or
##         "other" (evaluated in R).
## factor(x), as.factor(x) or as.character(x) of a column x that is text,
## an ENUM or an integer: a factor variable on that column, which the
## database groups by its text (xtx_sql_key()). Its levels are those R
## gives: of factor(x), the values in their order (of an integer column,
## numeric; of an ENUM, that of its type); of as.character(x), a character
## variable, which a model frame makes a factor with levels sorted as text
## (`text`). Other columns (a double has no text that is R's as.character()
## for every value; a logical has other labels) and other calls (with
## `levels` or `labels`) are evaluated in R. list(col, text), or NULL.
xtx_factor_call <- function(e, types) {
    if (!is.call(e) || length(e) != 2L || !is.name(e[[1L]]) || !is.name(e[[2L]])) return(NULL)
    if (!is.null(names(e)) && any(nzchar(names(e)))) return(NULL)
    fn <- as.character(e[[1L]])
    col <- as.character(e[[2L]])
    if (!fn %in% c("factor", "as.factor", "as.character") || !isTRUE(types[col] %in% c("character", "factor", "integer")))
        return(NULL)
    list(col = col, text = fn == "as.character")
}

xtx_remote_var <- function(name, expr, types, D, xlevels, formula_case) {
    v <- list(name = name, expr = expr, inputs = intersect(all.vars(expr), names(types)),
              kind = "other", reason = NULL, sql = NULL, col = NULL, levels = NULL)
    in_xlev <- name %in% names(xlevels)
    if (is.name(expr) && as.character(expr) %in% names(types)) {
        col <- as.character(expr)
        ty <- types[[col]]
        if (ty %in% c("character", "factor")) {
            if (formula_case || in_xlev) {
                v$kind <- "factor"
                v$col <- col
            } else {
                v$reason <- sprintf("`%s` is a character column but not a factor in the specification", col)
            }
            return(v)
        }
        if (ty == "logical") {
            if (formula_case || !in_xlev) {
                v$kind <- "logical"
                v$col <- col
            } else {
                v$reason <- sprintf("`%s` is a logical column but a factor in the specification", col)
            }
            return(v)
        }
    }
    ## factor(x), as.factor(x), as.character(x) of a text, ENUM or integer
    ## column: the factor of the column
    fc <- xtx_factor_call(expr, types)
    if (!is.null(fc) && (formula_case || in_xlev)) {
        v$kind <- "factor"
        v$col <- fc$col
        v$text <- fc$text
        return(v)
    }
    if (in_xlev) {
        v$reason <- sprintf("`%s` is a factor in the specification", name)
        return(v)
    }
    ## offset(e): a numeric variable in no term
    e <- if (is.call(expr) && identical(expr[[1L]], as.name("offset")) && length(expr) == 2L) expr[[2L]] else expr
    tr <- xtx_sql_num(e, types, xtx_con(D))
    if (is.null(tr$sql)) {
        v$reason <- tr$reason
    } else {
        v$kind <- "numeric"
        v$sql <- tr$sql
    }
    v
}

## SQL condition: is variable `v` present (not NULL, not NaN) in a row?
## Exact for variables with a SQL translation; a variable evaluated in R
## counts as present when all its input columns are, as do the columns of
## a "present" variable (FE variables, a factor response).
xtx_remote_var_ok <- function(D, v) {
    types <- D$src$types
    cond <- switch(v$kind,
        factor = , logical = sprintf("s.%s IS NOT NULL", xtx_qi(D, v$col)),
        numeric = sprintf("(%s) IS NOT NULL AND NOT isnan(%s)", v$sql, v$sql),
        other = , present = vapply(v$inputs, function(col) {
            q <- xtx_qi(D, col)
            if (types[[col]] == "double") sprintf("s.%s IS NOT NULL AND NOT isnan(s.%s)", q, q)
            else sprintf("s.%s IS NOT NULL", q)
        }, ""))
    cond <- unique(cond)
    if (!length(cond)) "TRUE" else paste(sprintf("(%s)", cond), collapse = " AND ")
}

## WHERE clause: the rows used. With na.omit, the rows with no missing model
## variable; with na.pass (a specification scoring remote data), all rows.
xtx_remote_filter <- function(D) {
    if (identical(D$na, "pass")) return("TRUE")
    cond <- unique(vapply(D$vars, function(v) xtx_remote_var_ok(D, v), ""))
    cond <- setdiff(cond, "TRUE")
    if (!length(cond)) "TRUE" else paste(sprintf("(%s)", cond), collapse = " AND ")
}

## Variables of a terms object, with their predvars.
xtx_terms_vars <- function(tt) {
    vars <- as.list(attr(tt, "variables"))[-1L]
    pv <- attr(tt, "predvars")
    exprs <- if (is.null(pv)) vars else as.list(pv)[-1L]
    list(names = vapply(vars, deparse2, ""), exprs = exprs, calls = vars)
}

## A remote design from a formula: learn levels and transformations from the
## remote data.
xtx_remote_formula <- function(formula, src, contrasts.arg, xlev, extra = list(), drop_unused = FALSE,
                              override = list()) {
    D <- list(src = src, na = "omit")
    tt <- stats::terms(formula, data = src$proto)
    tv <- xtx_terms_vars(tt)
    types <- src$types
    D$vars <- Map(function(nm, e) xtx_remote_var(nm, e, types, D, NULL, TRUE), tv$names, tv$exprs)
    names(D$vars) <- tv$names
    ## variables computed in R into columns of the source
    for (nm in names(override)) {
        D$vars[[nm]]$kind <- "numeric"
        D$vars[[nm]]$sql <- override[[nm]]$sql
        D$vars[[nm]]$reason <- NULL
        D$vars[[nm]]$inputs <- unique(c(D$vars[[nm]]$inputs, override[[nm]]$inputs))
    }
    ## arguments of a fit (weights, offset, response): variables in no term,
    ## which select rows
    for (nm in names(extra)) {
        D$vars[[nm]] <- list(name = nm, expr = NULL, inputs = extra[[nm]]$inputs,
                             kind = extra[[nm]]$kind %||% "numeric",
                             reason = NULL, col = NULL, levels = NULL, sql = extra[[nm]]$sql)
    }
    D$filter <- xtx_remote_filter(D)

    ## a sample of the rows used, to see what the R expressions return
    sample <- xtx_remote_sample(D)
    env <- environment(tt)
    predvars <- tv$exprs
    xlevels <- list()
    for (j in seq_along(D$vars)) {
        v <- D$vars[[j]]
        if (v$kind == "factor") {
            ## (as.character(x): a character variable, whose levels are the
            ## values that occur, sorted as text)
            x <- if (types[[v$col]] == "factor" && !drop_unused && !isTRUE(v$text)) src$proto[[v$col]]
                 else xtx_remote_distinct(D, v$col)[[1L]]
            if (isTRUE(v$text)) x <- as.character(x)
            x <- as.factor(x)
            if (drop_unused) x <- droplevels(x)
            xlevels[[v$name]] <- levels(x)
        } else if (v$kind == "other") {
            val <- tryCatch(eval(v$expr, sample, env), error = function(e) e)
            if (is.factor(val) || is.character(val)) {
                ## levels from the distinct values of the inputs
                inp <- xtx_remote_distinct(D, v$inputs)
                fx <- as.factor(eval(v$expr, inp, env))
                lv <- levels(if (drop_unused) droplevels(fx) else fx)
                xlevels[[v$name]] <- lv
            } else if (!is.logical(val)) {
                ## does model.frame() fix parameters from the data?
                changed <- inherits(val, "error") ||
                    !identical(stats::makepredictcall(val, tv$calls[[j]]), tv$calls[[j]])
                if (changed) {
                    full <- xtx_query(D, sprintf("SELECT %s FROM (%s) AS s",
                                                paste0("s.", xtx_qi(D, v$inputs), collapse = ", "),
                                                src$select))
                    predvars[[j]] <- stats::makepredictcall(eval(v$expr, full, env), tv$calls[[j]])
                }
            }
        }
    }
    for (nm in names(xlev)) xlevels[[nm]] <- xlev[[nm]]
    attr(tt, "predvars") <- as.call(c(list(as.name("list")), predvars))
    for (j in seq_along(predvars)) D$vars[[j]]$expr <- predvars[[j]]
    if (!nrow(sample))
        stop(structure(class = c("xtx_no_rows", "error", "condition"),
                       list(message = "The remote data has no complete rows.", call = NULL)))
    D$plan <- xtx_plan_formula(tt, sample, contrasts.arg = contrasts.arg,
                                  xlev = if (length(xlevels)) xlevels else NULL)
    D$terms_full <- tt
    resp <- attr(tt, "response")
    D$response <- if (resp > 0L) tv$names[[resp]] else NULL
    D
}

## A remote design from a specification: score the remote data.
xtx_remote_spec <- function(plan, src, na = "pass") {
    D <- list(src = src, plan = plan, terms_full = plan$terms, response = NULL, na = na)
    tv <- xtx_terms_vars(plan$terms)
    xl <- plan$xlevels
    xtx_remote_check_classes(plan, src, tv, xl)
    D$vars <- Map(function(nm, e) xtx_remote_var(nm, e, src$types, D, xl, FALSE), tv$names, tv$exprs)
    names(D$vars) <- tv$names
    D$filter <- xtx_remote_filter(D)
    D
}

## Columns of another type than in training (a numeric variable stored as
## text, say): the error of predict.lm(), before anything is built. Only
## plain columns are checked here; variables computed by R expressions are
## checked when the prototype rows are built (xtx_chunk_design()).
xtx_remote_check_classes <- function(plan, src, tv, xl) {
    cl <- plan$data_classes
    if (is.null(cl)) return(invisible())
    m <- list()
    for (j in seq_along(tv$exprs)) {
        e <- tv$exprs[[j]]
        nm <- tv$names[[j]]
        if (!is.name(e) || !as.character(e) %in% names(src$proto) || !nm %in% names(cl)) next
        x <- src$proto[[as.character(e)]]
        ## model.frame() turns a character variable with levels into a factor
        if (is.character(x) && nm %in% names(xl)) x <- factor(x)
        m[[nm]] <- x
    }
    if (length(m)) stats::.checkMFClasses(cl, m)
    invisible()
}

## Model variables that SQL cannot compute (kind "other": in a term, or the
## response) are evaluated in R, batch by batch, into columns of a temporary
## table, which becomes the source of the design (xtx_materialize()): a
## numeric variable is then a column (several for a matrix, as poly(x, 2)
## gives), a factor, character or logical one a column of that type. With
## engine = "sql" this is refused; otherwise a message names the first such
## term. (offset() terms in no term stay as they are: a fit or a prediction
## evaluates them.)
xtx_remote_evaluate <- function(D, engine) {
    tt <- D$plan$terms
    tv <- xtx_terms_vars(tt)$names
    fac <- attr(tt, "factors")
    labels <- attr(tt, "term.labels")
    other <- names(D$vars)[vapply(D$vars, function(v) v$kind == "other", NA)]
    in_term <- function(t) intersect(tv[fac[, t] > 0], other)
    used <- unique(unlist(lapply(seq_along(labels), in_term)))
    todo <- intersect(other, c(used, D$response, setdiff(other, tv[attr(tt, "offset")])))
    if (!length(todo)) return(D)
    first <- todo[1L]
    what <- "variable"
    if (length(used)) {
        t <- which(vapply(seq_along(labels), function(t) length(in_term(t)) > 0L, NA))[1L]
        first <- in_term(t)[1L]
        what <- "term"
        label <- labels[t]
    } else if (isTRUE(D$response %in% todo)) {
        first <- label <- D$response
        what <- "response"
    } else {
        label <- first
    }
    why <- sprintf("%s `%s` cannot be translated to SQL (%s)", what, label, D$vars[[first]]$reason)
    if (engine == "sql") stop("The SQL engine cannot compute this design: ", why, ".", call. = FALSE)
    message("Evaluated in R: ", why, ".")
    idx <- match(todo, names(D$vars))
    cols <- paste0("__xtx_var", idx)
    xl <- names(D$plan$xlevels)
    src <- xtx_materialize(D$src, stats::setNames(lapply(D$vars[idx], function(v) v$expr), cols),
                          environment(D$terms_full) %||% globalenv(), chunk_size = D$chunk_size,
                          model = stats::setNames(lapply(D$vars[idx], function(v) {
                              list(inputs = v$inputs, factor = v$name %in% xl, label = v$name,
                                   class = D$plan$data_classes[v$name][!is.na(names(D$plan$data_classes[v$name]))])
                          }), cols))
    D$src <- src
    info <- attr(src, "model")
    for (i in seq_along(idx)) {
        v <- D$vars[[idx[i]]]
        f <- info[[cols[i]]]
        v$evaluated <- v$reason
        v$reason <- NULL
        v$kind <- f$kind
        if (f$kind == "numeric") {
            v$sql <- sprintf("CAST(s.%s AS DOUBLE)", xtx_qi(D, f$cols))
            v[["k"]] <- length(f$cols)
        } else {
            v$col <- f$cols
        }
        D$vars[[idx[i]]] <- v
    }
    D$filter <- xtx_remote_filter(D)
    D$evaluated <- sprintf("%s `%s` is evaluated in R", what, label)
    D
}

## Levels, columns, term structures.
xtx_remote_finish <- function(D, engine) {
    D <- xtx_remote_evaluate(D, engine)
    plan <- D$plan
    xl <- plan$xlevels
    ## levels of the factors, and new levels (as model.frame() reports them)
    for (j in seq_along(D$vars)) {
        v <- D$vars[[j]]
        if (v$kind == "logical") D$vars[[j]]$levels <- c("FALSE", "TRUE")
        if (v$kind != "factor") next
        lv <- xl[[v$name]]
        D$vars[[j]]$levels <- lv
        x <- xtx_remote_distinct(D, v$col)[[1L]]
        if (isTRUE(v$text)) x <- as.character(x)
        xi <- if (is.factor(x)) x[, drop = TRUE] else as.factor(x)
        nxl <- levels(xi)
        m <- is.na(match(nxl, lv))
        if (any(m))
            stop(sprintf(ngettext(length(m), "factor %s has new level %s", "factor %s has new levels %s"),
                         v$name, paste(nxl[m], collapse = ", ")), call. = FALSE)
    }
    sample <- xtx_remote_sample(D)
    proto <- xtx_chunk_design(D, sample)$design
    D$assign <- if (is.null(proto$setup)) attr(proto$dense, "assign") else proto$setup$assign
    D$colnames <- colnames(proto)
    D$p <- length(D$colnames)
    D$n <- as.numeric(xtx_query(D, sprintf("SELECT COUNT(*) AS n FROM (%s) AS s WHERE %s",
                                           D$src$select, D$filter))$n)
    D$cache <- new.env(parent = emptyenv())
    ## identifies the design through copies with scaled rows
    D$id <- xtx_tmp_name("xtx")
    D$terms <- xtx_remote_terms(D)
    D$engine <- "sql"
    D <- xtx_sql_prepare(D, proto$setup)
    class(D) <- "xtx_remote"
    D
}

## Distinct values of columns, in the rows used.
xtx_remote_distinct <- function(D, cols) {
    xtx_query(D, sprintf("SELECT DISTINCT %s FROM (%s) AS s WHERE %s",
                        paste0("s.", xtx_qi(D, cols), collapse = ", "), D$src$select, D$filter))
}

## The columns that the model variables read.
xtx_remote_inputs <- function(D) {
    unique(unlist(lapply(D$vars, function(v) v$inputs)))
}

xtx_remote_sample <- function(D, n = 2000L) {
    cols <- xtx_remote_inputs(D)
    sel <- if (length(cols)) paste0("s.", xtx_qi(D, cols), collapse = ", ") else "1 AS xtx_one"
    xtx_query(D, sprintf("SELECT %s FROM (%s) AS s WHERE %s LIMIT %d", sel, D$src$select, D$filter, n))
}

## Terms of the design: their variables, and whether SQL computes them
## (`sql_ok`; otherwise `reason` says why a variable was evaluated in R).
xtx_remote_terms <- function(D) {
    tt <- D$plan$terms
    fac <- attr(tt, "factors")
    labels <- attr(tt, "term.labels")
    ## rows of `factors` are the variables, in order
    pv <- xtx_terms_vars(tt)$names
    out <- list()
    if (attr(tt, "intercept") == 1L && !isTRUE(D$plan$drop_intercept))
        out[[1L]] <- list(index = 0L, label = "(Intercept)", vars = integer(), sql_ok = TRUE, reason = NULL)
    for (t in seq_along(labels)) {
        vnames <- pv[fac[, t] > 0]
        L <- prod(vapply(D$vars[vnames], function(v) max(1L, length(v$levels)), 1L))
        if (L * sum(D$assign == t) > xtx_sql_max_cells)
            stop(sprintf("Term `%s` has %s combined levels: too many for a remote design (for a factor with many levels, use fixed effects: y ~ x | f).",
                         labels[t], format(L, big.mark = ",")), call. = FALSE)
        reason <- unlist(lapply(D$vars[vnames], function(v) v$evaluated))[1L]
        out[[length(out) + 1L]] <- list(index = t, label = labels[t], vars = match(vnames, names(D$vars)),
                                        sql_ok = is.null(reason), reason = if (!is.null(reason)) unname(reason))
    }
    out
}

#' @export
print.xtx_remote <- function(x, ...) {
    cat(sprintf("<remote design_matrix: %s rows x %d columns>\n", xtx_fmt_count(x$cache$n %||% x$n), x$p))
    cat("source: ", x$src$label, "\n", sep = "")
    if (!is.null(x$scale)) cat("rows scaled by: ", x$scale, "\n", sep = "")
    if (is.null(x$evaluated)) {
        cat("engine: sql (all terms pushed down)\n")
    } else {
        cat(sprintf("engine: sql, on columns computed in R (%s)\n", x$evaluated))
    }
    xtx_print_columns(x$colnames)
    invisible(x)
}

xtx_fmt_count <- function(n) format(n, big.mark = ",", scientific = FALSE, trim = TRUE)

#' @export
dim.xtx_remote <- function(x) c(x$cache$n %||% x$n, x$p)

#' @export
dimnames.xtx_remote <- function(x) list(NULL, x$colnames)


## Columns of the remote data that must be numeric (double, integer or
## logical).
xtx_numeric_cols <- function(D, cols) {
    for (col in cols) {
        if (!D$src$types[[col]] %in% c("double", "integer", "logical"))
            stop(sprintf("Column `%s` is not numeric.", col), call. = FALSE)
    }
}

## A column name argument (y, weights, group): check and quote.
xtx_remote_column <- function(D, x, what, several = FALSE) {
    if (several) checkmate::assert_character(x, min.len = 1L, any.missing = FALSE, .var.name = what)
    else checkmate::assert_string(x, .var.name = what)
    miss <- setdiff(x, names(D$src$types))
    if (length(miss))
        stop(sprintf("`%s`: no column %s in the remote data.", what, paste0("`", miss, "`", collapse = ", ")),
             call. = FALSE)
    x
}

## Columns of a remote design as SQL expressions over its source `s`, the
## one form in which the products take them: from the names of columns of
## the remote data, a DBI::SQL() expression over them (checked here), or
## remote columns of the fitting core (primitives.R; the design then has
## their nodes in its source). Returns list(D = the design to use, sql, names).
## `numeric`: named columns must be numeric.
xtx_remote_cols <- function(D, x, what, several = FALSE, numeric = TRUE) {
    if (xtx_is_column(x)) {
        return(list(D = xtx_design_with(D, x, short = FALSE), sql = unclass(x)$sql, names = NULL))
    }
    if (inherits(x, "SQL")) {
        checkmate::assert_string(as.character(x), .var.name = what)
        e <- sprintf("(%s)", as.character(x))
        tryCatch(xtx_query(D, sprintf("SELECT %s FROM (%s) AS s LIMIT 0", e, D$src$select)),
                 error = function(err) stop(sprintf("`%s`: invalid SQL expression: %s", what, conditionMessage(err)),
                                            call. = FALSE))
        return(list(D = D, sql = e, names = as.character(x)))
    }
    xtx_remote_column(D, x, what, several)
    if (numeric) xtx_numeric_cols(D, x)
    list(D = D, sql = sprintf("s.%s", xtx_qi(D, x)), names = x)
}

## The design over the same rows, with extra columns `exprs` (named SQL
## expressions over the columns of the source) added to the source.
xtx_remote_extend <- function(D, exprs, what = "expression") {
    con <- xtx_con(D)
    items <- sprintf("(%s) AS %s", exprs, xtx_qi(D, names(exprs)))
    select <- sprintf("SELECT s.*, %s FROM (%s) AS s", paste(items, collapse = ", "), D$src$select)
    proto <- tryCatch(DBI::dbGetQuery(con, sprintf("SELECT * FROM (%s) AS s LIMIT 0", select)),
                      error = function(e) stop(sprintf("`%s`: invalid SQL expression: %s", what,
                                                       conditionMessage(e)), call. = FALSE))
    ## (the nodes of the design are in the select extended)
    D$src$select <- select
    D$base <- NULL
    D$nodes <- NULL
    D$src$proto <- proto
    D$src$types <- vapply(proto, xtx_col_type, "")
    ## the extra columns are not in the base table: batches read a copy
    D$src$base_table <- NULL
    D
}

## X'WX, X'WY, Y'WY, the number of rows and the sum of the weights of an
## unscaled remote design: from its patterns when that pays (see
## *Compression*; not with `wy`, or with compress = FALSE), else one scan.
## `ys`, `w`, `wy`: SQL expressions over the source (xtx_remote_cols(),
## xtx_sql_moments()).
xtx_remote_moments <- function(x, ys = character(), w = NULL, gram = TRUE, wy = NULL, compress = TRUE) {
    cap <- if (compress && is.null(wy)) xtx_compress_cap(x, op = "products")
    out <- if (!is.null(cap)) xtx_sql_moments(x, ys, w, gram, cap = cap)
    if (is.null(out)) out <- xtx_sql_moments(x, ys, w, gram, wy = wy)
    out
}

## The product of row scales (SQL expressions, or NULL), or NULL.
xtx_sql_times <- function(...) {
    w <- c(...)
    if (length(w)) paste(sprintf("(%s)", w), collapse = " * ")
}

## X'X of a remote design (rows possibly scaled) and, with the remote
## columns Y, X'Y, from one scan: list(xtx, xty). (Without Y,
## crossprod(x), and xty NULL.)
xtx_remote_gram_xty <- function(x, Y = NULL) {
    if (is.null(Y) || (!xtx_is_column(Y) && !length(Y))) return(list(xtx = crossprod(x), xty = NULL))
    yc <- xtx_remote_cols(x, Y, "y")
    sc <- yc$D$scale
    if (is.null(sc) && length(yc$sql) == 1L)
        return(design_moments(yc$D, Y, NULL, centered = FALSE))
    xtx_remote_moments(xtx_remote_unscaled(yc$D), yc$sql, xtx_sql_times(sc, sc), wy = sc, compress = FALSE)[c("xtx", "xty")]
}

#' @export
print.xtx_remote_table <- function(x, ...) {
    cat(sprintf("<remote table \"%s\">: read it with as.data.frame(), or query it on attr(, \"connection\")\n",
                unclass(x)[1L]))
    invisible(x)
}

#' @export
as.data.frame.xtx_remote_table <- function(x, ...) {
    con <- attr(x, "connection")
    DBI::dbGetQuery(con, sprintf("SELECT * FROM %s", DBI::dbQuoteIdentifier(con, unclass(x)[1L])))
}

#' @export
`%*%.xtx_remote` <- function(x, y) {
    if (!inherits(x, "xtx_remote"))
        stop("Only `design %*% beta` is supported, with the design first.", call. = FALSE)
    ## a remote column (design_column), as a vector in memory
    k <- if (is.matrix(y)) ncol(y) else 1L
    lab <- if (k == 1L) ".fitted" else colnames(y)
    if (k > 1L && (is.null(lab) || anyNA(lab) || any(lab == "") || anyDuplicated(lab)))
        lab <- paste0(".fitted", seq_len(k))
    xtx_column_label(xtx_p_mult(x, y), lab)
}

## Scaled rows, as `w * D` for a design in memory: the scale is a SQL
## expression over the columns of the remote data (a number, a column
## name, or a DBI::SQL() expression), kept unevaluated in `x$scale`.
xtx_remote_scale_sql <- function(D, s, op) {
    if (inherits(s, "SQL")) {
        checkmate::assert_string(as.character(s), .var.name = "scale")
        e <- sprintf("(%s)", as.character(s))
    } else if (is.character(s)) {
        xtx_numeric_cols(D, xtx_remote_column(D, s, "scale"))
        e <- sprintf("CAST(%s AS DOUBLE)", xtx_qi(D, s))
    } else if (is.numeric(s) && length(s) == 1L && is.null(dim(s))) {
        e <- xtx_sql_double(s)
    } else {
        stop("A remote design can only be scaled by a number, the name of a column or a DBI::SQL() expression.",
             call. = FALSE)
    }
    if (op == "/") e <- sprintf("(CAST(1 AS DOUBLE) / %s)", e)
    out <- if (is.null(D$scale)) e else sprintf("%s * %s", D$scale, e)
    ## check the expression now
    xtx_remote_extend(D, c(`__xtx_scale` = out), "scale")
    out
}

#' @export
Ops.xtx_remote <- function(e1, e2) {
    op <- .Generic
    if (nargs() == 1L) {
        if (op == "+") return(e1)
        if (op == "-") {
            e1$scale <- xtx_remote_scale_sql(e1, -1, "*")
            return(e1)
        }
    } else if (op %in% c("*", "/")) {
        d1 <- inherits(e1, "xtx_remote")
        d2 <- inherits(e2, "xtx_remote")
        if (d1 && d2) stop(sprintf("`%s` between two remote designs is not supported.", op), call. = FALSE)
        if (op == "/" && d2) stop("Only `design / scale` is supported, with the design first.", call. = FALSE)
        D <- if (d1) e1 else e2
        s <- if (d1) e2 else e1
        ## (a remote column of the design, as `e * D` or `D / w`: the design
        ## moves to the column's frame, see xtx_p_scale(); tapx)
        if (xtx_is_column(s)) {
            if (xtx_ncols(s) != 1L) stop("A remote design is scaled by a column with one column.", call. = FALSE)
            return(xtx_p_scale(D, if (op == "/") 1 / s else s))
        }
        D$scale <- xtx_remote_scale_sql(D, s, op)
        return(D)
    }
    stop(sprintf("`%s` is not supported for remote designs.", op), call. = FALSE)
}

## `col * D` and `D * col` with a remote column (design_column) and a
## remote design: both classes have Ops methods, so R asks which one to use
## (base::chooseOpsMethod(), R >= 4.3): the design's, which scales rows.
#' @export
chooseOpsMethod.xtx_remote <- function(x, y, mx, my, cl, reverse) inherits(y, "xtx_column")

## The design without its row scale.
xtx_remote_unscaled <- function(D) {
    D$scale <- NULL
    D
}

## Is `b` the same remote design as `a` (up to the row scale)?
xtx_remote_same <- function(a, b) identical(a$id, b$id)

#' @export
as.matrix.xtx_remote <- function(x, force = FALSE, ...) as.matrix(xtx_remote_collect(x, force))

#' @method crossprod xtx_remote
#' @export
crossprod.xtx_remote <- function(x, y = NULL, ...) {
    if (!inherits(x, "xtx_remote"))
        stop("Only `crossprod(design)`, `crossprod(design, design)` and `crossprod(design, \"column\")` are supported, with a design first.",
             call. = FALSE)
    ## crossprod() is a primitive: its dispatch drops argument names, so
    ## that `crossprod(D, weights = "w")` would arrive as `y = "w"`
    nm <- names(sys.call())
    if (...length() || (!is.null(nm) && any(!nm %in% c("", "x", "y"))))
        stop("`crossprod()` of a remote design takes only `x` and `y`; for weights, scale the rows: `crossprod(D, \"w\" * D)`.",
             call. = FALSE)
    if (is.null(y) || inherits(y, "xtx_remote")) {
        ## X' diag(s_x s_y) X
        if (!is.null(y) && !xtx_remote_same(x, y))
            stop("`crossprod(x, y)` of remote designs needs the same design, possibly with scaled rows.", call. = FALSE)
        ## (scales can read nodes: the products of the fit)
        if (!is.null(y)) x <- xtx_design_with(x, y, short = FALSE)
        w <- xtx_sql_times(x$scale, if (is.null(y)) x$scale else y$scale)
        return(xtx_remote_moments(xtx_remote_unscaled(x), w = w)$xtx)
    }
    ## X' diag(s_x) Y for columns: names, a DBI::SQL() expression, or the
    ## remote columns of the fitting core (primitives.R)
    if (!xtx_is_column(y) && !inherits(y, "SQL") && !is.character(y))
        stop("For a remote design, `y` must be another remote design, column names or a DBI::SQL() expression.",
             call. = FALSE)
    yc <- xtx_remote_cols(x, y, "y", several = TRUE)
    out <- xtx_remote_moments(xtx_remote_unscaled(yc$D), yc$sql, yc$D$scale, gram = FALSE)$xty
    colnames(out) <- yc$names
    out
}

xtx_remote_colsums <- function(x, group, reorder, mean) {
    if (!isTRUE(reorder))
        stop("Remote designs support only `reorder = TRUE`.", call. = FALSE)
    g <- NULL
    if (!is.null(group)) {
        gc <- xtx_remote_cols(x, group, "group", numeric = FALSE)
        x <- gc$D
        g <- gc$sql
    }
    w <- x$scale
    x <- xtx_remote_unscaled(x)
    cap <- xtx_compress_cap(x, group = g, op = "colsums")
    out <- if (!is.null(cap)) xtx_sql_colsums(x, g, w, mean, cap = cap)
    if (is.null(out)) xtx_sql_colsums(x, g, w, mean) else out
}

setOldClass("xtx_remote")

#' @rdname design_remote
#' @param x,object A remote design.
#' @param na.rm,dims Only the defaults are supported.
#' @param ... Not used.
#' @param group For `colSums()` and `colMeans()`, the name of a column of the
#'   remote data whose values define the groups, or a SQL expression over
#'   its columns (`DBI::SQL()`).
#' @param reorder Only `TRUE` (groups sorted) is supported.
#' @export
setMethod("colSums", "xtx_remote", function(x, na.rm = FALSE, dims = 1, ..., group = NULL, reorder = TRUE) {
    xtx_defaults_only(na.rm, dims)
    xtx_remote_colsums(x, group, reorder, mean = FALSE)
})

#' @rdname design_remote
#' @export
setMethod("colMeans", "xtx_remote", function(x, na.rm = FALSE, dims = 1, ..., group = NULL, reorder = TRUE) {
    xtx_defaults_only(na.rm, dims)
    xtx_remote_colsums(x, group, reorder, mean = TRUE)
})

#' @rdname design_remote
#' @export
setMethod("rowSums", "xtx_remote", function(x, na.rm = FALSE, dims = 1, ...) {
    xtx_defaults_only(na.rm, dims)
    xtx_column_label(xtx_p_mult(x, rep(1, ncol(x))), ".rowsum")
})

#' @rdname design_remote
#' @export
setMethod("rowMeans", "xtx_remote", function(x, na.rm = FALSE, dims = 1, ...) {
    xtx_defaults_only(na.rm, dims)
    xtx_column_label(xtx_p_mult(x, rep(1 / ncol(x), ncol(x))), ".rowmean")
})

## The rows of a remote design where a logical remote column, or a
## DBI::SQL() predicate over its data, is TRUE (see xtx_p_subset()). Numeric
## indices are refused: remote rows have no order.
#' @export
`[.xtx_remote` <- function(x, i, j, ..., drop = FALSE) {
    if (!missing(j)) stop("Remote designs keep all their columns: use `D[i, ]`.", call. = FALSE)
    if (missing(i)) return(x)
    if (is.numeric(i) || is.logical(i))
        stop("Rows of a remote design are selected by a DBI::SQL() predicate or a logical remote column, not by position: remote rows have no order.",
             call. = FALSE)
    xtx_p_subset(x, i)
}

## A sparse matrix of the collected rows, with the guard of as.matrix().
setAs("xtx_remote", "CsparseMatrix", function(from) as(xtx_remote_collect(from), "CsparseMatrix"))
setAs("xtx_remote", "dgCMatrix", function(from) as(xtx_remote_collect(from), "CsparseMatrix"))

## The in-memory design of the rows of a remote design, refused above
## getOption("xtx.remote_max_rows", 1e5) rows unless force.
xtx_remote_collect <- function(x, force = FALSE) {
    max_rows <- getOption("xtx.remote_max_rows", 1e5)
    if (!isTRUE(force) && nrow(x) > max_rows)
        stop(sprintf("The remote design has %s rows (more than %s): use `force = TRUE`, or raise the option `xtx.remote_max_rows`, to collect it in memory.",
                     xtx_fmt_count(nrow(x)), xtx_fmt_count(max_rows)), call. = FALSE)
    sc <- x$scale
    x <- xtx_remote_unscaled(x)
    cols <- xtx_remote_inputs(x)
    if (!is.null(sc)) {
        x <- xtx_remote_extend(x, c(`__xtx_scale` = sc), "scale")
        cols <- c(cols, "__xtx_scale")
    }
    chunk <- xtx_query(x, xtx_chunk_sql(x, cols))
    cd <- xtx_chunk_design(x, chunk)
    if (is.null(sc)) cd$design else chunk$`__xtx_scale`[cd$keep] * cd$design
}

#' @rdname design_remote
#' @param sql For `summary()`: also print the generated SQL (the query that
#'   computes the model variables, and the aggregates of `crossprod()`).
#' @export
summary.xtx_remote <- function(object, sql = FALSE, ...) {
    x <- object
    checkmate::assert_flag(sql)
    tab <- data.frame(
        term = vapply(x$terms, function(t) t$label, ""),
        sql = vapply(x$terms, function(t) t$sql_ok, NA),
        reason = vapply(x$terms, function(t) t$reason %||% "", ""),
        stringsAsFactors = FALSE)
    attr(tab, "engine") <- if (is.null(x$evaluated)) "sql" else sprintf("sql, on columns computed in R (%s)", x$evaluated)
    attr(tab, "scale") <- x$scale
    if (sql) {
        attr(tab, "sql") <- c(paste0("-- model variables\n", xtx_sql_base(x)),
                              paste0("-- crossprod\n", xtx_sql_moments(x, sql_only = TRUE)))
    }
    class(tab) <- c("summary.xtx_remote", class(tab))
    tab
}

#' @export
print.summary.xtx_remote <- function(x, ...) {
    cat("engine: ", attr(x, "engine"), "\n", sep = "")
    if (!is.null(attr(x, "scale"))) cat("rows scaled by: ", attr(x, "scale"), "\n", sep = "")
    sql <- attr(x, "sql")
    tab <- x
    attr(tab, "engine") <- attr(tab, "sql") <- attr(tab, "scale") <- NULL
    class(tab) <- "data.frame"
    print(tab, row.names = FALSE)
    for (q in sql) cat("\n", q, "\n", sep = "")
    invisible(x)
}

## ---- the in-memory design of rows read from the database ----

## The query reading columns of the rows used.
xtx_chunk_sql <- function(D, cols) {
    sel <- if (length(cols)) paste0("s.", xtx_qi(D, cols), collapse = ", ") else "1 AS xtx_one"
    sprintf("SELECT %s FROM (%s) AS s WHERE %s", sel, D$src$select, D$filter)
}

## The in-memory design of rows read from the database (a sample for the
## prototype of a design, or all of them for as.matrix()): rows with a
## missing model variable (response included) are dropped, as na.omit()
## drops them, unless the design keeps them (na.pass).
## Returns list(design, keep = rows of the chunk kept, y = response).
xtx_chunk_design <- function(D, chunk) {
    plan <- D$plan
    na <- if (identical(D$na, "pass")) stats::na.pass else stats::na.omit
    mf <- design_frame(D$terms_full, chunk, xlev = plan$xlevels, na.action = na)
    ## variables of another type than in training: the error of predict.lm()
    if (!is.null(plan$data_classes)) stats::.checkMFClasses(plan$data_classes, mf)
    design <- design_matrix(plan$terms, mf, contrasts.arg = plan$contrasts)
    if (isTRUE(plan$drop_intercept)) design <- xtx_drop_intercept(design)
    omit <- attr(mf, "na.action")
    keep <- seq_len(nrow(chunk))
    if (length(omit)) keep <- keep[-omit]
    y <- if (is.null(D$response)) NULL else unname(stats::model.response(mf))
    list(design = design, keep = keep, y = y)
}
