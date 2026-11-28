#' Rows of a design, with the columns that go with them
#'
#' `design_subset()` keeps the rows of a design where a logical column is
#' `TRUE`, as `x[keep, ]` does, and carries other columns over the rows to
#' the new design: in memory their values on the rows kept; remotely the
#' same expressions over the new design (a subset is another design, and
#' its columns do not combine with those of the design it was made from;
#' only columns of that design can be carried over). With `keep = NULL` no
#' row is dropped and nothing changes: the design is returned as it is,
#' with the columns (a remote design and its columns always combine).
#'
#' @param x A design: a `design_matrix` or a remote design.
#' @param keep The rows to keep: a logical vector in memory (missing values
#'   count as `FALSE`), a logical [design_column()] of `x` remotely; or
#'   `NULL`, to keep every row (the design and the columns as they are).
#' @param cols Columns over the rows of `x`: a list, possibly nested, of
#'   columns (vectors, matrices or data frames in memory, [design_column()]s
#'   remotely) and `NULL`s. With `keep = NULL` it may also hold bare group
#'   codes (see [key_codes()]), which are carried over; otherwise groupings
#'   are refused, since codes change with the rows.
#' @return A list with `x` (the design) and `cols` (the columns, in the
#'   shape of `cols`).
#' @examples
#' d <- data.frame(y = rnorm(10), x = rnorm(10))
#' D <- design_matrix(y ~ x, d)
#' s <- design_subset(D, d$x > 0, list(y = d$y))
#' nrow(s$x) == length(s$cols$y)
#' @keywords internal
#' @export
design_subset <- function(x, keep = NULL, cols = list()) {
    checkmate::assert_multi_class(x, c("xtx_matrix", "xtx_remote"))
    if (is.null(keep)) return(list(x = x, cols = xtx_cols_map(cols, function(col) xtx_c_on(col, x))))
    x2 <- xtx_p_subset(x, keep)
    walk <- function(v) {
        if (is.list(v) && is.null(oldClass(v))) lapply(v, walk) else xtx_c_subset(v, keep, x2)
    }
    list(x = x2, cols = walk(cols))
}

#' The rows of a remote design and columns over them, as a SQL relation
#'
#' `design_relation()` gives what a query needs to read columns over the
#' rows of a remote design: the connection, the source to select from, the
#' filter of the rows and the SQL expression of each column. It is for
#' computations that stay in the database (tables of sums within the levels
#' of a variable with too many levels to read, say), built on the columns
#' of a design: products with the design (`x %*% b`), responses, weights
#' and keys ([design_key()]).
#'
#' The rows are those of
#' `SELECT ... FROM (<select>) AS s WHERE <filter>`, and each expression
#' names the columns of `s`. Temporary tables that the expressions or the
#' filter name must exist when the query runs.
#'
#' @param x A remote design.
#' @param cols A named list of columns over the rows of `x`
#'   ([design_column()]s), each with one or several SQL columns.
#' @return A list with `con` (the DBI connection), `select` and `filter`
#'   (strings) and `sql`, a list named as `cols`: the SQL expressions of
#'   each column, a character vector.
#' @examples
#' if (requireNamespace("duckdb", quietly = TRUE)) {
#'     con <- DBI::dbConnect(duckdb::duckdb())
#'     DBI::dbWriteTable(con, "d", data.frame(y = rnorm(10), x = rnorm(10), id = rep(1:5, 2)))
#'     D <- design_matrix(~ x, con, table = "d")
#'     r <- design_relation(D, list(xb = D %*% c(1, 2), id = design_key(D, "id")))
#'     DBI::dbGetQuery(r$con, sprintf("SELECT %s AS id, SUM(%s) AS s FROM (%s) AS s WHERE %s GROUP BY 1",
#'                                    r$sql$id, r$sql$xb, r$select, r$filter))
#'     DBI::dbDisconnect(con)
#' }
#' @export
design_relation <- function(x, cols) {
    checkmate::assert_class(x, "xtx_remote")
    checkmate::assert_list(cols, names = "unique", min.len = 1L)
    for (nm in names(cols))
        if (!xtx_is_column(cols[[nm]])) stop(sprintf("`cols$%s` must be a column of the design.", nm), call. = FALSE)
    all <- do.call(xtx_c_cbind, unname(cols))
    u <- unclass(all)
    fr <- xtx_frame_ready(u$frame, u$sql)
    k <- vapply(cols, xtx_ncols, 1L)
    list(con = xtx_con(fr), select = fr$src$select, filter = fr$filter,
         sql = unname(split(u$sql, rep(seq_along(cols), k))) |> stats::setNames(names(cols)))
}
