#import "/.calepin/calepin.typ" as calepin
#show: calepin.document

#set document(title: "The rows of a remote design and columns over them, as a SQL relation")
#metadata((
  title: "design_relation",
  pagetitle: "The rows of a remote design and columns over them, as a SQL relation — design_relation",
  aliases: ("design_relation",),
  "source-files": ("R/rows.R",),
)) <website-metadata>

#calepin.setup(echo: true, eval: true, results: "verbatim", warning: false, message: false)

#calepin.chunk(echo: false, results: "hide")[
```r
library(xtx)
```
]

= The rows of a remote design and columns over them, as a SQL relation

== Description

`design_relation()` gives what a query needs to read columns over the rows of a remote design: the connection, the source to select from, the filter of the rows and the SQL expression of each column\. It is for computations that stay in the database \(tables of sums within the levels of a variable with too many levels to read, say\), built on the columns of a design: products with the design \(`x %*% b`\), responses, weights and keys \(#link("design_key.typ")[`design_key()`]\)\.

== Usage

```r
#| eval: false
design_relation(x, cols)
```

== Arguments

#table(
  columns: 2,
  table.header([Argument], [Description]),
  [`x`], [A remote design\.],
  [`cols`], [A named list of columns over the rows of `x` \(#link("design_column.typ")[`design_column()`]s\), each with one or several SQL columns\.],
)

== Value

A list with `con` \(the DBI connection\), `select` and `filter` \(strings\) and `sql`, a list named as `cols`: the SQL expressions of each column, a character vector\.

== Details

The rows are those of `SELECT ... FROM (<select>) AS s WHERE <filter>`, and each expression names the columns of `s`\. Temporary tables that the expressions or the filter name must exist when the query runs\.

== Examples

```r
if (requireNamespace("duckdb", quietly = TRUE)) {
    con <- DBI::dbConnect(duckdb::duckdb())
    DBI::dbWriteTable(con, "d", data.frame(y = rnorm(10), x = rnorm(10), id = rep(1:5, 2)))
    D <- design_matrix(~ x, con, table = "d")
    r <- design_relation(D, list(xb = D %*% c(1, 2), id = design_key(D, "id")))
    DBI::dbGetQuery(r$con, sprintf("SELECT %s AS id, SUM(%s) AS s FROM (%s) AS s WHERE %s GROUP BY 1",
                                   r$sql$id, r$sql$xb, r$select, r$filter))
    DBI::dbDisconnect(con)
}
```
