#import "/.calepin/calepin.typ" as calepin
#show: calepin.document

#set document(title: [Standard errors])
#metadata((title: "Standard errors")) <website-metadata>

#import "/benchmarks/bench.typ": bench, bench-table, bench-code

#calepin.setup(echo: true, eval: true, results: "verbatim", warning: false, message: false)

#calepin.chunk(echo: false, results: "hide")[
```r
options(width = 200)
```
]

#title()

```r
library(xtx)
library(xtx.fit)
d <- data.frame(
    y = c(0.3, -1.2, 0.8, 1.5, -0.4, 0.9, -2.1, 0.2),
    x = c(-1.5, 0.2, 1.1, -0.3, 2.0, 0.7, -0.9, 1.4),
    f = factor(c("a", "b", "c", "a", "b", "c", "a", "b")),
    cl = c(1, 1, 2, 2, 3, 3, 4, 4))
```

`lm_xtx()` and `glm_xtx()` take `vcov = "iid"` (the default) or `vcov = ~cl`, for standard errors clustered by a variable. The last section computes heteroskedasticity-consistent and clustered covariance matrices for data that stay in a database.

= Cluster-robust standard errors by hand

The cluster-robust covariance matrix is the sandwich `B M B`, with bread `B = (X'X)^-1` and meat `M = S'S`, where row `c` of `S` sums `e_i x_i` over the observations in cluster `c` (without a small-sample correction). The same function works on the dense model matrix and on the design matrix.

```r
base_cluster_vcov <- function(X, e, cluster, bread) {
    S <- rowsum(e * X, cluster)
    bread %*% crossprod(S) %*% bread
}
fast_cluster_vcov <- function(D, e, cluster, bread) {
    S <- colSums(e * D, group = cluster)
    bread %*% crossprod(S) %*% bread
}
fit <- lm(y ~ x + f, d)
D <- design_matrix(fit, newdata = d)
bread <- vcov(fit) / sigma(fit)^2
a <- base_cluster_vcov(model.matrix(fit), residuals(fit), d$cl, bread)
b <- fast_cluster_vcov(D, residuals(fit), d$cl, bread)
all.equal(a, b, check.attributes = FALSE)
```

The benchmark runs on top of the #link("../benchmarks.html#setup")[benchmark setup], in a fresh R process, on 200,000 rows and 119 columns, on one thread.

#bench("cluster")

= Robust standard errors on a database

```r
library(xtx)
library(xtx.fit)
con <- DBI::dbConnect(duckdb::duckdb())
set.seed(1)
n <- 10000
d <- data.frame(
    id = seq_len(n),
    x = rnorm(n),
    z = sample(1:5, n, TRUE),
    f = sample(c("a", "b", "c"), n, TRUE),
    w = runif(n, 0.5, 2),
    cl = sample(1:50, n, TRUE))
d$y <- d$x + d$z + rnorm(n) * (1 + abs(d$x))
d$k <- rpois(n, exp(0.2 * d$z))
DBI::dbWriteTable(con, "d", d)
```

`residuals()` of a remote fit returns the remote design with two more columns, `.fitted` and `.resid`, which can scale its rows. The heteroskedasticity-consistent (HC1) and clustered (CR1, as in `sandwich::vcovCL()`) covariance matrices follow:

```r
fit <- lm_xtx(y ~ x + f, con, table = "d")
R <- residuals(fit)
bread <- vcov(fit) / sigma(fit)^2
n <- nobs(fit)
k <- fit$rank

meat <- crossprod(R, DBI::SQL('".resid" * ".resid"') * R)
V_hc1 <- n / (n - k) * bread %*% meat %*% bread

scores <- colSums(".resid" * R, group = "cl")
G <- nrow(scores)
V_cl <- G / (G - 1) * (n - 1) / (n - k) * bread %*% crossprod(scores) %*% bread
sqrt(cbind(HC1 = diag(V_hc1), CR1 = diag(V_cl)))
```

The same computation on the dense model matrix:

```r
l <- lm(y ~ x + f, d)
X <- model.matrix(l)
e <- residuals(l)
B <- solve(crossprod(X))
all.equal(V_hc1, n / (n - k) * B %*% crossprod(X, e^2 * X) %*% B, check.attributes = FALSE)
S <- rowsum(e * X, d$cl)
all.equal(V_cl, G / (G - 1) * (n - 1) / (n - k) * B %*% crossprod(S) %*% B, check.attributes = FALSE)
```

#calepin.chunk(echo: false, results: "hide")[
```r
DBI::dbDisconnect(con)
```
]
