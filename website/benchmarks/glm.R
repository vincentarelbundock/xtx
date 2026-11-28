library(xtx.fit)
big$yb <- as.integer(big$y + big$x > 0)
fb <- yb ~ x * f1 + f2
bench::mark(
    "glm()"     = coef(glm(fb, binomial, big)),
    "glm_xtx()" = coef(glm_xtx(fb, binomial, big)),
    iterations = iters(5), check = close)
