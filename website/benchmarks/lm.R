library(xtx.fit)
bench::mark(
    "lm()"     = coef(lm(fo, big)),
    "lm_xtx()" = coef(lm_xtx(fo, big)),
    iterations = iters(10), check = close)
