library(linmod2)
library(rbenchmark)
set.seed(1234)

cols = cols <- c("test", "replications", "elapsed", "relative")
reps = 5

m = 10000
n = 500
x = matrix(rnorm(m*n), m, n)
y = rnorm(m)


#benchmark(.lm.fit(x, y), .lm_fit(x, y), .lm_fit_minimal(x, y), replications=reps, columns=cols)
benchmark(lm.fit(x, y), lm_fit(x, y), replications=reps, columns=cols)
