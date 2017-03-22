library(linmod2)
library(rbenchmark)


m = 2000
n = 500
x = matrix(rnorm(m*n), m, n)
y = rnorm(m)

benchmark(.lm.fit(x, y), .lm_fit(x, y), replications=25)



m = 500
n = 2000
x = matrix(rnorm(m*n), m, n)
y = rnorm(m)

benchmark(.lm.fit(x, y), .lm_fit(x, y), replications=25)
