library(linmod2)
set.seed(1234)

cmp <- function(a, b) stopifnot(all.equal(as.matrix(a), b))

m = 20
n = 5
x = matrix(rnorm(m*n), m, n)
y = rnorm(m)



mdl1 = .lm.fit(x, y)
mdl2 = .lm_fit(x, y)

cmp(mdl1$coefficients, mdl2$coefficients)
cmp(mdl1$residuals, mdl2$residuals)
cmp(x %*% mdl1$coefficients, mdl2$fitted.values)



x = t(x)
y = y[1:n]
mdl1 = .lm.fit(x, y)
mdl2 = .lm_fit(x, y)

#### NOTE: the coefficients can't possibly be equal because R will fit a model
#### with only the first `rank==m` coefficients, while we use all `n`.
### cmp(mdl1$coefficients, mdl2$coefficients)
cmp(mdl1$residuals, mdl2$residuals)
cmp(x %*% mdl1$coefficients, mdl2$fitted.values)
