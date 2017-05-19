library(linmod2)
set.seed(1234)

cmp <- function(a, b) stopifnot(all.equal(as.matrix(a), b))
fit <- function(x, y)
{
  mdl1 <<- lm.fit(x, y)
  mdl2 <<- lm_fit(x, y)
}

m = 20
n = 5
x = matrix(rnorm(m*n), m, n)
y = rnorm(m)


# ---------------- Tall and skinny ----------------

# 1 RHS
fit(x, y)

cmp(mdl1$coefficients, mdl2$coefficients)
cmp(x %*% mdl1$coefficients, x %*% mdl2$coefficients)
cmp(mdl1$fitted.values, mdl2$fitted.values)
cmp(mdl1$residuals, mdl2$residuals)



# ---------------- Short and wide ----------------

#### NOTE: the coefficients can't possibly be equal because R will fit a model
#### with only the first 'rank==m' coefficients, while we use all 'n'.

# 1 RHS
x = t(x)
y = y[1:n]
fit(x, y)

cmp(mdl1$residuals, mdl2$residuals)
cmp(mdl1$fitted.values, mdl2$fitted.values)
