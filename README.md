# linmod2

* **Version:** 0.1-0
* **Status:** [![Build Status](https://travis-ci.org/wrathematics/linmod2.png)](https://travis-ci.org/wrathematics/linmod2)
* **License:** [![License](http://img.shields.io/badge/license-BSD%202--Clause-orange.svg?style=flat)](http://opensource.org/licenses/BSD-2-Clause)
* **Author:** Drew Schmidt
* Project home: https://github.com/wrathematics/linmod2
* Bug reports: https://github.com/wrathematics/linmod2/issues
* Email: wrathematics .AT. gmail .DOT. com
* Twitter: @wrathematics

Fast, portable, numerically stable linear and generalized linear model fitters for C and R.

Progress:
* lm fitter (~75% done)
* glm fitter (~5% done)


## Background

This is a ground-up rewrite of R's linear and generalized linear model fitters.  This roughly corresponds to replacing `lm.fit()` and `glm.fit()`. The package is also a rewrite of the [linmod package](https://github.com/wrathematics/linmod).

The primary goal of the package is to offer high quality lm/glm fitters in a standalone shared library written in C.  This will be bundled with the R package *linmod2* and given appropriate bindings that emulate (when reasonable) R's model fitter functions.  The **linmod** package actually mostly does this. However, it's written mostly in Fortran, which creates cross-platform issues for some people, and massive code modularity and build issues for me, the author (the package requires cmake, for example).  Additionally, there were some poor choices made in developing the original **linmod** package that are a bit hard to explain, but just trust me. This rewrite should solve all of these problems.

There are several reasons why you should be interested in such a thing.  First, R's lm fitter (which is the driver in its glm fitter) is not very sophisticated.  It uses the ancient LINPACK to fit the numerical solution.  This creates performance, as well as numerical issues.  What's more, R has no ability to take advantage of fitting models when you *know* that the model matrix is full rank; we can.  Additionally, R exclusively uses a QR decomposition, even when m>>n.  For numerical and performance reasons, it is arguably better to use a different but related factorization, LQ.  Finally, R's glm fitter is mostly written in R, which is not a good idea if performance is paramount.



## Installation

<!-- To install the R package, run:

```r
install.package("coop")
``` -->

The development version is maintained on GitHub, and can easily be installed by any of the packages that offer installations from GitHub:

```r
### Pick your preference
devtools::install_github("wrathematics/linmod2")
ghit::install_github("wrathematics/linmod2")
remotes::install_github("wrathematics/linmod2")
```



## Example Usage

```r
library(linmod2)

n <- 10
p <- 3
x <- matrix(rnorm(n*p), n, p)
y <- rnorm(n)

.lm_fit(x, y)
```



## Benchmarks

All of these benchmarks can be found in the source tree of this package, under `linmod2/inst/benchmarks`.  All tests performed using:

* R 3.3.3
* OpenBLAS
* gcc 6.2.0
* 4 cores of a Core i5-2500K CPU @ 3.30GHz


```r
library(linmod2)
library(rbenchmark)


m = 10000
n = 500
x = matrix(rnorm(m*n), m, n)
y = rnorm(m)

benchmark(.lm.fit(x, y), .lm_fit(x, y), replications=25)
##            test replications elapsed relative
## 2 .lm_fit(x, y)            5   1.133    1.000
## 1 .lm.fit(x, y)            5   9.497    8.382
```
