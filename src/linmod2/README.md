# linmod2

* **Version:** 0.1-0
* **License:** [![License](http://img.shields.io/badge/license-BSD%202--Clause-orange.svg?style=flat)](http://opensource.org/licenses/BSD-2-Clause)
* **Author:** Drew Schmidt
* Project home: https://github.com/wrathematics/linmod2
* Bug reports: https://github.com/wrathematics/linmod2/issues
* Email: wrathematics .AT. gmail .DOT. com
* Twitter: @wrathematics


**linmod2** is an efficient, portable, numerically stable C99 library for fitting linear and generalized linear models.  The primary functions `lm_fit()` and `glm_fit()` are meant to operate similarly to the R language's own `lm.fit()` and `glm.fit()` functions.

It is a from-the-ground-up rewrite of **linmod**, which is written in Fortran.  While the libraries are meant to work well with R (both in the API and numerical senses), they contain completely original code and are not affiliated with the R project in any way.

**linmod2** is also permissively licensed under the BSD 2-clause license.



## Installation and Use

The library is header-only (frequently against my better judgment), so no building or making is required. Simply copy the contents of **linmod2** and include the relevant top-level header file(s).

Currently the library lacks documentation, as its use in the identically named R package is taking priority at the moment.  As those bindings are developed, the API is subject to change.
