
<!-- README.md is generated from README.Rmd. Please edit that file -->
causalXtreme
============

<!-- badges: start -->
[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing) [![Travis build status](https://api.travis-ci.org/nicolagnecco/causalXtreme.svg?branch=master)](https://travis-ci.org/nicolagnecco/causalXtreme) [![codecov](https://codecov.io/gh/nicolagnecco/causalXtreme/branch/master/graph/badge.svg)](https://codecov.io/gh/nicolagnecco/causalXtreme) <!-- badges: end -->

The goal of causalXtreme is to provide an interface to perform causal discovery in linear structural equation models (SEM) with heavy-tailed noise. For more details see Gnecco et al. (2019 <https://arxiv.org/abs/1908.05097>).

Installation
------------

<!-- You can install the released version of causalXtreme from [CRAN](https://CRAN.R-project.org) with: -->
<!-- ``` r -->
<!-- install.packages("causalXtreme") -->
<!-- ``` -->
You can install the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("nicolagnecco/causalXtreme")
```

Example
-------

Let us generate 500 observations from a SEM of two Student-t variables, *X*<sub>1</sub> and *X*<sub>2</sub>, with 1.5 degrees of freedom (i.e., heavy-tailed). When the function `simulate_data` is called with the default values, it returns a list containing:

-   the simulated `dataset` represented as a matrix of size *n* × *p*. Here, *n* = 500 is the number of observations and *p* = 2 is the number of variables,
-   the underlying directed acyclic graph (DAG) represented as an adjacency matrix `dag` of size *p* × *p*.

``` r
library(causalXtreme)

# basic example code
set.seed(1)
sem <- simulate_data(n = 500, p = 2, prob_connect = 0.5,
                     distr = "student_t", tail_index = 1.5)
```

At this point, we can look at the randomly simulated DAG.

``` r
sem$dag
#>      [,1] [,2]
#> [1,]    0    1
#> [2,]    0    0
```

We interpret the adjacency matrix as follows. Loosely speaking, we say that variable *X*<sub>*i*</sub> causes variable *X*<sub>*j*</sub> if the entry (*i*, *j*) of the adjacency matrix is equal to 1. We see that the first variable *X*<sub>1</sub> causes the second variable *X*<sub>2</sub>, since the entry (1, 2) of the matrix `sem$dag` is equal to 1. We can plot the simulated dataset.

``` r
plot(sem$dataset, pch = 20,
     xlab = "X1", ylab = "X2")
```

<img src="man/figures/README-plotdata-1.png" width="100%" />

At this point, we can estimate the causal direction between *X*<sub>1</sub> and *X*<sub>2</sub> by computing the *causal tail coefficients* *Γ*<sub>12</sub> and *Γ*<sub>21</sub> (see Gnecco et al. 2019, Definition 1).

``` r
X1 <- sem$dataset[, 1]
X2 <- sem$dataset[, 2]

# gamma_12
causal_tail_coeff(X1, X2)
#> [1] 0.9523333

# gamma_21
causal_tail_coeff(X2, X1)
#> [1] 0.4816667
```

We see that the coefficient *Γ*<sub>12</sub> ≈ 1 (entry (1, 2) of the matrix) and *Γ*<sub>21</sub> &lt; 1 (entry (2, 1) of the matrix). This is evidence for a causal relationship from *X*<sub>1</sub> to *X*<sub>2</sub>.

We can also run the *extremal ancestral search* (EASE) algorithm, based on the causal tail coefficients (see Gnecco et al. 2019, Section 3.1). The algorithm estimates from the data a *causal order* of the DAG.

``` r
ease(dat = sem$dataset)
#> [1] 1 2
```

In this case, we can see that the estimated causal order is correct, since *X*<sub>1</sub> (the cause) is placed before *X*<sub>2</sub> (the effect).

References
----------

Gnecco, Nicola, Nicolai Meinshausen, Jonas Peters, and Sebastian Engelke. 2019. “Causal Discovery in Heavy-Tailed Models.” *arXiv Preprint arXiv:1908.05097*.
