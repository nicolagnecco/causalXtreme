library(Rcpp)

rm(list = ls())
sourceCpp("analysis/cpp/pwling_funcs.cpp")
set.seed(21)
X <- simulate_data(30, 5, 0.2, has_uniform_margins = FALSE)$dataset
direct_lingam_search(X)
