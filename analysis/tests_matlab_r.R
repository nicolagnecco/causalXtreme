# Tests ####
rm(list = ls())
library(tictoc)
library(profvis)

# Import data
source("analysis/kernel_ica_test.R")
sourceCpp("analysis/chol_gaussc.cpp")
sourceCpp("analysis/eigenc.cpp")

x <- read.csv("analysis/X_mat.csv", header = F)
x <- t(x)
x <- x[, 1:500]

# Try kernel ica
m <- NROW(x)
N <- NCOL(x)

contrast='kgv';

if (N < 1000){
  sigma=1;
  kappa=2e-2;

} else{
  sigma = 1/2;
  kappa = 2e-3;
}

kernel='gaussian';

mc=m;
kparam <- list()
kparam$kappas=kappa*rep(1,mc);
kparam$etas=kappa*1e-2*rep(1,mc);
kparam$neigs=N*rep(1,mc);
kparam$nchols=N*rep(1,mc);
kparam$kernel=kernel;
kparam$sigmas=sigma*rep(1,mc);

bench::mark(
  contrast_ica(contrast,x,kparam),
  contrast_icac(contrast,x,kparam)
)


etas <- kparam$etas
i <- 1
tic()
res <- chol_gauss(x[i,, drop = FALSE] / kparam$sigmas[i], 1, N * etas[i]);
toc()


# Tests Lingam ####
rm(list = ls())
library(tictoc)
library(profvis)
library(causalXtreme)
source("analysis/kernel_ica_test.R")
source("analysis/direct_lingam_test.R")
sourceCpp("analysis/chol_gaussc.cpp")
sourceCpp("analysis/eigenc.cpp")

n <- 2000
p <- 8
set.seed(42)
X <- simulate_data(n, p, .2, has_confounder = FALSE)

tic()
direct_lingam(X$dataset)
toc()

write.csv(X$dataset, "analysis/lingamX2.csv", row.names = FALSE)

# Exec time
profvis(direct_lingam(X$dataset))


# Performance
est_g <- causalXtreme:::caus_order_to_dag(order_ling)
est_cpdag <- causalXtreme:::dag_to_cpdag(est_g)
ling_res <- list(est_g = est_g, est_cpdag = est_cpdag)

tic()
ease_res <- causal_discovery(X$dataset, "ease")
toc()

ica_res <- causal_discovery(X$dataset, "order_lingam")


causal_metrics(X, ling_res)
causal_metrics(X, ease_res)
causal_metrics(X, ica_res)


## C++ tests ####
rm(list = ls())
library(bench)
library(Rcpp)
sourceCpp("analysis/chol_gaussc.cpp")
source("analysis/kernel_ica_test.R")

x <- rbind(1:3, 4:6)
chol_gaussc(x, 1, 1)
chol_gauss(x, 1, 1)

x <- rbind(1:1e6, 1:1e6)
bench::mark(
  chol_gaussc(x, 1, 1),
  chol_gauss(x, 1, 1)
)

sourceCpp("analysis/eigenc.cpp")
x <- matrix(rnorm(1e6), nrow = 1e3)

xx <- crossprod(x)
bench::mark(
  symm_eigen(xx),
  symm_eigen2(xx),
  eigen(xx), # FASTEST
  check = FALSE
)


sourceCpp("analysis/contrast_ica.cpp")
contrast_ica(xx, 1, 1, 1)
