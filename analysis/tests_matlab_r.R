# Tests ####
rm(list = ls())
source("analysis/kernel_ica_test.R")
x <- read.csv("analysis/X_mat.csv", header = F)
x <- t(x)
x <- x[, 1:100]

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

J = contrast_ica(contrast,x,kparam)

etas <- kparam$etas
i <- 1
tic()
res <- chol_gauss(x[i,1:100, drop = FALSE] / kparam$sigmas[i], 1, N * etas[i]);
toc()

