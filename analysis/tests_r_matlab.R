library(bench)
library(causalXtreme)
library(profvis)
library(Rcpp)
library(tictoc)
devtools::load_all(".")

# First tests ####
rm(list = ls())
sourceCpp("analysis/cpp/misc.cpp")

# Import data
x <- read.csv("analysis/temp_csv/X_mat.csv", header = F)
x <- t(x)
x <- x[, 1:30]

# Try kernel ica
m <- NROW(x)
N <- NCOL(x)

if (N < 1000){
  sigma <- 1
  kappa <- 2e-2

} else{
  sigma <- 1/2
  kappa <- 2e-3
}

eta <- kappa * 1e-2
i <- 1
tic()
res <- chol_gauss(x[i,, drop = FALSE] / sigma, 1, N * eta);
toc()


# Tests Lingam ####
rm(list = ls())
sourceCpp("analysis/cpp/misc.cpp")

# Simulate data
n <- 100
p <- 5
set.seed(42)
X <- simulate_data(n, p, 2/(p - 1), has_confounder = TRUE)


# Run direct lingam
tic()
order_ling <- direct_lingam_search(X$dataset)
toc()
est_g <- causalXtreme:::caus_order_to_dag(order_ling)
est_cpdag <- causalXtreme:::dag_to_cpdag(est_g)
ling_res <- list(est_g = est_g, est_cpdag = est_cpdag)

# Run EASE
tic()
ease_res <- causal_discovery(X$dataset, "ease")
toc()

# Run lingam ICA
tic()
ica_res <- causal_discovery(X$dataset, "order_lingam")
toc()

# Compare performances
causal_metrics(X, ling_res)
causal_metrics(X, ease_res)
causal_metrics(X, ica_res)


# Kernel ICA C++ tests ####
rm(list = ls())

# Test 1
sourceCpp("analysis/cpp/misc.cpp")
x <- rbind(1:3)
a <- chol_gaussc(x, 1, 1)
b <- chol_gauss(x, 1, 1)
all.equal(a, b)

x <- rbind(1:500)
bench::mark(
  chol_gaussc(x, 1, 1),
  chol_gauss(x, 1, 1),
  check = FALSE
)


# Test 2
sourceCpp("analysis/cpp/misc.cpp")
x <- rbind(1:3, 4:6)
a <- contrast_icac(x, 1, .3, 1)
chol_gaussc(x[2, , drop = FALSE], 1, .3*3)
G <- a$res$G
G <- scale(G, scale = FALSE)
eigen(crossprod(G))

xx <- crossprod(x)
eigen(xx)
symm_eigen(xx)

x <- c(1, NA, -1, 2)
which2(x, 0)


# Pairwise Lingam C++ tests ####
rm(list = ls())
sourceCpp("analysis/cpp/pwling_funcs.cpp")

set.seed(42)
x <- sample(1:1e3, 500)
y <- sample(1:1e3, 500)

bench::mark(
  setdiff(x, y),
  mysetdiff(x, y),
  check = FALSE
)


# Find index AND compute Residuals C++ tests ####
rm(list = ls())
sourceCpp("analysis/cpp/pwling_funcs.cpp")

n <- 1000
p <- 30
set.seed(42)
Xsim <- simulate_data(n, p, .8)
X <- t(Xsim$dataset)
M <- matrix(-1, ncol = p, nrow = p)
diag(M) <- 0
candidates <- 1:p
U_K <- 1:p

# Compute residuals
R <- computeRc(X, candidates - 1, U_K - 1, M)
R2 <- computeR(X, candidates, U_K, M)

bench::mark(
  computeRc(X, candidates - 1, U_K - 1, M),
  computeR(X, candidates, U_K, M)
)


# Find index
findindexc(X, candidates-1, U_K-1)
findindex(X, R, candidates, U_K)

bench::mark(
  findindexc(X, candidates-1, U_K-1),
  findindex(X, R, candidates, U_K)
)


# Tests Lingam C++ speed ####
rm(list = ls())
sourceCpp("analysis/cpp/pwling_funcs.cpp")
devtools::load_all(".")
set.seed(42)
X5 <- simulate_data(10000, 30, 0.2, has_uniform_margins = FALSE)
write.csv(X5$dataset, "analysis/temp_csv/X5.csv", row.names = FALSE)

tic()
direct_lingam_search(X5$dataset)
toc()

tic()
direct_lingam_search_r(X5$dataset)
toc()

profvis(direct_lingam_search(X5$dataset))
profvis(direct_lingam_search_r(X5$dataset))

