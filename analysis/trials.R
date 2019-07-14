# Understand causal tail coefficient with both tails
n <- 11
x1 <- runif(n)
r1 <- rank(x1)
x2 <- runif(n)
r2 <- rank(x2)
k <- 2
plot(r1 > n - k/2 | r1 <= k/2)
plot(abs(2 * r1 - n) > n - k)
plot(r1 > n - k)
sum(r1 > n - k / 2 | r1 <= k / 2)
table(abs(r2 - (n + 1)/2))
table(r2)
1 / (k * n) * sum(2 * abs(r2[r1 > n - k / 2 | r1 <= k / 2] - (n + 1)/2))

set.seed(1991)
n <- 5e2
X1 <- rt(n, df = 2.5)
X2 <- X1 + rt(n, df = 2.5)
X3 <- X1 + rt(n, df = 2.5)
dat <- cbind(X1, X2, X3)
compute_gamma_matrix(dat, both_tails = FALSE)
compute_gamma_matrix(dat, both_tails = TRUE)
pairs(dat)

# Rescale betas so that SNR = (1 - w)/w
p <- 3
g <- rbind(c(0, 1, 1), c(0, 0, 1), c(0, 0, 0))
A <- random_coeff(g, two_intervals = FALSE)
B <- t(A)
path <- solve(diag(p) - B)
has_parent <- apply(g, 2, sum) != 0
w <- .4
diag(path)[has_parent] <- sqrt(w)

path_temp <- path^2
diag(path_temp) <- 0
s <- apply(path_temp, 1, sum)
s[s != 0] <-  sqrt(1 / (s[s != 0]) * (1 - w))
s[s == 0] <- 1
mm <- matrix(rep(s, p), ncol = p)
diag(mm) <- 1
out <- path * mm
apply(out^2, 1, sum)
A_mod <- t(diag(p) - solve(out))
A
A_mod

# Check simulate data
n <- 1e4
p <- 50

X <- simulate_data(n, p, prob_connect = 1.5 / (p - 1), distr = "student_t", tail_index = 2,
                   has_confounder = T, has_uniform_margins = F)
pairs(X$dataset)



# Check Psi coefficient
n <- 1e3
p <- 15
X <- simulate_data(n, p, prob_connect = 1.5 / (p - 1), distr = "student_t", tail_index = 3.5,
                   has_confounder = T, is_nonlinear = T, has_uniform_margins = T)
dim(X$dataset)

g <- compute_gamma_matrix(X$dataset, both_tails = T)

est_dag <- caus_order_to_dag(minimax_search(g))
full_est_dag <- X$dag
full_est_dag[- X$pos_confounders, - X$pos_confounders] <- est_dag

compute_str_int_distance(X$dag, full_est_dag)

est_dag <- lingam_search(X$dataset)
full_est_dag <- X$dag
full_est_dag[- X$pos_confounders, - X$pos_confounders] <- est_dag
compute_str_int_distance(X$dag, full_est_dag)


# Plot graphs with different sparsity parameter
library(igraph)
X <- simulate_data(n, p, prob_connect = 1.5 / (p - 1), distr = "student_t", tail_index = 1.5,
                   has_confounder = T, is_nonlinear = T, has_uniform_margins = T)
g <- graph_from_adjacency_matrix(X$dag)
V(g)$color <- "white"
tkplot(g, layout = layout_in_circle(g))

# Compute theoretical Psi coefficient
library(lattice)
p <- 5
u <- sample(1:1e6, 1)
set.seed(u)
adj_mat <- random_coeff(random_dag(p, .3))
adj_mat
true_psi <- psi_matrix(adj_mat, tail_index = 3.5)

set.seed(u)
sim <- simulate_data(1e5, p, .3, tail_index = 3.5)
est_psi <- causal_tail_matrix(sim$dataset)

path_mat <- get_all_paths(adj_mat)

levelplot(abs(est_psi - true_psi),
          col.regions = heat.colors(100)[length(heat.colors(100)):1])
