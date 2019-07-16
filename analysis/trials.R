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


# Super gaussian
gaus_fam <- function(x, alpha){
  exp(-abs(x) ^ alpha)
}

x <- seq(-2, 2, length.out = 100)
alpha <- 2
plot(gaus_fam(x, alpha))
points(gaus_fam(x, 1.5), col = "blue")
lines(3*dt(x, 1.5))


# Change Lingam
dat <- X$dataset
t.k <- estLiNGAM(dat, only.perm = T, fun = "exp")$k
bb <- prune(t(dat), t.k)
bb$Bpruned


# Test Lingam and PC
X <- simulate_data(1e4, 50, 1.5 / 19, tail_index = 3.5, is_nonlinear = F)
lingam <- lingam_search(X$dataset, "logcosh")
lingam2 <- lingam_search(X$dataset, "exp")
pc <- pc_search(X$dataset, 5e-4)
rank_pc <- pc_rank_search(X$dataset, 5e-4)
greedy <- greedy_ancestral_search(X$dataset)
compute_str_int_distance(X$dag, lingam)
compute_str_int_distance(X$dag, lingam2)
compute_str_int_distance(X$dag, pc)
compute_str_int_distance(X$dag, rank_pc)
compute_str_int_distance(X$dag, caus_order_to_dag(greedy))

# Check simulated_data

for (i in 1:1000){
  set.seed(i)
  X <- simulate_data(1e3, 5, .5, has_confounder = T)
  if (dim(X$dataset)[2] == 0){break()}
}

undebug(simulate_data)
set.seed(i)
X <- simulate_data(1e3, 5, .5, has_confounder = T)


# Check DAG to CPDAG with confounders
set.seed(1991)
for(i in 1:1000){
  n <- sample(2:1e3, 1)
  p <- sample(2:10, 1)
  prob_connect <- runif(1)
  X <- simulate_data(n, p, prob_connect, has_confounder = T)
  if (length(X$pos_confounders) == 0){next()}
  if (length(X$pos_confounders) == 1){
    check <- sum(dag_to_cpdag(X$dag)[X$pos_confounders, ])
  } else {
    check <- apply(dag_to_cpdag(X$dag)[X$pos_confounders, ], 1, sum)
  }
  if(!all(check == 2)){break()}
}


n <- 1e3
p <- 5
prob_connect <- runif(1)
X <- simulate_data(n, p, prob_connect, has_confounder = T)
true_cpdag <- dag_to_cpdag(X$dag)
pc <- causal_discovery(X$dataset, "pc", alpha = 5e-4)
est_ext_cpdag <- X$dag
est_ext_cpdag[-X$pos_confounders, -X$pos_confounders] <- pc$est_cpdag
true_reduced_cpdag <- true_cpdag[-X$pos_confounders, -X$pos_confounders]

SID::hammingDist(true_reduced_cpdag, pc$est_cpdag)

true_dag <- X$dag
est_ext_g <- X$dag
est_ext_g[-X$pos_confounders, -X$pos_confounders] <- pc$est_g

SID::structIntervDist(true_dag, est_ext_g)
SID::structIntervDist(true_dag, dag_to_cpdag(est_ext_g))


est_cpdag <- rbind(c(0, 1, 0),
                   c(1, 0, 0),
                   c(1, 1, 0))
true_dag <- rbind(c(0, 1, 0),
                  c(0, 0, 0),
                  c(1, 1, 0))
SID::structIntervDist(true_dag, dag_to_cpdag(est_cpdag))

true_cpdag <- rbind(c(0, 1, 0),
                    c(0, 0, 1),
                    c(0, 0, 0))
est_cpdag <- rbind(c(0, 1, 1),
                   c(1, 0, 1),
                   c(1, 1, 0))
SID::hammingDist(true_cpdag, est_cpdag)

dag8 <- rbind(c(0, 0, 1, 0, 0, 1),
              c(0, 0, 1, 1, 0, 0),
              c(0, 0, 0, 0, 1, 0),
              rep(0, 6),
              c(rep(0, 5), 1),
              rep(0, 6))
dag_to_cpdag((dag8)[-c(1, 2), -c(1, 2)])
(dag_to_cpdag(dag8))[-c(1, 2), -c(1, 2)]


# Check causal metrics
n <- 1e4
p <- 50
prob_connect <- runif(1)
X <- simulate_data(n, p, 1.5/(p-1))
greedy <- causal_discovery(X$dataset, "greedy")
causal_metrics(X, greedy)
lingam <- causal_discovery(X$dataset, "lingam")
causal_metrics(X, lingam)
pc <- causal_discovery(X$dataset, "pc", alpha = 5e-4)
causal_metrics(X, pc)
pc_rank <- causal_discovery(X$dataset, "pc_rank", alpha = 5e-4)
causal_metrics(X, pc_rank)
undebug(causal_metrics)

# No Confounders
cpdag <- rbind(c(0, 0, 1, 0, 0, 1),
               c(0, 0, 1, 0, 0, 0),
               c(0, 0, 0, 0, 0, 0),
               c(0, 1, 0, 0, 0, 0),
               c(0, 0, 1, 0, 0, 1),
               c(0, 0, 0, 0, 0, 0))

dag <- rbind(c(0, 0, 0, 0, 0, 1),
             c(0, 0, 0, 0, 0, 0),
             c(0, 1, 0, 0, 0, 0),
             c(0, 0, 1, 0, 0, 0),
             c(1, 0, 1, 0, 0, 1),
             c(0, 0, 0, 0, 0, 0))
dag <- rbind(c(0, 0, 0, 0, 0, 1),
             c(0, 0, 0, 1, 0, 0),
             c(1, 1, 0, 0, 1, 0),
             c(0, 0, 0, 0, 0, 0),
             c(0, 0, 0, 0, 0, 1),
             c(0, 0, 0, 0, 0, 0))

dag_to_cpdag(dag) - dag
dag_to_cpdag(dag_to_cpdag(dag_to_cpdag(cpdag))) - dag_to_cpdag(dag_to_cpdag(cpdag))


# Confounders
confounded_true_dag <- rbind(c(0, 1, 0, 0, 0),
                             c(0, 0, 1, 0, 0),
                             c(0, 0, 0, 0, 0),
                             c(1, 0, 1, 0, 0),
                             c(1, 1, 0, 0, 0))

dag_to_cpdag(confounded_true_dag) - confounded_true_dag

est_ext_cpdag <- rbind(c(0, 1, 0, 0, 0),
                       c(1, 0, 1, 0, 0),
                       c(0, 1, 0, 0, 0),
                       c(1, 0, 1, 0, 0),
                       c(1, 1, 0, 0, 0))
dag_to_cpdag(est_ext_cpdag) - est_ext_cpdag
SID::structIntervDist(confounded_true_dag, est_ext_cpdag)
SID::structIntervDist(confounded_true_dag, dag_to_cpdag(est_ext_cpdag))
