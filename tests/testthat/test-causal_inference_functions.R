context("test-causal_inference_functions")

# Define variables
n <- sample(2:1e3, 1)
p <- sample(2:10, 1)
prob_connect <- runif(1)
X <- simulate_data(n, p, prob_connect)
dat <- X$dataset
dat2 <- matrix("bo", nrow = n, ncol = p)

true_dag <- rbind(c(0, 0, 1, 0, 0, 1),
                  c(0, 0, 1, 1, 0, 0),
                  c(0, 0, 0, 0, 1, 0),
                  c(0, 0, 0, 0, 0, 0),
                  c(0, 0, 0, 0, 0, 1),
                  c(0, 0, 0, 0, 0, 0))
true_cpdag <- rbind(c(0, 0, 1, 0, 0, 1),
                    c(0, 0, 1, 1, 0, 0),
                    c(0, 0, 0, 0, 1, 0),
                    c(0, 1, 0, 0, 0, 0),
                    c(0, 0, 0, 0, 0, 1),
                    c(0, 0, 0, 0, 0, 0))
est_dag1 <- rbind(c(0, 0, 0, 0, 0, 1),
                  c(0, 0, 0, 0, 0, 0),
                  c(0, 1, 0, 0, 0, 0),
                  c(0, 0, 1, 0, 0, 0),
                  c(1, 0, 1, 0, 0, 1),
                  c(0, 0, 0, 0, 0, 0))
est_cpdag1 <- rbind(c(0, 0, 0, 0, 1, 1),
                    c(0, 0, 0, 0, 0, 0),
                    c(0, 1, 0, 0, 0, 0),
                    c(0, 0, 1, 0, 0, 0),
                    c(1, 0, 1, 0, 0, 1),
                    c(1, 0, 0, 0, 1, 0))
est_cpdag2 <-rbind(c(0, 0, 1, 0, 0, 1),
                   c(0, 0, 1, 1, 0, 0),
                   c(1, 1, 0, 0, 1, 0),
                   c(0, 1, 0, 0, 0, 0),
                   c(0, 0, 1, 0, 0, 1),
                   c(0, 0, 0, 0, 0, 0))

# Run tests
test_that("causal discovery works", {
  expect_error(causal_discovery(dat, "greedy", foobar = 20, both_tails = T))
  expect_error(causal_discovery(dat, "greedy", k = 20, foobar = T))
  expect_error(causal_discovery(dat, "greedy", foobar = 20, foobar = T))
  expect_length(causal_discovery(dat, "greedy", k = 20, both_tails = T), 2)
  expect_length(causal_discovery(dat, "greedy", k = 20), 2)
  expect_length(causal_discovery(dat, "greedy", both_tails = T), 2)
  expect_length(causal_discovery(dat, "greedy"), 2)
  expect_error(causal_discovery(dat, "lingam", foobar = "exp"))
  expect_length(causal_discovery(dat, "lingam", contrast_fun = "exp"), 2)
  expect_length(causal_discovery(dat, "lingam"), 2)
  expect_length(causal_discovery(dat2, "lingam"), 2)
  expect_error(causal_discovery(dat, "pc", foobar = 5e-3))
  expect_length(causal_discovery(dat, "pc", alpha = 5e-3), 2)
  expect_length(causal_discovery(dat, "pc"), 2)
  expect_length(causal_discovery(dat2, "pc"), 2)
  expect_error(causal_discovery(dat, "pc_rank", foobar = 5e-3))
  expect_length(causal_discovery(dat, "pc_rank", alpha = 5e-3), 2)
  expect_length(causal_discovery(dat, "pc_rank"), 2)
  expect_length(causal_discovery(dat2, "pc_rank"), 2)
  expect_error(causal_discovery(dat, "random", foobar = 5e-3))
  expect_length(causal_discovery(dat, "random"), 2)
  expect_error(causal_discovery(dat))
  expect_error(causal_discovery(dat, "foobar"))
})

test_that("causal metrics work",{
  # No confounders
  # DAG and CPDAG
  sim_data <- list(NA,
                   dag = true_dag,
                   pos_confounders = integer(0))
  est_graphs <- list(est_g = est_dag1,
                     est_cpdag = est_cpdag1)

  out <- list(sid = compute_str_int_distance(true_dag, est_dag1),
          shd = compute_str_ham_distance(true_cpdag, est_cpdag1))
  expect_equal(causal_metrics(sim_data, est_graphs), out)

  # Both CPDAGs
  sim_data <- list(NA,
                   dag = true_dag,
                   pos_confounders = integer(0))
  est_graphs <- list(est_g = est_cpdag2,
                     est_cpdag = est_cpdag2)
  out <- list(sid = compute_str_int_distance(true_dag, est_cpdag2),
              shd = compute_str_ham_distance(true_cpdag, est_cpdag2))
  expect_equal(causal_metrics(sim_data, est_graphs), out)

  # Confounders
  # !!! continue
})
