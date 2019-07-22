context("test-causal_inference_functions")

# Define variables
set.seed(1991)
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
est_cpdag2 <- rbind(c(0, 0, 1, 0, 0, 1),
                   c(0, 0, 1, 1, 0, 0),
                   c(1, 1, 0, 0, 1, 0),
                   c(0, 1, 0, 0, 0, 0),
                   c(0, 0, 1, 0, 0, 1),
                   c(0, 0, 0, 0, 0, 0))
true_dag3 <- rbind(c(0, 1, 0, 0, 0),
                   c(0, 0, 1, 0, 0),
                   c(0, 0, 0, 0, 0),
                   c(1, 0, 1, 0, 0),
                   c(1, 1, 0, 0, 0))
true_red_cpdag3 <- rbind(c(0, 1, 0),
                         c(1, 0, 1),
                         c(0, 1, 0))
est_dag4 <- rbind(c(0, 1, 1),
                  c(0, 0, 0),
                  c(0, 1, 0))
est_ext_dag4 <- rbind(c(0, 1, 1, 0, 0),
                      c(0, 0, 0, 0, 0),
                      c(0, 1, 0, 0, 0),
                      c(1, 0, 1, 0, 0),
                      c(1, 1, 0, 0, 0))
est_cpdag4 <- rbind(c(0, 1, 1),
                    c(1, 0, 1),
                    c(1, 1, 0))
est_cpdag5 <- true_red_cpdag3
est_ext_cpdag5 <- rbind(c(0, 1, 1, 0, 0),
                        c(1, 0, 1, 0, 0),
                        c(1, 1, 0, 0, 0),
                        c(1, 0, 1, 0, 0),
                        c(1, 1, 0, 0, 0))

# Run tests
test_that("causal discovery works", {
  expect_error(causal_discovery(dat, "ease", list(foobar = 20, both_tails = T)))
  expect_error(causal_discovery(dat, "ease", list(k = 20, foobar = T)))
  expect_error(causal_discovery(dat, "ease", list(foobar = 20, foobar = T)))
  expect_error(causal_discovery(dat, "ease",
                                list(k = 20, both_tails = T, foobar = 3)))
  expect_length(causal_discovery(dat, "ease", list(k = 20, both_tails = T)), 2)
  expect_length(causal_discovery(dat, "ease", list(k = 20)), 2)
  expect_length(causal_discovery(dat, "ease", list(both_tails = T)), 2)
  expect_length(causal_discovery(dat, "ease"), 2)
  expect_error(causal_discovery(dat, "lingam", list(foobar = "exp")))
  expect_error(causal_discovery(dat, "lingam",
                                list(contrast_fun = "logcosh", foobar = "exp")))
  expect_length(causal_discovery(dat, "lingam", list(contrast_fun = "exp")), 2)
  expect_length(causal_discovery(dat, "lingam"), 2)
  expect_length(causal_discovery(dat2, "lingam"), 2)
  expect_error(causal_discovery(dat, "pc", list(foobar = 5e-3)))
  expect_error(causal_discovery(dat, "pc", list(alpha = 5e-3, foobar = 2)))
  expect_length(causal_discovery(dat, "pc", list(alpha = 5e-3)), 2)
  expect_length(causal_discovery(dat, "pc"), 2)
  expect_length(causal_discovery(dat2, "pc"), 2)
  expect_error(causal_discovery(dat, "pc_rank", list(foobar = 5e-3)))
  expect_error(causal_discovery(dat, "pc_rank", list(alpha = 5e-3, foobar = 2)))
  expect_length(causal_discovery(dat, "pc_rank", list(alpha = 5e-3)), 2)
  expect_length(causal_discovery(dat, "pc_rank"), 2)
  expect_length(causal_discovery(dat2, "pc_rank"), 2)
  expect_error(causal_discovery(dat, "random", list(foobar = 5e-3)))
  expect_length(causal_discovery(dat, "random"), 2)
  expect_error(causal_discovery(dat))
  expect_error(causal_discovery(dat, "foobar"))
})

test_that("causal metrics work", {
  ## No confounders
  # (Output from EASE, Lingam, Random) DAG and CPDAG
  sim_data <- list(NA,
                   dag = true_dag,
                   pos_confounders = integer(0))
  est_graphs <- list(est_g = est_dag1,
                     est_cpdag = est_cpdag1)

  out <- list(sid = compute_str_int_distance(true_dag, est_dag1),
          shd = compute_str_ham_distance(true_cpdag, est_cpdag1))
  expect_equal(causal_metrics(sim_data, est_graphs), out)

  # (Output from PC and PC Rank) Both CPDAGs
  sim_data <- list(NA,
                   dag = true_dag,
                   pos_confounders = integer(0))
  est_graphs <- list(est_g = est_cpdag2,
                     est_cpdag = est_cpdag2)
  out <- list(sid = compute_str_int_distance(true_dag, est_cpdag2),
              shd = compute_str_ham_distance(true_cpdag, est_cpdag2))
  expect_equal(causal_metrics(sim_data, est_graphs), out)

  ## Confounders
  # (Output from EASE, Lingam, Random) DAG and CPDAG
  sim_data <- list(NA,
                   dag = true_dag3,
                   pos_confounders = c(4, 5))
  est_graphs <- list(est_g = est_dag4,
                     est_cpdag = est_cpdag4)

  out <- list(sid = compute_str_int_distance(true_dag3, est_ext_dag4),
              shd = compute_str_ham_distance(true_red_cpdag3, est_cpdag4))
  expect_equal(causal_metrics(sim_data, est_graphs), out)

  # (Output from PC and PC Rank) Both CPDAGs
  sim_data <- list(NA,
                   dag = true_dag3,
                   pos_confounders = c(4, 5))
  est_graphs <- list(est_g = est_cpdag5,
                     est_cpdag = est_cpdag5)
  out <- list(sid = compute_str_int_distance(true_dag3, est_ext_cpdag5),
              shd = compute_str_ham_distance(true_red_cpdag3, est_cpdag5))
  expect_equal(causal_metrics(sim_data, est_graphs), out)

  ## Estimated graphs are NA
  sim_data <- list(NA,
                   dag = true_dag3,
                   pos_confounders = c(4, 5))
  est_graphs <- list(est_g = NA,
                     est_cpdag = est_cpdag5)
  out <- list(sid = NA,
              shd = NA)
  expect_equal(causal_metrics(sim_data, est_graphs), out)

  sim_data <- list(NA,
                   dag = true_dag3,
                   pos_confounders = c(4, 5))
  est_graphs <- list(est_g = est_cpdag5,
                     est_cpdag = NA)
  out <- list(sid = NA,
              shd = NA)
  expect_equal(causal_metrics(sim_data, est_graphs), out)

  sim_data <- list(NA,
                   dag = true_dag3,
                   pos_confounders = c(4, 5))
  est_graphs <- list(est_g = NA,
                     est_cpdag = NA)
  out <- list(sid = NA,
              shd = NA)
  expect_equal(causal_metrics(sim_data, est_graphs), out)
})
