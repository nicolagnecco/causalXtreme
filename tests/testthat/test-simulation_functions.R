context("test-simulation_functions")

# Define variables
set.seed(1991)
n <- sample(2:1e3, 1)
p <- sample(2:50, 1)
tail_index <- sample(1:10, 1)
prob_connect <- runif(1)

# Run tests

test_that("simulate data works", {

  # No confounders
  for (distr in c("student_t", "gaussian", "log_normal")){
    out <- simulate_data(n, p, prob_connect, distr, has_confounder = F)

    expect_equal(NCOL(out$dataset), p)
    expect_equal(NCOL(out$dag), p)
    expect_equal(length(out$pos_confounders), 0)
  }

  # With confounders
  for (distr in c("student_t", "gaussian", "log_normal")){
    out <- simulate_data(n, p, prob_connect, distr, has_confounder = T)

    expect_equal(NCOL(out$dataset), p)
    expect_gte(NCOL(out$dag), p)
    expect_gte(length(out$pos_confounders), 0)
  }

  # With non-linear SCM
  for (distr in c("student_t", "gaussian", "log_normal")){
    out <- simulate_data(n, p, prob_connect, distr, is_nonlinear = T)

    expect_equal(NCOL(out$dataset), p)
    expect_gte(NCOL(out$dag), p)
    expect_gte(length(out$pos_confounders), 0)
  }

  # With uniform margins
  for (distr in c("student_t", "gaussian", "log_normal")){
    out <- simulate_data(n, p, prob_connect, distr, has_uniform_margins = T)

    expect_equal(NCOL(out$dataset), p)
    expect_gte(NCOL(out$dag), p)
    expect_gte(length(out$pos_confounders), 0)
  }

  expect_error(simulate_data(1, 2, prob_connect = .4))
  expect_error(simulate_data(2, 1, prob_connect = .4))
  expect_error(simulate_data(1, 1, prob_connect = .4))
})
