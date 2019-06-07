context("test-simulation_helpers")

# Define variables
u <- sample(1e6, 1)

# Run tests
test_that("pick elements works", {

  # vec has length 0
  expect_equal(pick_elements(vec = numeric(0), prob = 0), numeric(0))
  expect_equal(pick_elements(vec = numeric(0), prob = 1), numeric(0))
  expect_equal(pick_elements(vec = numeric(0), prob = 0.3), numeric(0))

  # vec has length 1
  set.seed(u)
  vec <- c(1)
  proba <- 0.3
  r <- rbinom(n = length(vec), size = 1, proba)
  picked <- vec[r == 1]
  expect_equal(pick_elements(vec = vec, prob = 0), numeric(0))
  expect_equal(pick_elements(vec = vec, prob = 1), c(1))
  set.seed(u)
  expect_equal(pick_elements(vec = vec, prob = proba), picked)

  # vec has length > 1
  set.seed(u)
  vec <- c(4, 0, 3, 1)
  proba <- 0.3
  r <- rbinom(n = length(vec), size = 1, proba)
  picked <- vec[r == 1]
  expect_equal(pick_elements(vec = vec, prob = 0), numeric(0))
  expect_equal(pick_elements(vec = vec, prob = 1), vec)
  set.seed(u)
  expect_equal(pick_elements(vec = vec, prob = proba), picked)

})

test_that("inverse mirror uniform works", {
  lb <- runif(1)
  ub <- lb + runif(1)
  expect_equal(inverse_mirror_uniform(prob = 0, min = lb, max = ub), -ub)
  expect_equal(inverse_mirror_uniform(prob = 1, min = lb, max = ub), ub)
  expect_equal(inverse_mirror_uniform(prob = 1 / 2, min = lb, max = ub), lb)
  expect_equal(inverse_mirror_uniform(prob = 1 / 4, min = lb, max = ub),
               2 * 1 / 4 * (ub - lb) - ub)
  expect_equal(inverse_mirror_uniform(prob = 3 / 4, min = lb, max = ub),
              (2 * 3 / 4 - 1) * (ub - lb) + lb)

  expect_error(inverse_mirror_uniform(prob = runif(1), min = ub, max = lb))
  expect_error(inverse_mirror_uniform(prob = runif(1), min = -ub, max = lb))
  expect_error(inverse_mirror_uniform(prob = runif(1), min = ub, max = -lb))
  expect_error(inverse_mirror_uniform(prob = runif(1), min = -ub, max = -lb))

})

test_that("sample uniform works", {
  n <- 7
  lb <- runif(1)
  ub <- lb + runif(1)

  # mirror == FALSE
  set.seed(u)
  r <- runif(n, min = lb, max = ub)
  set.seed(u)
  expect_equal(sample_uniform(n = n, min = lb, max = ub), r)

  set.seed(u)
  r <- runif(n, min = -ub, max = -lb)
  set.seed(u)
  expect_equal(sample_uniform(n = n, min = -ub, max = -lb), r)

  expect_error(sample_uniform(n = n, min = ub, max = lb))

  # mirror == TRUE
  set.seed(u)
  r <- sapply(runif(n), inverse_mirror_uniform, min = lb, max = ub)
  set.seed(u)
  expect_equal(sample_uniform(n = n, min = lb, max = ub, mirror = TRUE), r)

  expect_error(sample_uniform(n = n, min = -ub, max = -lb, mirror = TRUE))
  expect_error(sample_uniform(n = n, min = -ub, max =  lb, mirror = TRUE))
  expect_error(sample_uniform(n = n, min =  ub, max = -lb, mirror = TRUE))
  expect_error(sample_uniform(n = n, min = ub, max = lb, mirror = TRUE))

})

test_that("random dag works", {
  # caus_order has length 0
  expect_error(random_dag(p = 0, prob = 0, caus_order = numeric(0)))
  expect_error(random_dag(p = 0, prob = 1, caus_order = numeric(0)))
  expect_error(random_dag(p = 0, prob = 0.4, caus_order = numeric(0)))

  # caus_order has length 1
  caus_order <- c(1)

  expect_equal(random_dag(p = 1, prob = 0, caus_order = caus_order),
               matrix(0))
  expect_equal(random_dag(p = 1, prob = 1, caus_order = caus_order),
               matrix(0))
  expect_equal(random_dag(p = 1, prob = 0.5, caus_order = caus_order),
               matrix(0))

  # caus_order has length > 1
  p <- 4
  caus_order <- sample(p)

  dag_empty <- matrix(0, nrow = length(caus_order), ncol = length(caus_order))
  expect_equal(random_dag(p = p, prob = 0, caus_order = caus_order),
               dag_empty)

  dag_full <-  caus_order_to_adjmat(caus_order)
  expect_equal(random_dag(p = p, prob = 1, caus_order = caus_order),
               dag_full)

  proba <- runif(1)
  dag <- matrix(0, nrow = length(caus_order), ncol = length(caus_order))

  set.seed(u)
  elms <- pick_elements(caus_order[2:p], proba)
  dag[caus_order[1], elms] <- 1
  elms <- pick_elements(caus_order[3:p], proba)
  dag[caus_order[2], elms] <- 1
  elms <- pick_elements(caus_order[4:p], proba)
  dag[caus_order[3], elms] <- 1

  set.seed(u)
  expect_equal(random_dag(p = p, prob = proba, caus_order = caus_order),
               dag)


  # mismatch between length of caus_order and p
  expect_error(random_dag(p = 1, prob = 0, caus_order = sample(p)))
  expect_error(random_dag(p = 1, prob = 1, caus_order = sample(p)))
  expect_error(random_dag(p = 1, prob = 0.4, caus_order = sample(p)))



})

test_that("random coefficients works", {
  lb <- runif(1)
  ub <- lb + runif(1)


  # DAG has one variable
  expect_equal(random_coeff(matrix(0), lb = lb, ub = ub,
                            two_intervals = FALSE), matrix(0))
  expect_equal(random_coeff(matrix(0), lb = lb, ub = ub,
                            two_intervals = TRUE), matrix(0))

  # DAG has > 1 variable
  p <- 7
  g <- random_dag(p = p, prob = runif(1))

  num_coeff <- sum(g)

  # two_intervals == TRUE
  g_coeff <- matrix(0, nrow = p, ncol = p)
  set.seed(u)
  g_coeff[g == 1] <- sample_uniform(n = num_coeff, min = lb, max = ub,
                                           mirror = TRUE)
  set.seed(u)
  expect_equal(random_coeff(g, lb = lb, ub = ub, two_intervals = TRUE),
               g_coeff)

  # two_intervals == FALSE
  g_coeff <- matrix(0, nrow = p, ncol = p)
  set.seed(u)
  g_coeff[g == 1] <- sample_uniform(n = num_coeff, min = lb, max = ub,
                                    mirror = FALSE)
  set.seed(u)
  expect_equal(random_coeff(g, lb = lb, ub = ub, two_intervals = FALSE),
               g_coeff)

  # not adjacency matrix
  expect_error(random_coeff(g * 0.1, lb = lb, ub = ub, two_intervals = FALSE))
  expect_error(random_coeff(g * 0.1, lb = lb, ub = ub, two_intervals = TRUE))

  # lb > ub
  expect_error(random_coeff(g, lb = ub, ub = lb, two_intervals = FALSE))
  expect_error(random_coeff(g, lb = ub, ub = lb, two_intervals = TRUE))

  # if two_intervals == TRUE, lb > 0 and ub > 0
  expect_error(random_coeff(g, lb = -ub, ub =  lb, two_intervals = TRUE))
  expect_error(random_coeff(g, lb =  ub, ub = -lb, two_intervals = TRUE))
  expect_error(random_coeff(g, lb = -ub, ub = -lb, two_intervals = TRUE))

})
