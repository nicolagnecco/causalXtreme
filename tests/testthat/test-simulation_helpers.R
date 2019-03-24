context("test-simulation_helpers")

# Define variables
u <- sample(1e6, 1)
n_iter <- 2

g1 <- rbind(c(0, 0, 1, 0),
            c(1, 0, 0, 0),
            c(0, 0, 0, 1),
            c(0, 0, 0, 0))

g2 <- rbind(c(0, 1.5, 0, 0),
            c(0, 0, -1.2, 2),
            c(0, 0, 0, 0.3),
            c(0, 0, 0, 0))

g3 <- rbind(c(0, 1, 0),
            c(0, 0, 0),
            c(0, 0, 0))

lb_vec <- runif(10)

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




test_that("random dag works", {
  # prob_connect = 1, caus_order = given
  for (i in 1:n_iter){
    for (p in 3:7){
      caus_order <- sample(p, p, replace = FALSE)
      dag <- caus_order_to_adjmat(caus_order)
      dag_sparse <- Matrix::Matrix(nrow = p, ncol = p, 0, sparse = TRUE)
      dag_sparse[dag == 1] <- 1

      expect_equal(random_dag(p, 1, caus_order),
                   dag)
      expect_equal(random_dag(p, 1, caus_order, sparse = TRUE),
                   dag_sparse)
    }
  }

  # prob_connect = 1, caus_order = not given
  for (i in 1:n_iter){
    for (p in 3:7){
      u <- sample(1e4, 1)
      set.seed(u)
      caus_order <- sample(p, p, replace = FALSE)[p:1]
      dag <- caus_order_to_adjmat(caus_order)
      dag_sparse <- Matrix::Matrix(nrow = p, ncol = p, 0, sparse = TRUE)
      dag_sparse[dag == 1] <- 1

      set.seed(u)
      expect_equal(random_dag(p, 1), dag)
      set.seed(u)
      expect_equal(random_dag(p, 1, sparse = TRUE),
                   dag_sparse)
    }
  }

  # prob_connect = 0, caus_order = given
  for (i in 1:n_iter){
    for (p in 3:7){
      caus_order <- sample(p, p, replace = FALSE)
      dag <- matrix(0, ncol = p, nrow = p)
      dag_sparse <- Matrix::Matrix(nrow = p, ncol = p, 0, sparse = TRUE)

      expect_equal(random_dag(p, 0, caus_order),
                   dag)
      expect_equal(random_dag(p, 0, caus_order, sparse = TRUE),
                   dag_sparse)
    }
  }

  # prob_connect = 0, caus_order = not given
  for (i in 1:n_iter){
    for (p in 3:7){
      dag <- matrix(0, ncol = p, nrow = p)
      dag_sparse <- Matrix::Matrix(nrow = p, ncol = p, 0, sparse = TRUE)

      expect_equal(random_dag(p, 0), dag)
      expect_equal(random_dag(p, 0, sparse = TRUE),
                   dag_sparse)
    }
  }

  # 0 < prob_connect < 1, caus_order = given
  for (i in 1:n_iter){
    for (p in 3:7){
      caus_order <- sample(p, p, replace = FALSE)
      dag <- caus_order_to_adjmat(caus_order)
      dag_sparse <- Matrix::Matrix(nrow = p, ncol = p, 0, sparse = TRUE)
      dag_sparse[dag == 1] <- 1

      expect_true(all(dag - random_dag(p, runif(1), caus_order) >= 0))
      expect_true(all(dag_sparse -
                        random_dag(p, runif(1),
                                   caus_order, sparse = TRUE) >= 0))
    }
  }

  # 0 < prob_connect < 1, caus_order = not given
  for (i in 1:n_iter){
    for (p in 3:7){
      u <- sample(1e4, 1)
      set.seed(u)
      caus_order <- sample(p, p, replace = FALSE)[p:1]
      dag <- caus_order_to_adjmat(caus_order)
      dag_sparse <- Matrix::Matrix(nrow = p, ncol = p, 0, sparse = TRUE)
      dag_sparse[dag == 1] <- 1

      proba <- runif(1)
      set.seed(u)
      expect_true(all(dag - random_dag(p, proba) >= 0))
      proba <- runif(1)
      set.seed(u)
      expect_true(all(dag_sparse -
                        random_dag(p, proba, caus_order, sparse = TRUE) >= 0))
    }
  }

  # errors
  for (p in 1:2){
    expect_error(random_dag(p, 0.3))
    expect_error(random_dag(p, 0.4, caus_order = 1:4))
    expect_error(random_dag(p, 0.4, caus_order = 1:4, sparse = TRUE))
    expect_error(random_dag(p, 0.4, sparse = TRUE))
  }

  for (p in 3:5){
    expect_error(random_dag(p, 0.4, caus_order = 1:7))
    expect_error(random_dag(p, 0.4, caus_order = 1:7, sparse = TRUE))
  }

})

test_that("random weighted dag works", {

  # when num_coeff > 1
  for (lb in lb_vec){
    for (two_intervals in c(T, F)){

      ub <- lb + runif(1)

      set.seed(1)
      num_coeff <- sum(g1)
      g <- matrix(0, nrow = NROW(g1), ncol = NCOL(g1))
      rand_coeff <- sample(c(-1, 1), num_coeff, replace = TRUE) ^
        (two_intervals) * runif(num_coeff, min = lb, max = ub)
      g[g1 != 0] <- rand_coeff

      set.seed(1)
      expect_equal(random_b(g1, lb, ub, two_intervals), g)
    }
  }

  # when num_coeff == 1
  for (lb in lb_vec){
    for (two_intervals in c(T, F)){

      ub <- lb + runif(1)

      set.seed(1)
      num_coeff <- sum(g3)
      g <- matrix(0, nrow = NROW(g3), ncol = NCOL(g3))
      rand_coeff <- sample(c(-1, 1), num_coeff, replace = TRUE) ^
        (two_intervals) * runif(num_coeff, min = lb, max = ub)
      g[g3 != 0] <- rand_coeff

      set.seed(1)
      expect_equal(random_b(g3, lb, ub, two_intervals), g)
    }
  }

  # error when g is not adjacency matrix
  lb <- runif(1)
  ub <- lb + runif(1)
  expect_error(random_b(g2, lb, ub, FALSE))
  expect_error(random_b(g2, lb, ub, TRUE))

  # error when lb > ub
  lb <- runif(1)
  ub <- max(lb - runif(1), 0)
  expect_error(random_b(g1, lb, ub, FALSE))
  expect_error(random_b(g1, lb, ub, TRUE))

  # error when two_intervals == TRUE and lb < 0 and ub < 0
  expect_error(random_b(g1, -0.3, -0.1, TRUE))
  expect_error(random_b(g1, -0.3,  0.1, TRUE))
  expect_error(random_b(g1,  0.3, -0.1, TRUE))

})
