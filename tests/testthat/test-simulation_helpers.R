context("test-simulation_helpers")

# Define variables


# Run tests
test_that("random dag works", {
  # prob_connect = 1, caus_order = given
  for (i in 1:3){
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
  for (i in 1:3){
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
  for (i in 1:3){
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
  for (i in 1:3){
    for (p in 3:7){
      dag <- matrix(0, ncol = p, nrow = p)
      dag_sparse <- Matrix::Matrix(nrow = p, ncol = p, 0, sparse = TRUE)

      expect_equal(random_dag(p, 0), dag)
      expect_equal(random_dag(p, 0, sparse = TRUE),
                   dag_sparse)
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

test_that("random coefficients works", {
  # !!! write tests for the function
})
