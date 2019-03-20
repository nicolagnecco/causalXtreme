context("test-causal_extreme_functions")

# Define variables
v1     <- c(0, 2, 6, 1, 9, 8, 7, 5, 3, 4)
v1_rep <- c(0, 2, 2, 1, 2, 8, 7, 5, 3, 4)
v2     <- c(5, 6, 4, 2, 0, 8, 9, 1, 7, 3)
v2_rep <- c(5, 6, 4, 2, 0, 8, 8, 8, 7, 3)
mat <- cbind(v1, v1_rep, v2)

A <- rbind(c(0, 1.5, 0, 0),
           c(0, 0, 1.2, 2),
           c(0, 0, 0, 0.3),
           c(0, 0, 0, 0))
alpha <- 1.5
noise_w <- c(1, 0.5 ^ (1 / alpha),
             0.5 ^ (1 / alpha), 0.5 ^ (1 / alpha))

# Run tests
test_that("gamma coefficient works", {
  expect_equal(compute_gamma(v2, v1), (4 + 9 + 8) / 11 * 1 / 3)
  expect_equal(compute_gamma(v2, v1, k = 2), (8 + 9) / 11 * 1 / 2)
  expect_equal(compute_gamma(v1, v2), (10 + 9 + 1) / 11 * 1 / 3)
  expect_equal(compute_gamma(v2_rep, v1), (9 + 8 + 6) / 11 * 1 / 3)
  expect_equal(compute_gamma(v2_rep, v1, k = 4),
               (9 + 8 + 6 + 4) / 11 * 1 / 4)
  expect_equal(compute_gamma(v1_rep, v2_rep),
               (8 + 9 + 10) / 11 * 1 / 3)
  expect_equal(compute_gamma(v2_rep, v1_rep, k = 4),
               (6 + 8 + 9 + 10) / 11 * 1 / 4)
})

test_that("gamma matrix works", {
  expect_equal(compute_gamma_matrix(mat),
               rbind(c(NA,
                       compute_gamma(v1, v1_rep),
                       compute_gamma(v1, v2)),
                     c(compute_gamma(v1_rep, v1),
                       NA,
                       compute_gamma(v1_rep, v2)),
                     c(compute_gamma(v2, v1),
                       compute_gamma(v2, v1_rep),
                       NA)))
  expect_equal(compute_gamma_matrix(mat, k = 4),
               rbind(c(NA,
                       compute_gamma(v1, v1_rep, k = 4),
                       compute_gamma(v1, v2, k = 4)),
                     c(compute_gamma(v1_rep, v1, k = 4),
                       NA,
                       compute_gamma(v1_rep, v2, k = 4)),
                     c(compute_gamma(v2, v1, k = 4),
                       compute_gamma(v2, v1_rep, k = 4),
                       NA)))
})

test_that("theoretical gamma works", {
  gamma_21 <- 1 - 0.5 *
    ( (noise_w[2] ^ alpha) /
       (noise_w[2] ^ alpha + noise_w[1] ^ alpha * A[1, 2] ^ alpha))

  gamma_31 <- 1 - 0.5 *
    ( (noise_w[3] ^ alpha + noise_w[2] ^ alpha * A[2, 3] ^ alpha) /
       (noise_w[3] ^ alpha +
          noise_w[2] ^ alpha * A[2, 3] ^ alpha +
          noise_w[1] ^ alpha * A[1, 2] ^ alpha * A[2, 3] ^ alpha))

  gamma_32 <- 1 - 0.5 *
    ( (noise_w[3] ^ alpha) /
       (noise_w[3] ^ alpha +
          noise_w[2] ^ alpha * A[2, 3] ^ alpha +
          noise_w[1] ^ alpha * A[1, 2] ^ alpha * A[2, 3] ^ alpha))

  gamma_41 <- 1 - 0.5 *
    ( (noise_w[4] ^ alpha +
        noise_w[3] ^ alpha * A[3, 4] ^ alpha +
        noise_w[2] ^ alpha * (A[2, 4] + A[2, 3] * A[3, 4]) ^ alpha) /
       (noise_w[4] ^ alpha +
          noise_w[3] ^ alpha * A[3, 4] ^ alpha +
          noise_w[2] ^ alpha * (A[2, 4] + A[2, 3] * A[3, 4]) ^ alpha +
          noise_w[1] ^ alpha * ( (A[2, 4] + A[2, 3] * A[3, 4]) *
                                  A[1, 2]) ^ alpha))

  gamma_42 <- 1 - 0.5 *
    ( (noise_w[4] ^ alpha +
        noise_w[3] ^ alpha * A[3, 4] ^ alpha) /
       (noise_w[4] ^ alpha +
          noise_w[3] ^ alpha * A[3, 4] ^ alpha +
          noise_w[2] ^ alpha * (A[2, 4] + A[2, 3] * A[3, 4]) ^ alpha +
          noise_w[1] ^ alpha * ( (A[2, 4] + A[2, 3] * A[3, 4]) *
                                  A[1, 2]) ^ alpha))

  gamma_43 <- 1 - 0.5 *
    ( (noise_w[4] ^ alpha) /
       (noise_w[4] ^ alpha +
          noise_w[3] ^ alpha * A[3, 4] ^ alpha +
          noise_w[2] ^ alpha * (A[2, 4] + A[2, 3] * A[3, 4]) ^ alpha +
          noise_w[1] ^ alpha * ( (A[2, 4] + A[2, 3] * A[3, 4]) *
                                  A[1, 2]) ^ alpha))

  expect_equal(compute_gamma_theo(A, 1, 2, alpha, noise_w),
               1)
  expect_equal(compute_gamma_theo(A, 1, 3, alpha, noise_w),
               1)
  expect_equal(compute_gamma_theo(A, 1, 4, alpha, noise_w),
               1)
  expect_equal(compute_gamma_theo(A, 2, 3, alpha, noise_w),
               1)
  expect_equal(compute_gamma_theo(A, 2, 4, alpha, noise_w),
               1)
  expect_equal(compute_gamma_theo(A, 3, 4, alpha, noise_w),
               1)
  expect_equal(compute_gamma_theo(A, 4, 3, alpha, noise_w),
               gamma_43)
  expect_equal(compute_gamma_theo(A, 4, 2, alpha, noise_w),
               gamma_42)
  expect_equal(compute_gamma_theo(A, 4, 1, alpha, noise_w),
               gamma_41)
  expect_equal(compute_gamma_theo(A, 3, 2, alpha, noise_w),
               gamma_32)
  expect_equal(compute_gamma_theo(A, 3, 1, alpha, noise_w),
               gamma_31)
  expect_equal(compute_gamma_theo(A, 2, 1, alpha, noise_w),
               gamma_21)
  expect_error(compute_gamma_theo(-A, 2, 1, alpha, noise_w))
})

test_that("theoretical gamma matrix works", {
  expect_equal(compute_gamma_theo_matrix(A, alpha, noise_w),
               rbind(c(NA,
                       compute_gamma_theo(A, 1, 2, alpha, noise_w),
                       compute_gamma_theo(A, 1, 3, alpha, noise_w),
                       compute_gamma_theo(A, 1, 4, alpha, noise_w)),
                     c(compute_gamma_theo(A, 2, 1, alpha, noise_w),
                       NA,
                       compute_gamma_theo(A, 2, 3, alpha, noise_w),
                       compute_gamma_theo(A, 2, 4, alpha, noise_w)),
                     c(compute_gamma_theo(A, 3, 1, alpha, noise_w),
                       compute_gamma_theo(A, 3, 2, alpha, noise_w),
                       NA,
                       compute_gamma_theo(A, 3, 4, alpha, noise_w)),
                     c(compute_gamma_theo(A, 4, 1, alpha, noise_w),
                       compute_gamma_theo(A, 4, 2, alpha, noise_w),
                       compute_gamma_theo(A, 4, 3, alpha, noise_w),
                       NA)))
  expect_error(compute_gamma_theo_matrix(-A, alpha, noise_w))
})
