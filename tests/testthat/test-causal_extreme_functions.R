context("test-causal_extreme_functions")

# Define variables
v1 <- c(0, 2, 6, 1, 9, 8, 7, 5, 3, 4)
r1 <- c(1, 3, 7, 2, 10, 9, 8, 6, 4, 5)
v1_rep <- c(0, 2, 2, 1, 2, 8, 7, 5, 3, 4)
r1_rep <- c(1, 3, 4, 2, 5, 10, 9, 8, 6, 7)
v2 <- c(5, 6, 4, 2, 0, 8, 9, 1, 7, 3)
r2 <- c(6, 7, 5, 3, 1, 9, 10, 2, 8, 4)
v2_rep <- c(5, 6, 4, 2, 0, 8, 8, 8, 7, 3)
r2_rep <- c(5, 6, 4, 2, 1, 8, 9, 10, 7, 3)

mat <- cbind(v1, v1_rep, v2)

A <- rbind(c(0, 1.5, 0, 0),
           c(0, 0, 1.2, 2),
           c(0, 0, 0, 0.3),
           c(0, 0, 0, 0))
alpha <- 1.5
noise_w <- c(1, 0.5 ^ (1 / alpha),
             0.5 ^ (1 / alpha), 0.5 ^ (1 / alpha))

# Run tests
test_that("causal tail coefficient works", {
  expect_equal(causal_tail_coeff(v2, v1, k = 3, both_tails = FALSE),
               (4 + 9 + 8) / 10 * 1 / 3)
  expect_equal(causal_tail_coeff(v2, v1, k = 2, both_tails = FALSE),
               (8 + 9) / 10 * 1 / 2)
  expect_equal(causal_tail_coeff(v1, v2, k = 3, both_tails = FALSE),
               (10 + 9 + 1) / 10 * 1 / 3)
  expect_equal(causal_tail_coeff(v2_rep, v1, k = 3, both_tails = FALSE),
               (9 + 8 + 6) / 10 * 1 / 3)
  expect_equal(causal_tail_coeff(v2_rep, v1, k = 4, both_tails = FALSE),
               (9 + 8 + 6 + 4) / 10 * 1 / 4)
  expect_equal(causal_tail_coeff(v1_rep, v2_rep, k = 3, both_tails = FALSE),
               (8 + 9 + 10) / 10 * 1 / 3)
  expect_equal(causal_tail_coeff(v2_rep, v1_rep, k = 4, both_tails = FALSE),
               (6 + 8 + 9 + 10) / 10 * 1 / 4)
  expect_equal(causal_tail_coeff(r2, r1, k = 2, to_rank = FALSE,
                                 both_tails = FALSE),
               (8 + 9) / 10 * 1 / 2)
  expect_equal(causal_tail_coeff(r1, r2, to_rank = FALSE, k = 3,
                                 both_tails = FALSE),
               (10 + 9 + 1) / 10 * 1 / 3)
  expect_equal(causal_tail_coeff(v1_rep, v2_rep, k = 3, both_tails = TRUE),
               1 / (2 * 10) * 2 * (abs(5 - 5.5) + abs(8 - 5.5)))
  expect_equal(causal_tail_coeff(v2, v1, k = 5, both_tails = TRUE),
               1 / (4 * 10) * 2 * (abs(10 - 5.5) + abs(9 - 5.5) +
                                     abs(8 - 5.5) + abs(6 - 5.5)))
  expect_error(causal_tail_coeff(v2, v1, k = 1))
  expect_error(causal_tail_coeff(v1, v2, k = 10))
})

test_that("causal tail matrix works", {
  expect_equal(causal_tail_matrix(mat, both_tails = FALSE),
               rbind(c(NA,
                       causal_tail_coeff(v1, v1_rep, both_tails = FALSE),
                       causal_tail_coeff(v1, v2, both_tails = FALSE)),
                     c(causal_tail_coeff(v1_rep, v1, both_tails = FALSE),
                       NA,
                       causal_tail_coeff(v1_rep, v2, both_tails = FALSE)),
                     c(causal_tail_coeff(v2, v1, both_tails = FALSE),
                       causal_tail_coeff(v2, v1_rep, both_tails = FALSE),
                       NA)))
  expect_equal(causal_tail_matrix(mat, k = 4, both_tails = TRUE),
               rbind(c(NA,
                       causal_tail_coeff(v1, v1_rep, k = 4, both_tails = TRUE),
                       causal_tail_coeff(v1, v2, k = 4, both_tails = TRUE)),
                     c(causal_tail_coeff(v1_rep, v1, k = 4, both_tails = TRUE),
                       NA,
                       causal_tail_coeff(v1_rep, v2, k = 4, both_tails = TRUE)),
                     c(causal_tail_coeff(v2, v1, k = 4, both_tails = TRUE),
                       causal_tail_coeff(v2, v1_rep, k = 4, both_tails = TRUE),
                       NA)))
})

test_that("psi coefficient works", {
  psi_21 <- 1 - 0.5 *
    ( (1) /
        (1 + 1 * abs(A[1, 2]) ^ alpha))

  psi_31 <- 1 - 0.5 *
    ( (1 + 1 * abs(A[2, 3]) ^ alpha) /
        (1 +
           1 * abs(A[2, 3]) ^ alpha +
           1 * abs(A[1, 2]) ^ alpha * abs(A[2, 3]) ^ alpha))

  psi_32 <- 1 - 0.5 *
    ( (1) /
        (1 +
           1 * abs(A[2, 3]) ^ alpha +
           1 * abs(A[1, 2]) ^ alpha * abs(A[2, 3]) ^ alpha))

  psi_41 <- 1 - 0.5 *
    ( (1 +
         1 * abs(A[3, 4]) ^ alpha +
         1 * (abs(A[2, 4]) + abs(A[2, 3]) * abs(A[3, 4])) ^ alpha) /
        (1 +
           1 * abs(A[3, 4]) ^ alpha +
           1 * (abs(A[2, 4]) + abs(A[2, 3]) * abs(A[3, 4])) ^ alpha +
           1 * ( (abs(A[2, 4]) + abs(A[2, 3]) * abs(A[3, 4])) *
                   abs(A[1, 2])) ^ alpha))

  psi_42 <- 1 - 0.5 *
    ( (1 +
         1 * abs(A[3, 4]) ^ alpha) /
        (1 +
           1 * abs(A[3, 4]) ^ alpha +
           1 * (abs(A[2, 4]) + abs(A[2, 3]) * abs(A[3, 4])) ^ alpha +
           1 * ( (abs(A[2, 4]) + abs(A[2, 3]) * abs(A[3, 4])) *
                   abs(A[1, 2])) ^ alpha))

  psi_43 <- 1 - 0.5 *
    ( (1) /
        (1 +
           1 * abs(A[3, 4]) ^ alpha +
           1 * (abs(A[2, 4]) + abs(A[2, 3]) * abs(A[3, 4])) ^ alpha +
           1 * ( (abs(A[2, 4]) + abs(A[2, 3]) * abs(A[3, 4])) *
                   abs(A[1, 2])) ^ alpha))

  expect_equal(psi_coefficient(A, 1, 2, alpha),
               1)
  expect_equal(psi_coefficient(A, 1, 3, alpha),
               1)
  expect_equal(psi_coefficient(A, 1, 4, alpha),
               1)
  expect_equal(psi_coefficient(A, 2, 3, alpha),
               1)
  expect_equal(psi_coefficient(A, 2, 4, alpha),
               1)
  expect_equal(psi_coefficient(A, 3, 4, alpha),
               1)
  expect_equal(psi_coefficient(A, 4, 3, alpha),
               psi_43)
  expect_equal(psi_coefficient(A, 4, 2, alpha),
               psi_42)
  expect_equal(psi_coefficient(A, 4, 1, alpha),
               psi_41)
  expect_equal(psi_coefficient(A, 3, 2, alpha),
               psi_32)
  expect_equal(psi_coefficient(A, 3, 1, alpha),
               psi_31)
  expect_equal(psi_coefficient(A, 2, 1, alpha),
               psi_21)
})

test_that("psi matrix works", {
  expect_equal(psi_matrix(A, alpha),
               rbind(c(NA,
                       psi_coefficient(A, 1, 2, alpha),
                       psi_coefficient(A, 1, 3, alpha),
                       psi_coefficient(A, 1, 4, alpha)),
                     c(psi_coefficient(A, 2, 1, alpha),
                       NA,
                       psi_coefficient(A, 2, 3, alpha),
                       psi_coefficient(A, 2, 4, alpha)),
                     c(psi_coefficient(A, 3, 1, alpha),
                       psi_coefficient(A, 3, 2, alpha),
                       NA,
                       psi_coefficient(A, 3, 4, alpha)),
                     c(psi_coefficient(A, 4, 1, alpha),
                       psi_coefficient(A, 4, 2, alpha),
                       psi_coefficient(A, 4, 3, alpha),
                       NA)))
})
