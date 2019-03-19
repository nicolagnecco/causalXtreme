context("test-causal_extreme_functions")

# Define variables
v1     <- c(0, 2, 6, 1, 9, 8, 7, 5, 3, 4)
v1_rep <- c(0, 2, 2, 1, 2, 8, 7, 5, 3, 4)
v2     <- c(5, 6, 4, 2, 0, 8, 9, 1, 7, 3)
v2_rep <- c(5, 6, 4, 2, 0, 8, 8, 8, 7, 3)
mat <- cbind(v1, v1_rep, v2)

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

node1 <- 2
node2 <- 4
alpha <- 1.5
A <- rbind(c(0, 1.5, 0, 0),
           c(0, 0, 1.2, 2),
           c(0, 0, 0, 0.3),
           c(0, 0, 0, 0))
noise_w <- matrix(c(1, (1/2)^(1/alpha), (1/2)^(1/alpha), (1/2)^(1/alpha)), ncol = 1)


p <- NROW(A)
noise_w_alpha <- noise_w ^ alpha

ancestor <- get_ancestors(A)
an_node1 <- ancestor[, node1]; an_node1
an_node2 <- ancestor[, node2]; an_node2
nan_node2 <- abs(ancestor[, node2] - 1); nan_node2

v_num <- an_node1 * nan_node2 * noise_w_alpha
v_denom <- an_node1 * noise_w_alpha
e <- numeric(p)
e[node1] <- 1

p_ij <- (e %*% solve(diag(p) - t(A))^alpha %*% v_num) /
  (e %*% solve(diag(p) - t(A))^alpha %*% v_denom)
gamma <- drop(1 - 0.5 * p_ij); gamma

gamma_21 <- 1 - 0.5 *
  ((noise_w[2] ^ alpha) /
     (noise_w[2] ^ alpha + noise_w[1] ^ alpha * A[1, 2] ^ alpha))
gamma_31 <- 1 - 0.5 *
  ((noise_w[3] ^ alpha + noise_w[2] ^ alpha * A[2, 3] ^ alpha) /
     (noise_w[3] ^ alpha +
        noise_w[2] ^ alpha * A[2, 3] ^ alpha +
        noise_w[1] ^ alpha * A[1, 2] ^ alpha * A[2, 3] ^ alpha))
gamma_32 <- 1 - 0.5 *
  ((noise_w[3] ^ alpha) /
     (noise_w[3] ^ alpha +
        noise_w[2] ^ alpha * A[2, 3] ^ alpha +
        noise_w[1] ^ alpha * A[1, 2] ^ alpha * A[2, 3] ^ alpha))
gamma_41 <- 1 - 0.5 *
  ((noise_w[4] ^ alpha +
      noise_w[3] ^ alpha * A[3, 4] ^ alpha +
      noise_w[2] ^ alpha * (A[2, 4] + A[2, 3] * A[3, 4]) ^ alpha) /
     (noise_w[4] ^ alpha +
        noise_w[3] ^ alpha * A[3, 4] ^ alpha +
        noise_w[2] ^ alpha * (A[2, 4] + A[2, 3] * A[3, 4]) ^ alpha +
        noise_w[1] ^ alpha * ((A[2, 4] + A[2, 3] * A[3, 4]) *
                                A[1, 2]) ^ alpha))
gamma_42 <- 1 - 0.5 *
  ((noise_w[4] ^ alpha +
      noise_w[3] ^ alpha * A[3, 4] ^ alpha) /
     (noise_w[4] ^ alpha +
        noise_w[3] ^ alpha * A[3, 4] ^ alpha +
        noise_w[2] ^ alpha * (A[2, 4] + A[2, 3] * A[3, 4]) ^ alpha +
        noise_w[1] ^ alpha * ((A[2, 4] + A[2, 3] * A[3, 4]) *
                                A[1, 2]) ^ alpha))

gamma_43 <- 1 - 0.5 *
  ((noise_w[4] ^ alpha) /
     (noise_w[4] ^ alpha +
        noise_w[3] ^ alpha * A[3, 4] ^ alpha +
        noise_w[2] ^ alpha * (A[2, 4] + A[2, 3] * A[3, 4]) ^ alpha +
        noise_w[1] ^ alpha * ((A[2, 4] + A[2, 3] * A[3, 4]) *
                                A[1, 2]) ^ alpha))

gamma_34 <- 1

gamma_24 <-1

gamma_23 <- 1

gamma_12 <- gamma_13 <- gamma_14 <- 1
