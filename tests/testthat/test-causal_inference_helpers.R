context("test-causal_inference_helpers")

# Define variables
g <- rbind(c(0, 1, 0, 0),
           c(0, 0, 1, 1),
           c(0, 0, 0, 1),
           c(0, 0, 0, 0))
gamma <- rbind(c(NA, 1, 1, 1),
               c(.89, NA, 1, 1),
               c(.83, .93, NA, 1),
               c(.86, .96, .97, NA))
delta <- gamma - t(gamma)
n <- 1e2
set.seed(1)
X1 <- rnorm(n)
X2 <- X1 + rnorm(n)
X3 <- X1 + rnorm(n)
dat1 <- cbind(X1, X2, X3)
dat2 <- cbind(X1, X1, X1)
dat3 <- matrix("NA", nrow = 3, ncol = 10)

# Run tests
test_that("greedy search works", {
  expect_equal(greedy_perm_search(delta),
               list(order = c(1, 2, 3, 4),
                    score = sum(delta[upper.tri(delta)])))
  expect_equal(greedy_perm_search(delta, silent = TRUE),
               list(order = c(1, 2, 3, 4),
                    score = sum(delta[upper.tri(delta)])))
  expect_error(greedy_perm_search(delta[1:2, 1:2], silent = TRUE))
  expect_error(greedy_perm_search(delta[1:2, 1:2], silent = FALSE))
})

test_that("fast search works", {
  expect_equal(fast_perm_search(delta, mode = "sum", silent = FALSE),
               list(order = c(1, 2, 3, 4),
                    score = sum(delta[upper.tri(delta)])))
  expect_equal(fast_perm_search(delta, mode = "sum", silent = TRUE),
               list(order = c(1, 2, 3, 4),
                    score = sum(delta[upper.tri(delta)])))
  expect_equal(fast_perm_search(delta, mode = "maxmin", silent = FALSE),
               list(order = c(1, 2, 3, 4),
                    score = sum(delta[upper.tri(delta)])))
  expect_equal(fast_perm_search(delta, mode = "maxmin", silent = TRUE),
               list(order = c(1, 2, 3, 4),
                    score = sum(delta[upper.tri(delta)])))
  expect_error(fast_perm_search(delta, mode = "abc", silent = FALSE))
  expect_error(fast_perm_search(delta, mode = "abc", silent = TRUE))
})

test_that("full search works", {
  expect_equal(full_perm_search(delta),
              list(order = c(1, 2, 3, 4),
                   score = sum(delta[upper.tri(delta)])))
  expect_equal(full_perm_search(delta, silent = TRUE),
              list(order = c(1, 2, 3, 4),
                   score = sum(delta[upper.tri(delta)])))
  expect_error(full_perm_search(matrix(0, nrow = 11, ncol = 11)))
})

test_that("random search works", {
  set.seed(1)
  order <- sample(1:NROW(g))
  set.seed(1)

  expect_equal(random_perm_search(g), order)
})

test_that("minimax search works", {
  expect_equal(minimax_search(gamma), c(1, 2, 3, 4))
})

test_that("oracle search works", {
  expect_equal(oracle_search(g), c(1, 2, 3, 4))
})

test_that("lingam search works", {
  set.seed(1)
  out <- pcalg::lingam(dat1)
  dag <- t( (out$Bpruned != 0) * 1)
  set.seed(1)

  expect_equal(lingam_search(dat1), dag)
  expect_equal(lingam_search(dat2), NA)
})

test_that("pc search works", {
  out <- pcalg::pc(suffStat = list(C = cor(dat1), n = NROW(dat1)),
                   indepTest = pcalg::gaussCItest,
                   p = NCOL(dat1),
                   alpha = 0.05,
                   u2pd = "retry",
                   skel.method = "stable")
  cpdag <- as(out@graph, "matrix")

  expect_equal(pc_search(dat1), cpdag)

  expect_equal(pc_search(dat3), NA)
  })
