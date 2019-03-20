context("test-causal_inference_helpers")

gamma <- rbind(c(NA, 1, 1, 1),
               c(.89, NA, 1, 1),
               c(.83, .93, NA, 1),
               c(.86, .96, .97, NA))
delta <- gamma - t(gamma)

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

test_that("fast search works", {})

test_that("fast search works", {})

test_that("random search works", {})

test_that("minimax search works", {})

test_that("oracle search works", {})

test_that("lingam search works", {})

test_that("pc search works", {})
