context("test-lingam_helpers")

# Define variables
set.seed(42)
x <- rnorm(100)
# write.csv(x, "analysis/temp_csv/x.csv", row.names = FALSE)

X1 <- matrix(rnorm(2 * 100), nrow = 2)
X1 <- X1[2:1, ]
# write.csv(X1, "analysis/temp_csv/X1.csv", row.names = FALSE)

X2 <- matrix(rnorm(3 * 100), nrow = 3)

X3 <- simulate_data(200, 5, 0.2)
# write.csv(X3$dataset, "analysis/temp_csv/X3.csv", row.names = FALSE)

X4 <- simulate_data(200, 15, 0.2)
# write.csv(X4$dataset, "analysis/temp_csv/X4.csv", row.names = FALSE)

set.seed(42)
X5 <- simulate_data(10000, 30, 0.2, has_uniform_margins = FALSE)
# write.csv(X5$dataset, "analysis/temp_csv/X5.csv", row.names = FALSE)

# Run tests
test_that("mentappr works", {
  expect_equal(mentappr(x), 1.441529478667)
})

test_that("pwling works", {
  expect_equal(pwling(X1), (-0.00250585943443604)^2)
  expect_error(pwling(X2))
})

test_that("direct lingam works", {
  expect_equal(direct_lingam_search(X3$dataset), c(2, 5, 4, 3, 1))
  expect_equal(direct_lingam_search(X4$dataset),
               c(2, 10, 3, 13, 14, 11, 8, 1, 12, 4, 6, 15, 7, 9, 5))
  expect_equal(direct_lingam_search(X5$dataset),
               c(17, 25, 10, 9, 4, 30, 18, 7, 1, 29, 22, 15, 14, 12, 16, 3, 11, 5,
                 28, 27, 26, 23, 13, 20, 8, 19, 2, 24, 6, 21))
})


