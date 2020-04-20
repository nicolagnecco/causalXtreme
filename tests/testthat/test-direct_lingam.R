context("test-direct_lingam")

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

X5 <- simulate_data(10000, 30, 0.2, has_uniform_margins = FALSE)
# write.csv(X5$dataset, "analysis/temp_csv/X5.csv", row.names = FALSE)

X6 <- simulate_data(1000, 30, 0.2, has_uniform_margins = FALSE)
# write.csv(X6$dataset, "analysis/temp_csv/X6.csv", row.names = FALSE)

X7 <- simulate_data(100, 3, 0.1, has_confounder = TRUE,
                    has_uniform_margins = FALSE)
# write.csv(X7$dataset, "analysis/temp_csv/X7.csv", row.names = FALSE)


X8 <- simulate_data(1000, 20, 0.2, has_confounder = TRUE,
                    has_uniform_margins = FALSE)
# write.csv(X8$dataset, "analysis/temp_csv/X8.csv", row.names = FALSE)


# Run tests
test_that("direct lingam works", {
  expect_equal(direct_lingam_search(X3$dataset), c(2, 5, 4, 3, 1))
  expect_equal(direct_lingam_search(X4$dataset),
               c(2, 10, 3, 13, 14, 11, 8, 1, 12, 4, 6, 15, 7, 9, 5))
  expect_equal(direct_lingam_search(X5$dataset),
               c(16, 28, 22, 5, 3, 13, 7, 11, 12, 10, 14, 23, 25, 15, 20, 9, 6,
                 19, 21, 26, 1, 27, 30, 8, 17, 4, 24, 29, 18, 2))
  expect_equal(direct_lingam_search(X6$dataset),
               c(9, 28, 11, 1, 16, 14, 22, 3, 24, 19, 8, 25, 5, 20, 30, 2, 4,
                 17, 26, 12, 18, 6, 29, 7, 21, 15, 27, 23, 13, 10))
  expect_equal(direct_lingam_search(X7$dataset),
               c(2, 3, 1))
  expect_equal(direct_lingam_search(X8$dataset),
               c(18, 10, 7, 11, 16, 20, 4, 13, 5, 3, 2, 6, 9, 15,
                 1, 8, 14, 12, 19, 17))
  expect_equal(direct_lingam_search(3), NA)

})


