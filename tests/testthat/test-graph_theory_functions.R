context("test-graph_theory_functions")

# Define variables
g1 <- rbind(c(0, 0, 1, 0),
            c(1, 0, 0, 0),
            c(0, 0, 0, 1),
            c(0, 0, 0, 0))

g2 <- rbind(c(0, 0, 0, 0),
            c(0, 0, 0, 0),
            c(0, 0, 0, 1),
            c(1, 0, 0, 0))

g3 <- rbind(c(0, 1.5, 0, 0),
            c(0, 0, -1.2, 2),
            c(0, 0, 0, 0.3),
            c(0, 0, 0, 0))

g1_anc <- rbind(c(1, 0, 1, 1),
                c(1, 1, 1, 1),
                c(0, 0, 1, 1),
                c(0, 0, 0, 1))

g2_anc <- rbind(c(1, 0, 0, 0),
                c(0, 1, 0, 0),
                c(1, 0, 1, 1),
                c(1, 0, 0, 1))

g3_anc <- rbind(c(1, 1, 1, 1),
                c(0, 1, 1, 1),
                c(0, 0, 1, 1),
                c(0, 0, 0, 1))

g3_path_count <- rbind(c(1, 1, 1, 2),
                       c(0, 1, 1, 2),
                       c(0, 0, 1, 1),
                       c(0, 0, 0, 1))

g3_path_weighted <- rbind(c(1, 1.5, -1.5 * 1.2, 1.5 * (-1.2 * 0.3 + 2)),
                          c(0, 1, -1.2, -1.2 * 0.3 + 2),
                          c(0, 0, 1, 0.3),
                          c(0, 0, 0, 1))

g4 <- rbind(c(0, 0, 1),
            c(0, 0, 0),
            c(0, 1, 0))

g5 <- rbind(c(0, 0, 1),
            c(0, 0, 1),
            c(0, 0, 0))

g6 <- rbind(c(0, 0, 1),
            c(1, 0, 1),
            c(0, 0, 0))

g7 <- rbind(c(0, 1, 0),
            c(0, 0, 0),
            c(1, 1, 0))


# Run tests
test_that("the computed causal order is correct", {
  expect_equal(compute_caus_order(g1), c(2, 1, 3, 4))
  expect_equal(compute_caus_order(g2), c(2, 3, 4, 1))
  expect_error(compute_caus_order(g3))
})

test_that("the checked causal order is correct", {
  expect_equal(check_caus_order(c(1, 2, 3, 4), g1), FALSE)
  expect_equal(check_caus_order(c(2, 1, 3, 4), g1), TRUE)
  expect_equal(check_caus_order(c(2, 3, 4, 1), g2), TRUE)
  expect_equal(check_caus_order(c(3, 2, 4, 1), g2), TRUE)
  expect_equal(check_caus_order(c(3, 4, 2, 1), g2), TRUE)
  expect_equal(check_caus_order(c(3, 4, 1, 2), g2), TRUE)
  expect_equal(check_caus_order(c(3, 1, 4, 2), g2), FALSE)
  expect_equal(check_caus_order(c(1, 4, 2, 3), g2), FALSE)
  expect_equal(check_caus_order(c(4, 3, 2, 1), g2), FALSE)
  expect_error(check_caus_order(c(4, 3, 2, 1), g3))
  expect_error(check_caus_order(c(NA, 1), g2))
  expect_error(check_caus_order(NA, g2))
  expect_error(check_caus_order(c(NA, NA), g2))
})

test_that("ancestors are correct", {
  expect_equal(get_ancestors(g1), g1_anc)
  expect_equal(get_ancestors(g2), g2_anc)
  expect_error(get_ancestors(g3))
})

test_that("descendants are correct", {
  expect_equal(get_descendants(g1), t(g1_anc))
  expect_equal(get_descendants(g2), t(g2_anc))
  expect_error(get_descendants(g3))
})

test_that("parents are correct", {
  expect_equal(get_parents(g1), g1)
  expect_equal(get_parents(g2), g2)
  expect_error(get_parents(g3))
})

test_that("children are correct", {
  expect_equal(get_children(g1), t(g1))
  expect_equal(get_children(g2), t(g2))
  expect_error(get_children(g3))
})

test_that("paths are correct", {
  expect_equal(get_all_paths(g3, type = "count"), g3_path_count)
  expect_equal(get_all_paths(g3, type = "weighted"), g3_path_weighted)
  expect_equal(get_all_paths(g3), g3_path_count)
  expect_error(get_all_paths(g3, type = "foo"))
})

test_that("ancestral distance is correct", {
  expect_equal(compute_ancestral_distance(g1, c(2, 1, 3, 4)), 0)
  expect_equal(compute_ancestral_distance(g1, c(2, 1, 4, 3)), 1 / (4 * 3 / 2))
  expect_equal(compute_ancestral_distance(g1, c(1, 2, 4, 3)), 2 / (4 * 3 / 2))
  expect_equal(compute_ancestral_distance(g1, c(1, 4, 2, 3)), 3 / (4 * 3 / 2))
  expect_equal(compute_ancestral_distance(g1, c(4, 3, 1, 2)), 6 / (4 * 3 / 2))
  expect_error(compute_ancestral_distance(g1, c(NA, NA, NA, NA)))
  expect_error(compute_ancestral_distance(g1, c(NA, NA, 2, NA)))
  expect_error(compute_ancestral_distance(g1, NA))
  expect_equal(compute_ancestral_distance(g2, c(4, 3, 2, 1)), 1 / (4 * 3 / 2))
  expect_equal(compute_ancestral_distance(g2, c(4, 2, 3, 1)), 1 / (4 * 3 / 2))
  expect_equal(compute_ancestral_distance(g2, c(3, 4, 2, 1)), 0 / (4 * 3 / 2))
  expect_error(compute_ancestral_distance(g3, c(4, 3, 2, 1)), 6 / (4 * 3 / 2))
})

test_that("structural intervention distance is correct", {
  expect_equal(compute_str_int_distance(g4, g5), 3 / 6)
  expect_equal(compute_str_int_distance(g4, g6), 5 / 6)
  expect_false(compute_str_int_distance(g5, g6) ==
                 compute_str_int_distance(g6, g5))
  expect_error(compute_str_int_distance(g4, NA))
  expect_error(compute_str_int_distance(g3, g5))
  expect_error(compute_str_int_distance(g4, g3))
})

test_that("converting causal order into DAG works", {
  expect_equal(caus_order_to_adjmat(c(3, 1, 2)), g7)
  expect_equal(caus_order_to_adjmat(c(2, 1, 3)), g6)
  expect_error(caus_order_to_adjmat(c(NA, NA)))
  expect_error(caus_order_to_adjmat(c(NA, 2)))
  expect_error(caus_order_to_adjmat(NA))
})
