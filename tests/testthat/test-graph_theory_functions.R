context("test-graph_theory_functions")

# Define variables
dag1 <- rbind(c(0, 0, 1, 0),
              c(1, 0, 0, 0),
              c(0, 0, 0, 1),
              c(0, 0, 0, 0))

dag2 <- rbind(c(0, 0, 0, 0),
              c(0, 0, 0, 0),
              c(0, 0, 0, 1),
              c(1, 0, 0, 0))

adj_mat3 <- rbind(c(0, 1.5, 0, 0),
                  c(0, 0, -1.2, 2),
                  c(0, 0, 0, 0.3),
                  c(0, 0, 0, 0))

dag1_anc <- rbind(c(1, 0, 1, 1),
                  c(1, 1, 1, 1),
                  c(0, 0, 1, 1),
                  c(0, 0, 0, 1))

dag2_anc <- rbind(c(1, 0, 0, 0),
                  c(0, 1, 0, 0),
                  c(1, 0, 1, 1),
                  c(1, 0, 0, 1))

dag3_anc <- rbind(c(1, 1, 1, 1),
                  c(0, 1, 1, 1),
                  c(0, 0, 1, 1),
                  c(0, 0, 0, 1))

dag3_path_count <- rbind(c(1, 1, 1, 2),
                         c(0, 1, 1, 2),
                         c(0, 0, 1, 1),
                         c(0, 0, 0, 1))

dag3_path_weighted <- rbind(c(1, 1.5, -1.5 * 1.2, 1.5 * (-1.2 * 0.3 + 2)),
                            c(0, 1, -1.2, -1.2 * 0.3 + 2),
                            c(0, 0, 1, 0.3),
                            c(0, 0, 0, 1))

dag4 <- rbind(c(0, 0, 1),
              c(0, 0, 0),
              c(0, 1, 0))

dag5 <- rbind(c(0, 0, 1),
              c(0, 0, 1),
              c(0, 0, 0))

dag6 <- rbind(c(0, 0, 1),
              c(1, 0, 1),
              c(0, 0, 0))

dag7 <- rbind(c(0, 1, 0),
              c(0, 0, 0),
              c(1, 1, 0))

dag8 <- rbind(c(0, 0, 1, 0, 0, 1),
              c(0, 0, 1, 1, 0, 0),
              c(0, 0, 0, 0, 1, 0),
              rep(0, 6),
              c(rep(0, 5), 1),
              rep(0, 6))

cpdag8 <- rbind(c(0, 0, 1, 0, 0, 1),
                c(0, 0, 1, 1, 0, 0),
                c(0, 0, 0, 0, 1, 0),
                c(0, 1, 0, 0, 0, 0),
                c(rep(0, 5), 1),
                rep(0, 6))

dag9 <- rbind(c(0, 0, 1, 0, 0, 1),
              c(0, 0, 1, 0, 0, 0),
              c(0, 0, 0, 0, 1, 0),
              c(0, 1, 0, 0, 0, 0),
              c(rep(0, 5), 1),
              rep(0, 6))

cpdag9 <- rbind(c(0, 0, 1, 0, 0, 1),
                c(0, 0, 1, 1, 0, 0),
                c(0, 0, 0, 0, 1, 0),
                c(0, 1, 0, 0, 0, 0),
                c(rep(0, 5), 1),
                rep(0, 6))

dag10 <- rbind(c(0, 0, 1, 0, 0, 1),
               c(0, 0, 1, 1, 0, 0),
               c(0, 0, 0, 0, 1, 0),
               rep(0, 6),
               rep(0, 6),
               c(0, 0, 0, 0, 1, 0))

cpdag10 <- rbind(c(0, 0, 1, 0, 0, 1),
                 c(0, 0, 1, 1, 0, 0),
                 c(0, 0, 0, 0, 1, 0),
                 c(0, 1, 0, 0, 0, 0),
                 rep(0, 6),
                 c(1, 0, 0, 0, 1, 0))

# Run tests
test_that("the computed causal order is correct", {
  expect_equal(compute_caus_order(dag1), c(2, 1, 3, 4))
  expect_equal(compute_caus_order(dag2), c(2, 3, 4, 1))
  expect_error(compute_caus_order(adj_mat3))
})

test_that("the checked causal order is correct", {
  expect_equal(check_caus_order(c(1, 2, 3, 4), dag1), FALSE)
  expect_equal(check_caus_order(c(2, 1, 3, 4), dag1), TRUE)
  expect_equal(check_caus_order(c(2, 3, 4, 1), dag2), TRUE)
  expect_equal(check_caus_order(c(3, 2, 4, 1), dag2), TRUE)
  expect_equal(check_caus_order(c(3, 4, 2, 1), dag2), TRUE)
  expect_equal(check_caus_order(c(3, 4, 1, 2), dag2), TRUE)
  expect_equal(check_caus_order(c(3, 1, 4, 2), dag2), FALSE)
  expect_equal(check_caus_order(c(1, 4, 2, 3), dag2), FALSE)
  expect_equal(check_caus_order(c(4, 3, 2, 1), dag2), FALSE)
  expect_error(check_caus_order(c(4, 3, 2, 1), adj_mat3))
  expect_error(check_caus_order(c(NA, 1), dag2))
  expect_error(check_caus_order(NA, dag2))
  expect_error(check_caus_order(c(NA, NA), dag2))
})

test_that("ancestors are correct", {
  expect_equal(get_ancestors(dag1), dag1_anc)
  expect_equal(get_ancestors(dag2), dag2_anc)
  expect_error(get_ancestors(adj_mat3))
})

test_that("descendants are correct", {
  expect_equal(get_descendants(dag1), t(dag1_anc))
  expect_equal(get_descendants(dag2), t(dag2_anc))
  expect_error(get_descendants(adj_mat3))
})

test_that("parents are correct", {
  expect_equal(get_parents(dag1), dag1)
  expect_equal(get_parents(dag2), dag2)
  expect_error(get_parents(adj_mat3))
})

test_that("children are correct", {
  expect_equal(get_children(dag1), t(dag1))
  expect_equal(get_children(dag2), t(dag2))
  expect_error(get_children(adj_mat3))
})

test_that("paths are correct", {
  expect_equal(get_all_paths( (adj_mat3 != 0) * 1), dag3_path_count)
  expect_equal(get_all_paths(adj_mat3), dag3_path_weighted)
  expect_equal(get_all_paths( (adj_mat3 != 0) * 1), dag3_path_count)
  expect_error(get_all_paths(adj_mat3, type = "foo"))
})

test_that("converting causal order into DAG works", {
  expect_equal(caus_order_to_dag(c(3, 1, 2)), dag7)
  expect_equal(caus_order_to_dag(c(2, 1, 3)), dag6)
  expect_error(caus_order_to_dag(c(NA, NA)))
  expect_error(caus_order_to_dag(c(NA, 2)))
  expect_error(caus_order_to_dag(NA))
})

test_that("converting DAG into CPDAG works", {
  expect_equal(dag_to_cpdag(dag8), cpdag8)
  expect_equal(dag_to_cpdag(dag9), cpdag9)
  expect_equal(dag_to_cpdag(dag10), cpdag10)
  expect_error(dag_to_cpdag(adj_mat3))
})

test_that("converting causal order into CPDAG works", {
  p <- sample(1:20, 1)
  dag_temp <- caus_order_to_dag(1:p)
  expect_equal(caus_order_to_cpdag(1:p),  dag_temp + t(dag_temp))
  expect_error(caus_order_to_cpdag(NA))
})

test_that("ancestral distance is correct", {
  expect_equal(compute_ancestral_distance(dag1, c(2, 1, 3, 4)), 0)
  expect_equal(compute_ancestral_distance(dag1, c(2, 1, 4, 3)), 1 / (4 * 3 / 2))
  expect_equal(compute_ancestral_distance(dag1, c(1, 2, 4, 3)), 2 / (4 * 3 / 2))
  expect_equal(compute_ancestral_distance(dag1, c(1, 4, 2, 3)), 3 / (4 * 3 / 2))
  expect_equal(compute_ancestral_distance(dag1, c(4, 3, 1, 2)), 6 / (4 * 3 / 2))
  expect_error(compute_ancestral_distance(dag1, c(NA, NA, NA, NA)))
  expect_error(compute_ancestral_distance(dag1, c(NA, NA, 2, NA)))
  expect_error(compute_ancestral_distance(dag1, NA))
  expect_equal(compute_ancestral_distance(dag2, c(4, 3, 2, 1)), 1 / (4 * 3 / 2))
  expect_equal(compute_ancestral_distance(dag2, c(4, 2, 3, 1)), 1 / (4 * 3 / 2))
  expect_equal(compute_ancestral_distance(dag2, c(3, 4, 2, 1)), 0 / (4 * 3 / 2))
  expect_error(compute_ancestral_distance(adj_mat3, c(4, 3, 2, 1)),
               6 / (4 * 3 / 2))
})

test_that("structural intervention distance is correct", {
  expect_equal(compute_str_int_distance(dag4, dag5), 3 / 6)
  expect_equal(compute_str_int_distance(dag4, dag6), 5 / 6)
  expect_false(compute_str_int_distance(dag5, dag6) ==
                 compute_str_int_distance(dag6, dag5))
  expect_error(compute_str_int_distance(dag4, NA))
  expect_error(compute_str_int_distance(dag4, adj_mat3))
  expect_error(compute_str_int_distance(adj_mat3, dag5))
  expect_error(compute_str_int_distance(NA, dag5))
})

test_that("structural Hamming distance is correct", {
  expect_equal(compute_str_ham_distance(cpdag8, cpdag9), 0 / 15)
  expect_equal(compute_str_ham_distance(cpdag8, cpdag10), 2 / 15)
  expect_true(compute_str_ham_distance(cpdag8, cpdag9) ==
                compute_str_ham_distance(cpdag9, cpdag8))
  expect_error(compute_str_ham_distance(cpdag9, NA))
  expect_error(compute_str_ham_distance(cpdag10, adj_mat3))
  expect_error(compute_str_ham_distance(NA, cpdag9))
  expect_error(compute_str_ham_distance(adj_mat3, cpdag9))

})
