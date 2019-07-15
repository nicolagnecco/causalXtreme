context("test-causal_inference_functions")

# Define variables
method_list <-  c("fast", "full", "greedy",
                  "lingam", "maxmin", "minimax",
                  "oracle", "pc", "random")

g <- rbind(c(0, 1, 0, 0),
           c(0, 0, 1, 1),
           c(0, 0, 0, 1),
           c(0, 0, 0, 0))
est_g <- rbind(c(0, 1, 0, 0),
               c(0, 0, 1, 0),
               c(0, 0, 0, 0),
               c(0, 0, 1, 0))
order <- sample(1:NROW(g))
gamma <- rbind(c(NA, 1, 1, 1),
               c(.89, NA, 1, 1),
               c(.83, .93, NA, 1),
               c(.86, .96, .97, NA))
delta <- gamma - t(gamma)
n <- 1e2
set.seed(1)
X1 <- rnorm(n)
X2 <- X1 + rnorm(n)
X3 <- X2 + rnorm(n)
X4 <- X2 + X3 + rnorm(n)
mat <- cbind(X1, X2, X3, X4)
mat2 <- cbind(X1, X1, X1, X1)
mat3 <-  matrix("NA", nrow = 3, ncol = 10)

# Run tests
# test_that("causal discovery works", {
#
#   order <- fast_perm_search(delta, "sum")$order
#   est_g <- caus_order_to_dag(order)
#   ll <- list(est_g = est_g, order = order)
#   expect_equal(causal_discovery("fast", delta = delta), ll)
#   expect_equal(causal_discovery("fast", g = g, gamma = gamma, mat = mat,
#                                 delta = delta), ll)
#
#   order <- full_perm_search(delta)$order
#   est_g <- caus_order_to_dag(order)
#   ll <- list(est_g = est_g, order = order)
#   expect_equal(causal_discovery("full", delta = delta), ll)
#   expect_equal(causal_discovery("full", gamma = gamma, delta = delta), ll)
#
#   order <- greedy_perm_search(delta)$order
#   est_g <- caus_order_to_dag(order)
#   ll <- list(est_g = est_g, order = order)
#   expect_equal(causal_discovery("greedy", delta = delta), ll)
#   expect_equal(causal_discovery("greedy", delta = delta, g = g), ll)
#
#   order <- fast_perm_search(delta, "maxmin")$order
#   est_g <- caus_order_to_dag(order)
#   ll <- list(est_g = est_g, order = order)
#   expect_equal(causal_discovery("maxmin", delta = delta), ll)
#   expect_equal(causal_discovery("maxmin", mat = mat, g = g, delta = delta), ll)
#
#   est_g <- lingam_search(mat)
#   order <- compute_caus_order(est_g)
#   ll <- list(est_g = est_g, order = order)
#   expect_equal(causal_discovery("lingam", mat = mat), ll)
#   expect_equal(causal_discovery("lingam", gamma = gamma, mat = mat), ll)
#
#   est_g <- lingam_search(mat2)
#   order <- NA
#   ll <- list(est_g = est_g, order = order)
#   expect_equal(causal_discovery("lingam", mat = mat2), ll)
#   expect_equal(causal_discovery("lingam", gamma = gamma, mat = mat2), ll)
#
#   est_g <- pc_search(mat)
#   dag <- est_g * (t(est_g) == 0)
#   order <- compute_caus_order(dag)
#   ll <- list(est_g = est_g, order = order)
#   expect_equal(causal_discovery("pc", mat = mat), ll)
#   expect_equal(causal_discovery("pc", delta = delta, mat = mat), ll)
#
#   est_g <- pc_search(mat3)
#   dag <- est_g * (t(est_g) == 0)
#   order <- NA
#   ll <- list(est_g = est_g, order = order)
#   expect_equal(causal_discovery("pc", mat = mat3), ll)
#   expect_equal(causal_discovery("pc", delta = delta, mat = mat3), ll)
#
#   order <- oracle_search(g)
#   est_g <- caus_order_to_dag(order)
#   ll <- list(est_g = est_g, order = order)
#   expect_equal(causal_discovery("oracle", g = g), ll)
#   expect_equal(causal_discovery("oracle", delta = delta, g = g), ll)
#
#   set.seed(1)
#   order <- random_perm_search(g)
#   est_g <- caus_order_to_dag(order)
#   ll <- list(est_g = est_g, order = order)
#   set.seed(1)
#   expect_equal(causal_discovery("random", g = g), ll)
#   set.seed(1)
#   expect_equal(causal_discovery("random", g = g, mat = mat, gamma = gamma), ll)
#
#   order <- minimax_search(gamma)
#   est_g <- caus_order_to_dag(order)
#   ll <- list(est_g = est_g, order = order)
#   expect_equal(causal_discovery("minimax", gamma = gamma), ll)
#
#   expect_equal(causal_discovery("minimax", mat = mat, gamma = gamma, g = g), ll)
#
#   expect_error(causal_discovery("miniiiis", gamma = gamma))
#
#   for (m in method_list){
#     expect_error(causal_discovery(m))
#   }
#
# })
#
# test_that("causal metrics work", {
#   expect_equal(causal_metrics(g, est_g, order),
#                list(str_int_dist = compute_str_int_distance(g, est_g),
#                     ancestral_dist = compute_ancestral_distance(g, order)))
#   expect_equal(causal_metrics(g, NA, order),
#                list(str_int_dist = NA,
#                     ancestral_dist = compute_ancestral_distance(g, order)))
#   expect_equal(causal_metrics(g, est_g, NA),
#                list(str_int_dist = compute_str_int_distance(g, est_g),
#                     ancestral_dist = NA))
#   expect_equal(causal_metrics(g, NA, NA),
#                list(str_int_dist = NA,
#                     ancestral_dist = NA))
#   expect_error(causal_metrics(NA, est_g, order))
# })
