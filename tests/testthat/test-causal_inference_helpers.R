context("test-causal_inference_helpers")

# Define variables
u <- 1321
set.seed(u)
n <- 500
X1 <- rt(n, 2.5)
X3 <- X1 + rt(n, 2.5)
X2 <- X1 + X3 + rt(n, 2.5)
dat1 <- cbind(X1, X2, X3)
dat2 <- cbind(X1, X1, X1)
dat3 <- matrix("NA", nrow = 3, ncol = 10)
dat4 <- simulate_data(n, sample(2:20, 1), prob_connect = runif(1))$dataset

# Run tests
test_that("random search works", {
  set.seed(u)
  dag <- random_dag(p = NCOL(dat1), prob_connect = runif(1))
  set.seed(u)
  expect_equal(random_search(dat1), dag)
})

test_that("greedy ancestral search works", {
  set.seed(u)
  expect_equal(greedy_ancestral_search(dat1), c(1, 3, 2))
  expect_equal(greedy_ancestral_search(dat2), c(1, 2, 3))
})

test_that("lingam search works", {
  set.seed(u)
  out <- pcalg::lingam(dat1)
  dag <- t( (out$Bpruned != 0) * 1)
  set.seed(u)

  expect_equal(lingam_search(dat1, "logcosh"), dag)
  expect_equal(lingam_search(dat2), NA)

  set.seed(u)
  t.k <- estLiNGAM(dat4, only.perm = T, fun = "exp")$k
  out <- prune(t(dat4), k = t.k)$Bpruned
  dag2 <- (t(out) != 0) * 1
  set.seed(u)
  expect_equal(lingam_search(dat4, contrast_fun = "exp"), dag2)
})

test_that("pc search works", {
  out <- pcalg::pc(suffStat = list(C = cor(dat1), n = NROW(dat1)),
                   indepTest = pcalg::gaussCItest,
                   p = NCOL(dat1),
                   alpha = 0.05,
                   u2pd = "retry",
                   skel.method = "stable")
  cpdag <- as(out@graph, "matrix")

  expect_equal(pc_search(dat1, 5e-2), cpdag)

  expect_equal(pc_search(dat3, 5e-2), NA)
})

test_that("pc rank search works", {
  suff_stat <-  list(C = 2 * sin(cor(dat4, method = "spearman") * pi / 6),
                     n = nrow(dat4))
  out <- pcalg::pc(suffStat = suff_stat,
                   indepTest = pcalg::gaussCItest,
                   p = NCOL(dat4),
                   alpha = 0.05,
                   u2pd = "retry",
                   skel.method = "stable")
  cpdag <- as(out@graph, "matrix")

  expect_equal(pc_rank_search(dat4, 5e-2), cpdag)
  expect_equal(pc_rank_search(dat3, 5e-2), NA)
})
