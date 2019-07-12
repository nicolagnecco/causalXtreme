# To do !!!

# CPDAG -> DAG
# possible to extend
test1 <- pcalg::pdag2dag(truecpdag)
test2 <- pcalg::pdag2dag(cpdag1)

# impossible to extend.
impcpdag <- as(rbind(c(0, 1, 0, 1),
                     c(1, 0, 1, 0),
                     c(0, 1, 0, 1),
                     c(1, 0, 1, 0)), "graphNEL")
test3 <- pcalg::pdag2dag(impcpdag)

if (test3$success){
  as(test3$graph, "matrix")
}

# PC algorithm

# NORMAL PC
# pc.fit <- pc(suffStat = suffStat, indepTest = gaussCItest, p = p, alpha = 5e-3,
#              u2pd = "retry", skel.method = "stable")
#

# RANK PC
# suffStat.data <- list(C = 2 * sin(cor(X, method = "spearman") *
#                                     pi/6), n = nrow(X))
# setOptions$indepTest <- pcalg::gaussCItest
# rankpc.fit <- pc(suffStat = suffStat, indepTest = gaussCItest, p = p, alpha = 5e-3,
#              u2pd = "retry", skel.method = "stable")
