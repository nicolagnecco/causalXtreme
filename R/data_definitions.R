# Data definitions


# dag is square binary_matrix that can be made strictly upper triangular
# interp. a directed acyclic graph (DAG)

dag <- rbind(c(0, 0, 0, 1),
             c(1, 0, 1, 0),
             c(0, 0, 0, 1),
             c(0, 0, 0, 0))


# cpdag is square binary_matrix
# interp. a complete partially directed acyclic graph (CPDAG)

cpdag <-  rbind(c(0, 1, 0, 1),
                c(1, 0, 1, 0),
                c(0, 0, 0, 1),
                c(1, 0, 0, 0))


# causal_order is numeric_vector
# interp. the causal order of an underlying DAG

caus_order <- c(2, 1, 3, 4)


# adjacency_matrix is square numeric_matrix that can be made strictly upper triangular
# interp. the weighted adjacency matrix of an underlying DAG

adj_mat <- rbind(c(0, 0, 0, 0.1),
                 c(-1, 0, 2.3, 0),
                 c(0, 0, 0, -0.8),
                 c(0, 0, 0, 0))


# coeff_matrix is square numeric_matrix that can be made strictly upper triangular
# interp. the transpose of an adjacency_matrix

coeff_mat <- t(adj_mat)


# dataset is numeric_matrix
# interp. a dataset with observations in the rows and variables in the columns

dat <- matrix(rt(2 * 100, df = 1.5), nrow = 100, ncol = 2)


# variable is numeric_vector
# interp. one variable, i.e., column, of a dataset

v1 <- dat[, 1]
v2 <- dat[, 2]


# observation is numeric_vector
# interp. one observation, i.e., row, of a dataset

o1 <- dat[17, ]
o2 <- dat[42, ]
