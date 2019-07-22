#' Estimate causal tail coefficient
#'
#' Estimates the causal tail coefficient between two variables
#' \code{v1} and \code{v2}, given the threshold \code{k}.
#'
#' The causal tail coefficient is defined in
#' the paper "Causality in heavy-tailed models" and has two formulations.
#' \enumerate{
#' \item The first formulation is defined in the paper as the
#' \eqn{\Gamma}-coefficient, and it considers only the upper tails of the
#' variables. To use this formulation set \code{both_tails = FALSE}.
#' \item The second formulation is defined in the paper as the
#' \eqn{\Psi}-coefficient, and considers both tails of the variables.
#' To use this formulation set \code{both_tails = FALSE}.
#' }
#' In general, the causal tail coefficient is not symmetric, i.e.,
#' \code{causal_tail_coeff(v1, v2) != causal_tail_coeff(v2, v1)}.
#'
#' @param v1,v2 Numeric vectors. Two vectors with \code{n} observations.
#' @param k Positive integer. The number of extreme observations used to
#' compute the causal tail coefficient. Set by default to
#' \code{k = floor(2 * n^0.4)}. It must be greater than 1 and smaller
#' than \code{n}.
#' @param to_rank Boolean. Are the vectors \code{v1} and \code{v2}
#' already sorted?
#' Set by default to \code{TRUE}.
#' @param both_tails Boolean. Do you consider both tails when computing the
#' causal tail coefficient? If \code{both_tails = TRUE}, then you use the
#' \eqn{\Psi} formulation. If \code{both_tails = FALSE},
#' then you use the \eqn{\Gamma} formulation.
#' Set by default to \code{TRUE}.
#' @return Numeric --- between 0 and 1.
#' The causal tail coefficient between \code{v1} and \code{v2}.
#' @export
causal_tail_coeff <- function(v1, v2, k = floor(2 * n ^ 0.4), to_rank = TRUE,
                              both_tails = TRUE){
  # number of observations
  n <- NROW(v1)

  # check k
  if (k <= 1 | k >= n) {
    stop("k must be greater than 1 and smaller than n.")
  }

  # rank variables?
  if (to_rank){
    r1 <- rank(v1, ties.method = "first") # ranks of v1
    r2 <- rank(v2, ties.method = "first") # ranks of v2
  } else{
    r1 <- v1
    r2 <- v2
  }

  # compute causal tail coefficient
  if (both_tails){
    k <- (k %/% 2) * 2
    1 / (k * n) * sum(2 * abs(r2[r1 > n - k / 2 | r1 <= k / 2] - (n + 1) / 2))
  } else{
    1 / (k * n) * sum(r2[r1 > n - k])
  }
}


#' Estimate causal tail matrix
#'
#' Estimates the matrix with causal tail coefficients,
#' given a dataset \code{dat}, and a threshold \code{k}.
#'
#' For more information see the documentation of
#' \code{\link{causal_tail_coeff}}.
#'
#' @param dat Numeric matrix. Dataset matrix with \code{n}
#' rows (observations) and \code{p} columns (variables).
#' @inheritParams causal_tail_coeff
#' @return Numeric matrix. The entry \eqn{(i, j)} contains the
#' causal tail coefficient between the \eqn{i}-th and the \eqn{j}-th
#' column of \code{dat}. The values on the main diagonal are set to \code{NA}.
#' @export
causal_tail_matrix <- function(dat, k = floor(2 * n ^ 0.4),
                               both_tails = TRUE){

  # get number of observations and variables
  n   <- NROW(dat)
  p   <- NCOL(dat)

  # rank variables
  ranked_dat <- apply(dat, 2, rank, ties.method = "first")

  # compute causal tail coefficient for all pairs of indices
  causal_mat <- sapply(1:p,
                       function(j){
                         sapply(1:p,
                                function(i){
                                  if (i == j){
                                    NA
                                  } else {
                                    causal_tail_coeff(ranked_dat[, i],
                                                      ranked_dat[, j], k,
                                                      to_rank = FALSE,
                                                      both_tails = both_tails)
                                  }
                                })
                       })

  # return data
  return(causal_mat)
}


#' Compute \eqn{\Psi}-coefficient
#'
#' Computes the theoretical \eqn{\Psi}-coefficient defined
#' in the paper "Causality in heavy-tailed models".
#'
#' For more information see the documentation of
#' \code{\link{psi_matrix}}.
#'
#' @param adj_mat Square numeric matrix. The adjacency matrix of a DAG.
#' If \code{both_tails = FALSE}, all the entries must be non-negative.
#' @param i,j Integer. The indices of the nodes to consider to compute
#' the theoretical gamma coefficient.
#' @param tail_index Numeric. The tail-index of the noise distribution.
#' @return Numeric --- between 0 and 1. The theoretical psi coefficient
#' between node \code{i} and \code{j} in the DAG \code{adj_mat}.
#' @export
psi_coefficient <- function(adj_mat, i, j, tail_index){

  # compute DAG
  dag <- (adj_mat != 0) * 1

  # compute theoretical psi
  p <- NROW(adj_mat)

  ancestor <- get_ancestors(dag)
  an_node1 <- ancestor[, i]
  an_node2 <- ancestor[, j]

  e_ij <- an_node1 * an_node2
  e_an_i <- an_node1
  e <- numeric(p)
  e[i] <- 1

  path_matrix <- t(get_all_paths(adj_mat))
  frac <- (e %*% abs(path_matrix) ^ tail_index %*% e_ij) /
    (e %*% abs (path_matrix) ^ tail_index %*% e_an_i)
  psi <- drop(1 / 2 + 1 / 2 * frac)

  return(psi)
}


#' Compute \eqn{\Psi}-matrix
#'
#' Computes the theoretical \eqn{\Psi}-matrix defined
#' in the paper "Causality in heavy-tailed models".
#'
#' The causal tail coefficient is defined in
#' the paper "Causality in heavy-tailed models" and has two formulations.
#' This function considers the \eqn{\Psi}-formulation, and assumes that the
#' scales of the noise variables are equal to one.
#'
#' @inheritParams psi_coefficient
#' @return Numeric matrix. The entry \eqn{(i, j)} contains the
#' \eqn{\Psi} coefficient between node \eqn{i} and node \eqn{j}.
#' The values on the main diagonal are set to \code{NA}.
#' @export
psi_matrix <- function(adj_mat, tail_index){

  # number of variables
  p   <- NCOL(adj_mat)

  # compute psi coefficient for all pairs of indices
  psi_mat <- sapply(1:p,
                      function(j){
                        sapply(1:p,
                               function(i){
                                 if (i == j){
                                   NA
                                 } else {
                                   psi_coefficient(adj_mat, i, j, tail_index)
                                 }
                               })
                      })

  # return data
  return(psi_mat)
}
