#' Causal gamma coefficient
#'
#' Computes the causal gamma coefficient between \code{v1} and \code{v2},
#' given the threshold \code{k}. In general, the gamma coefficient is not
#' symmetric, i.e., \code{compute_gamma(v1, v2) != compute_gamma(v2, v1)}.
#'
#' @param v1,v2 Numeric vectors. Two vectors with \code{n} observations.
#' @param k Integer. The number of extreme observations used to
#' compute the gamma coefficient. Set by default to \code{k = floor(2 * n^0.4)}.
#' @param to_rank Boolean. Are the vectors \code{v1} and \code{v2} already sorted?
#' Set by default to \code{TRUE}.
#' @param both_tails Boolean. Do you consider both tails when computing gamma?
#' Set by default to \code{TRUE}.
#' @return Numeric --- between 0.5 and 1. The gamma causal coefficient between \code{v1} and
#' \code{v2}.
compute_gamma <- function(v1, v2, k = floor(2 * n ^ 0.4), to_rank = TRUE,
                          both_tails = TRUE){
  n <- NROW(v1)   # number of observations

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


#' Causal gamma matrix
#'
#' Computes the gamma coefficient matrix given a dataset \code{dat},
#' and a threshold \code{k}.
#'
#' @param dat Numeric matrix. Dataset matrix with \code{n}
#' rows (observations) and \code{p} columns (variables).
#' @inheritParams compute_gamma
#' @return Numeric matrix. The entry \eqn{(i, j)} contains the
#' gamma coefficient between the \eqn{i}-th and the \eqn{j}-th
#' column of \code{dat}. If \eqn{i = j}, the value of the entry \eqn{(i, j)}
#' is set to \code{NA}.
compute_gamma_matrix <- function(dat, k = floor(2 * n ^ 0.4), both_tails = TRUE){
  n   <- NROW(dat) # number of observations
  p   <- NCOL(dat) # number of variables

  # rank variables
  ranked_dat <- apply(dat, 2, rank, ties.method = "first")

  # compute gamma for all combinations of indices
  gamma_mat <- sapply(1:p,
                      function(j){
                        sapply(1:p,
                               function(i){
                                 if (i == j){
                                   NA
                                 } else {
                                   compute_gamma(ranked_dat[, i],
                                                 ranked_dat[, j], k,
                                                 to_rank = FALSE,
                                                 both_tails = both_tails)
                                 }
                               })
                      })

  # return data
  return(gamma_mat)
}


#' Pairwise theoretical gamma coefficient
#'
#' Computes the theoretical gamma coefficient between node \code{i}
#' and node \code{j}, given the adjacency matrix \code{adj_mat},
#' the tail-index \code{alpha} and the scale of the
#' noise variables \code{noise_w}.
#'
#' @param adj_mat Numeric matrix. The adjacency matrix of a DAG.
#' All the entries must be non-negative.
#' @param i,j Integer. The indices of the nodes to consider to compute
#' the theoretical gamma coefficient.
#' @param alpha Numeric. The tail-index of the noise distribution
#' @param noise_w Numeric vector. A \eqn{p}-dimensional vector containing
#' the scale coefficients of the noise variables. Note that \eqn{p} is the
#' number of variables (nodes) in the DAG. By default, all entries are
#' set equal to one.
#' @return Numeric --- between 0.5 and 1. The theoretical gamma coefficient
#' between node \code{i} and \code{j} in the DAG \code{adj_mat}.
compute_gamma_theo <- function(adj_mat, i, j, alpha, noise_w = rep(1, p)){

  # check that weighted adjacency matrix has no negative entries
  if (min(adj_mat, na.rm = T) < 0){
    stop("The weighted adjacency matrix must have non-negative weights.")
    # !!! allow for positive and negative betas
  }

  # compute adjacency matrix
  g <- (adj_mat != 0) * 1

  # compute theoretical gamma
  p <- NROW(adj_mat)
  noise_w <- matrix(noise_w, ncol = 1)
  noise_w_alpha <- noise_w ^ alpha

  ancestor <- get_ancestors(g)
  an_node1 <- ancestor[, i]
  an_node2 <- ancestor[, j]
  nan_node2 <- abs(an_node2 - 1)

  v_num <- an_node1 * nan_node2 * noise_w_alpha
  v_denom <- an_node1 * noise_w_alpha
  e <- numeric(p)
  e[i] <- 1

  path_matrix <- t(get_all_paths(adj_mat))
  p_ij <- (e %*% path_matrix ^ alpha %*% v_num) /
    (e %*% path_matrix ^ alpha %*% v_denom)
  gamma <- drop(1 - 0.5 * p_ij)

  return(gamma)
}


#' Theoretical gamma coefficient of a DAG
#'
#' Computes the theoretical gamma coefficient matrix given the weighted
#' adjacency matrix \code{adj_mat} of an underlying DAG,
#' the tail-index \code{alpha} and the scales of the noise variables
#' \code{noise_w}.
#'
#' @inheritParams compute_gamma_theo
#' @return Numeric matrix. The entry \eqn{(i, j)} contains the theoretical
#' gamma coefficient between the \eqn{i}-th and the \eqn{j}-th node
#' of the weighted adjacency matrix \code{adj_mat}.
#' If \eqn{i = j}, the value of the entry \eqn{(i, j)} is set
#' to \code{NA}.
compute_gamma_theo_matrix <- function(adj_mat, alpha, noise_w){

  if (min(adj_mat, na.rm = T) < 0){
    stop("The weighted adjacency matrix must have non-negative weights.")
    # !!! allow for positive and negative betas
  }

  p   <- NCOL(adj_mat) # number of variables

  # compute theoretical gamma for all combinations of indices
  gamma_mat <- sapply(1:p,
                      function(j){
                        sapply(1:p,
                               function(i){
                                 if (i == j){
                                   NA
                                 } else {
                                   compute_gamma_theo(adj_mat, i, j,
                                                      alpha, noise_w)
                                 }
                               })
                      })

  # return data
  return(gamma_mat)
}
