#' Causal gamma coefficient
#'
#' Computes the causal gamma coefficient between \code{v1} and \code{v2},
#' given the threshold \code{k}. In general, the gamma coefficient is not
#' symmetric, i.e., \code{compute_gamma(v1, v2) != compute_gamma(v2, v1)}.
#'
#' @param v1,v2 Numeric vectors. Two vectors with \code{n} observations.
#' @param k Integer. The number of upper-order statistics used to
#' compute the gamma coefficient. Set by default to \code{k = floor(sqrt(n))}.
#' @return Numeric --- between 0 and 1. The gamma causal coefficient between \code{v1} and
#' \code{v2}.
compute_gamma <- function(v1, v2, k = floor(sqrt(n))){
  n <- NROW(v1)   # number of observations
  r1 <- rank(v1, ties.method = "first") # ranks of v1
  r2 <- rank(v2, ties.method = "first") # ranks of v2

  1 / (k * (n + 1)) * sum(r2[r1 > n - k + 1 / 2]) # gamma coefficient
}


#' Causal gamma matrix
#'
#' Computes the gamma coefficient matrix given a dataset \code{mat},
#' and a threshold \code{k}.
#'
#' @param mat Numeric matrix. Observation matrix with \code{n}
#' rows (observations) and \code{p} columns (variables).
#' @inheritParams compute_gamma
#' @return Numeric matrix. The entry \eqn{(i, j)} contains the
#' gamma coefficient between the \eqn{i}-th and the \eqn{j}-th
#' column of \code{mat}. If \eqn{i = j}, the value of the entry \eqn{(i, j)}
#' is set to \code{NA}.
compute_gamma_matrix <- function(mat, k = floor(sqrt(n))){
  n   <- NROW(mat) # number of observations
  p   <- NCOL(mat) # number of variables

  # compute gamma for all combinations of indices
  gamma_mat <- sapply(1:p,
              function(j){
                sapply(1:p,
                       function(i){
                         if (i == j){
                           NA
                         } else {
                           compute_gamma(mat[, i], mat[, j], k)
                         }
                       })
              })

  # return data
  return(gamma_mat)
}


#' Pairwise theoretical gamma coefficient
#'
#' Computes the theoretical gamma coefficient between node \code{i}
#' and node \code{j}, given the weighted adjacency matrix \code{wg},
#' the tail-index \code{alpha} and the scales of the
#' noise variables \code{noise_w}.
#'
#' @param wg Numeric matrix. The weighted adjacency matrix of a DAG.
#' All the entries must be non-negative.
#' @param i,j Integer. The indices of the nodes to consider to compute
#' the theoretical gamma coefficient.
#' @param alpha Numeric. The tail-index of the noise distribution
#' @param noise_w Numeric vector. A \eqn{p}-dimensional vector containing
#' the scaling coefficients of the noise variables. Note that \eqn{p} is the
#' number of variables (nodes) in the DAG. By default, all entries are
#' set equal to one.
#' @return Numeric --- between 0.5 and 1. The theoretical gamma coefficient
#' between node \code{i} and \code{j} in the DAG \code{wg}.
compute_gamma_theo <- function(wg, i, j, alpha, noise_w){

  # check that weighted adjacency matrix has no negative entries
  if (min(wg, na.rm = T) < 0){
    stop("The weighted adjacency matrix must have non-negative weights.")
  }

  # compute adjacency matrix
  g <- (wg != 0) * 1

  # compute theoretical gamma
  p <- NROW(wg)
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

  path_matrix <- t(get_all_paths(wg, "weighted"))
  p_ij <- (e %*% path_matrix ^ alpha %*% v_num) /
    (e %*% path_matrix ^ alpha %*% v_denom)
  gamma <- drop(1 - 0.5 * p_ij)

  return(gamma)
}


#' Theoretical gamma coefficient of a DAG
#'
#' Computes the theoretical gamma coefficient matrix given the DAG \code{g},
#' the tail-index \code{alpha} and the scales of the noise variables
#' \code{noise_w}.
#'
#' @inheritParams compute_gamma_theo
#' @return Numeric matrix. The entry \eqn{(i, j)} contains the theoretical
#' gamma coefficient between the \eqn{i}-th and the \eqn{j}-th node
#' of \code{g}. If \eqn{i = j}, the value of the entry \eqn{(i, j)} is set
#' to \code{NA}.
compute_gamma_theo_matrix <- function(g, alpha, noise_w){

  if (min(g, na.rm = T) < 0){
    stop("The weighted adjacency matrix must have non-negative weights.")
  }

  p   <- NCOL(g) # number of variables

  # compute theoretical gamma for all combinations of indices
  gamma_mat <- sapply(1:p,
                      function(j){
                        sapply(1:p,
                               function(i){
                                 if (i == j){
                                   NA
                                 } else {
                                   compute_gamma_theo(g, i, j, alpha, noise_w)
                                 }
                               })
                      })

  # return data
  return(gamma_mat)
}
