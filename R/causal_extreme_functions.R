#' Causal gamma coefficient
#'
#' Computes the causal gamma coefficient between \code{v1} and \code{v2},
#' given the threshold \code{k}. In general, the gamma coefficient is not
#' symmetric, i.e., \code{compute_gamma(v1, v2) != compute_gamma(v2, v1)}.
#'
#' @param v1,v2 Numeric vectors. Two vectors with \code{n} observations.
#' @param k Integer. The number of upper-order statistics used to
#' compute the gamma coefficient. Set by default to \code{k = floor(sqrt(n))}.
#' @return Numeric. The gamma causal coefficient between \code{v1} and
#' \code{v2}, which lies between 0.5 and 1.
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
#' @param mat Numeric matrix. Matrix with \code{n} rows (observations)
#' and \code{p} columns (variables).
#' @inheritParams compute_gamma
#' @return Numeric matrix. The entry \eqn{(i, j)} contains the
#' gamma coefficient between the \eqn{i}-th and the \eqn{j}-th
#' column of \code{mat}. If \eqn{i = j}, the entry \eqn{(i, j)} = \code{NA}.
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
