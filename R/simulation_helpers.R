#' Simulate random DAG
#'
#' Simulates a directed acyclic graph (DAG) and returns its adjacency matrix.
#' Copyright (c) 2013 Jonas Peters \email{peters@@math.ku.dk}.
#' All rights reserved.
#'
#' @param p Integer --- greater than 2. Number of nodes.
#' @param prob_connect Numeric --- between 0 and 1. The probability that an edge
#' \eqn{i {\rightarrow} j} is added to the DAG.
#' @param caus_order Numeric vector. The causal order of the DAG.
#' If the argument is not provided it is generated randomly.
#' @param sparse Logical. Should the adjacency matrix be sparse
#' (as a data structure)?
#'
#' @return Numeric matrix. The adjacency matrix of the simulated DAG.
random_dag <- function(p, prob_connect,
                       caus_order = sample(p, p, replace = FALSE),
                       sparse = FALSE){

  if (p < 3){
    stop("The number of nodes p must be greater than 2.")
  }

  if (!missing(caus_order)){
    if (length(caus_order) != p){
      stop(paste("The number of nodes p does not match with the length of",
                 "the causal order."))
    }
  }

  if (!missing(caus_order)){
    caus_order <- caus_order[p:1]
  }

  if (sparse){
    dag <- Matrix::Matrix(nrow = p, ncol = p, 0, sparse = TRUE)
  } else {
    dag <- matrix(nrow = p, ncol = p, 0)
  }

  for (i in 1:(p - 2)) {
    node <- caus_order[i]
    possible_parents <- caus_order[(i + 1):p]
    number_parents <- rbinom(n = 1, size = (p - i), prob = prob_connect)
    parents <- sample(x = possible_parents, size = number_parents,
                      replace = FALSE)
    dag[parents, node] <- rep(1, number_parents)
  }

  # Sample does not work properly when choosing from sets with one element.
  # We thus consider the last case separately.
  node <- caus_order[p - 1]
  parent_yes_no <- rbinom(n = 1, size = 1, prob = prob_connect)
  dag[caus_order[p], node] <- parent_yes_no

  return(dag)
}


#' Sample random coefficients
#'
#' Sample random coefficients from uniform distribution
#' for the given DAG \code{g}.
#' Copyright (c) 2013 Jonas Peters \email{peters@@math.ku.dk}.
#' All rights reserved.
#'
#' @inheritParams compute_caus_order
#' @param lB Numeric. The lower bound for the uniform distribution.
#' @param uB Numeric. The upper bound for the uniform distribution.
#' It must be stricly greater than \code{lB}.
#' @param two_intervals Logical. Should the coefficient be sampled
#' from two symmetric uniform distributions?
#' If \code{two_intervals == TRUE}, \code{lB} and \code{uB} must be
#' positive.
#' @return Numeric matrix. The weighted adjacency matrix of the
#' DAG \code{g}.
random_b <- function(g, lB = 0.1, uB = 0.9, two_intervals = FALSE){

  if (!(uB > lB)) {
    stop("The upper bound uB must be stricly greater than the lower bound lB")
  }

  if (two_intervals == TRUE){
    if (!(lB > 0 & uB > 0)){
      stop(paste("If two_intervals == TRUE, lB and uB must be positive."))
    }
  }

  numCoeff <- sum(g)
  B <- t(g)
  if (numCoeff == 1) {
    coeffs <- sample(c(-1, 1), size = numCoeff,
                     TRUE) ^ (two_intervals) * runif(1, lB, uB)
  } else {
    coeffs <- diag(sample(c(-1, 1), size = numCoeff,
                          TRUE) ^ (two_intervals)) %*% runif(numCoeff, lB, uB)
  }
  B[B == 1] <- coeffs
  return(B)
}
