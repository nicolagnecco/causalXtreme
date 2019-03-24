#' Randomly pick elements
#' Select each element of \code{vec} with probability \code{prob} and
#' produce a vector with the selected elements.
#'
#' @param vec Numeric vector. Vector containing the elements to select.
#' @param prob Numeric --- between 0 and 1. The probability of selecting
#' each element of \code{vec}.
#'
#' @return Numeric vector. A vector with the selected elements
#'
pick_elements <- function(vec, prob){
  r <- rbinom(n = length(vec), size = 1, prob = prob)
  vec[r == 1]
}



#' Inverse mirror uniform
#'
#' Produces the quantile of a mirrored uniform distribution
#' associated to the probability \code{prob}. The mirrored uniform
#' distribution has support [-max, -min] \eqn{\cup} [min, max].
#'
#' @param prob Numeric --- between 0 and 1. The probability associated
#' to the desired quantile.
#' @param min Numeric --- above 0. The lower bound of the positive half of
#' the support.
#' @param max Numeric --- above 0. The upper bound of the positive half of
#' the support. Note that \code{max} must be strictly greater than
#' \code{min}.
#'
#' @return Numeric --- between \code{-max} and \code{max}. The quantile
#' associated to the probability \code{prob}.
inverse_mirror_uniform <- function(prob, min, max){

  # check if max < min
  if (!(max > min)){
    stop(paste("The maximum value, max, must be strictly greater",
               "than the minimum value, min."))
  }

  # check if min > 0 and max > 0
  if (!(min > 0 & max > 0)){
    stop(paste("Both the minimum and the maximum values, min and max,",
               "must be positive."))
  }

  if (prob == 0){
    return(-max)
  } else if (prob > 0 & prob < 1 / 2) {
    return(2 * prob * (max - min) - max)
  } else if (prob >= 1 / 2 & prob < 1){
    return( (2 * prob - 1) * (max - min) + min)
  } else if (prob == 1){
    return(max)
  } else {
    stop("The probability, prob, must be between 0 and 1.")
  }
}



#' Sample from uniform family
#'
#' Sample n elements from a uniform distribution with lower and upper
#' bound equal to \code{min} and \code{max}, respectively.
#' If \code{mirror == TRUE}, the elements are sampled from a mirrored
#' uniform distribution, see \code{\link{inverse_mirror_uniform}}.
#'
#' @param n Integer. The number of elements to sample.
#' @param min Numeric. The lower bound of the distribution support.
#' If \code{mirror == TRUE}, this represents the lower bound of
#' the positive half of the support.
#' @param max Numeric. The upper bound of the distribution support.
#' If \code{mirror == TRUE}, this represents the upper bound of
#' the positive half of the support.
#' Note that \code{max} must be strictly greater than \code{min}.
#' @param mirror Logical. Should the elements be sampled from a
#' mirrored uniform distribution? If \code{mirror == TRUE},
#' \code{min} and \code{max} must be strictly greater than 0.
#'
#' @return Numeric vector. A vector with the sampled elements.
sample_uniform <- function(n, min, max, mirror = FALSE){

  # check if max < min
  if (!(max > min)){
    stop(paste("The maximum value, max, must be strictly greater",
               "than the minimum value, min."))
  }

  if (mirror == TRUE){

    purrr::map_dbl(runif(n), inverse_mirror_uniform, min = min, max = max)

  } else {

    runif(n = n, min = min, max = max)
  }
}



#' Simulate random DAG
#'
#' Simulates a directed acyclic graph (DAG) and returns its adjacency matrix.
#' Copyright (c) 2013 Jonas Peters \email{peters@@math.ku.dk}.
#' All rights reserved.
#'
#' @param p Integer --- greater than 0. Number of nodes.
#' @param prob Numeric --- between 0 and 1. The probability that an edge
#' \eqn{i {\rightarrow} j} is added to the DAG.
#' @param caus_order Numeric vector. The causal order of the DAG.
#' If the argument is not provided it is generated randomly.
#'
#' @return Numeric matrix. The adjacency matrix of the simulated DAG.
random_dag <- function(p, prob, caus_order = sample(p, p, replace = FALSE)){

  # check inputs
  if (p < 1){
    stop("The number of nodes p must be greater than 0.")
  }

  if (!missing(caus_order)){
    if (length(caus_order) != p){
      stop(paste("The number of nodes p does not match with the length of",
                 "the causal order."))
    }
  }

  # get random dag
  dag <- matrix(0, nrow = p, ncol = p)

  if (p > 1){
    for (i in 1:(p - 1)){
      elms <- pick_elements(caus_order[(i + 1):p], prob)
      dag[caus_order[i], elms] <- 1
    }
  }
  dag
}



#' Sample random coefficients
#'
#' Sample random coefficients from uniform distribution
#' for the given DAG \code{g}.
#' Copyright (c) 2013 Jonas Peters \email{peters@@math.ku.dk}.
#' All rights reserved.
#'
#' @inheritParams compute_caus_order
#' @param lb Numeric. The lower bound of the coefficient that can be
#' sampled.
#' @param ub Numeric. The lower bound of the coefficient that can be
#' sampled. It must be stricly greater than \code{lb}.
#' @param two_intervals Logical. Should the coefficient be sampled
#' from two symmetric uniform distributions?
#' If \code{two_intervals == TRUE}, \code{lb} and \code{ub} must be
#' positive.
#' @return Numeric matrix. The weighted adjacency matrix of the
#' DAG \code{g}.
random_coeff <- function(g, lb = 0.1, ub = 0.9, two_intervals = FALSE){

  # check if g is a (non-weighted) adjacency matrix
  if (!all(g %in% c(0, 1))){
    stop("The entries of g must be either 0 or 1.")
  }

  g_coeff <- matrix(0, nrow = NROW(g), ncol = NCOL(g))
  num_coeff <- sum(g)
  g_coeff[g == 1] <- sample_uniform(num_coeff, lb, ub, two_intervals)
  g_coeff
}
