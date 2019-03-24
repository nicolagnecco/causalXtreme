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
#' @param lb Numeric. The lower bound for the uniform distribution.
#' @param ub Numeric. The upper bound for the uniform distribution.
#' It must be stricly greater than \code{lb}.
#' @param two_intervals Logical. Should the coefficient be sampled
#' from two symmetric uniform distributions?
#' If \code{two_intervals == TRUE}, \code{lb} and \code{ub} must be
#' positive.
#' @return Numeric matrix. The weighted adjacency matrix of the
#' DAG \code{g}.
random_b <- function(g, lb = 0.1, ub = 0.9, two_intervals = FALSE){

  # check if g is a (non-weighted) adjacency matrix
  if (!all(g %in% c(0, 1))){
    stop("The entries of g must be either 0 or 1.")
  }

  # check if the lower and upper bound are well specified
  if (!(ub > lb)) {
    stop("The upper bound ub must be stricly greater than the lower bound lb")
  }

  if (two_intervals == TRUE){
    if (!(lb > 0 & ub > 0)){
      stop(paste("If two_intervals == TRUE, lb and ub must be positive."))
    }
  }

  numCoeff <- sum(g)
  if (numCoeff == 1) {
    coeffs <- sample(c(-1, 1), size = numCoeff,
                     TRUE) ^ (two_intervals) * runif(1, lb, ub)
  } else {
    coeffs <- diag(sample(c(-1, 1), size = numCoeff,
                          TRUE) ^ (two_intervals)) %*% runif(numCoeff, lb, ub)
  }
  g[g == 1] <- coeffs
  return(g)
}


addRandomConfounders <- function(N, probConfouder, conf_w, alpha){
  ## numeric_matrix numeric numeric df -> numeric_matrix
  ## produces a noise matrix where each pair of variables (i, j) in N
  ## is confounded with probability probConfounder;
  ## conf_w is the signal proportion of the confounder versus the noise;
  ## alpha is the tail index of the distribution
  ## ASSUME:
  ## 1. conf_w is between 0 and 1

  # Check inputs
  if(conf_w > 1 | conf_w < 0){stop("The signal weight must be between 0 and 1.")}

  # Get variables
  n <- NROW(N)
  p <- NCOL(N)

  # Consider all possible pairs of nodes in G
  M <- which(upper.tri(diag(p)), arr.ind=TRUE)
  nPossibleConf  <- NROW(M)

  # Choose the confounders randomly
  numberConf  <- rbinom(n=1, size=nPossibleConf, prob=probConfounder)
  Confounders <- sample(x = 1:nPossibleConf, size = numberConf, replace = FALSE)
  Confounders <- M[Confounders, ]

  if(numberConf == 0){
    return(N)
  }

  # Sample the confounder weights
  S <- matrix(0, nrow = numberConf, ncol = p) # preallocate matrix with confounder-variable pairs
  confvar_pairs <- cbind(rep(1:numberConf, each=2), c(t(Confounders))) # confounder-variable pairs
  S[confvar_pairs] <- 1 # populate matrix with confounder-variable pairs
  S <- abs(randomB(t(S), lB = lB, uB = uB)) # sample random weights of confounders

  # Rescale confounder weights
  S <- scale_adjmat(S, conf_w, alpha)
  noise_w <- 1 - conf_w
  N <- scale_noise(N, S, noise_w, alpha)

  # Generate confounders
  C <- simulateNoise(n, numberConf, distr, alpha)

  # Return confounded noise
  return(N + C %*% S)
}

