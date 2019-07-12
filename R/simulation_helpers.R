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
#' @param prob_connect Numeric --- between 0 and 1. The probability that an edge
#' \eqn{i {\rightarrow} j} is added to the DAG.
#' @param caus_order Numeric vector. The causal order of the DAG.
#' If the argument is not provided it is generated randomly.
#'
#' @return Square binary matrix. A matrix representing the random DAG.
random_dag <- function(p, prob_connect,
                       caus_order = sample(p, p, replace = FALSE)){

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
      elms <- pick_elements(caus_order[(i + 1):p], prob_connect)
      dag[caus_order[i], elms] <- 1
    }
  }
  dag
}


#' Sample random coefficients
#'
#' Sample random coefficients from uniform distribution
#' for the given DAG \code{dag}.
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
#' @return Square numeric matrix. The adjacency matrix of the underlying
#' DAG \code{dag}.
random_coeff <- function(dag, lb = 0.1, ub = 0.9, two_intervals = TRUE){

  # check if dag is a (non-weighted) adjacency matrix
  if (!all(dag %in% c(0, 1))){
    stop("The entries of dag must be either 0 or 1.")
  }

  adj_mat <- matrix(0, nrow = NROW(dag), ncol = NCOL(dag))
  num_coeff <- sum(dag)
  adj_mat[dag == 1] <- sample_uniform(num_coeff, lb, ub, two_intervals)
  adj_mat
}


#' Add confounders to a DAG
#'
#' Add confounders (i.e., hidden variables) at random to the DAG \code{dag}.
#' Each confounder is a source node and acts only on two observed variables.
#' The probability that a pair of observed variables is confounded is set by
#' \code{prob_confound}.
#'
#' @inheritParams compute_caus_order
#' @param prob_confound Numeric --- between 0 and 1. The probability that
#' a pair of observed nodes is confounded.
#' @return List. The list is made of:
#' \itemize{
#' \item \code{dag_confounders} --- Square binary matrix. Represents the full
#' DAG made of observed and hidden variables.
#' \item \code{pos_confounders} --- Integer vector. Represents the position
#' of confounders (rows and columns) in dag_confounders.
#' }
#'
add_random_confounders <- function(dag, prob_confound){

  # check if dag is a (non-weighted) adjacency matrix
  if (!all(dag %in% c(0, 1))){
    stop("The entries of dag must be either 0 or 1.")
  }

  # Consider all possible pairs of nodes in the DAG
  p <- NROW(dag)
  node_pairs <- which(upper.tri(dag), arr.ind = TRUE)
  max_n_confounders  <- NROW(node_pairs)

  # Choose the confounders randomly
  n_confounders  <- rbinom(n = 1, size = max_n_confounders,
                           prob = prob_confound)
  confounders <- sample(x = 1:max_n_confounders, size = n_confounders,
                        replace = FALSE)

  # Exit function if there are no confounders
  if (n_confounders == 0){
    return(list (dag_confounders = dag,
                 pos_confounders = integer(0)))
  }

  # Add confounders to the DAG
  dag_confounders <- dag

  for (i in 1:n_confounders){
    confounded.pair <- node_pairs[confounders[i], ]
    col_to_append <- rep(0, NROW(dag_confounders))
    row_to_append <- rep(0, NROW(dag_confounders) + 1)
    row_to_append[confounded.pair] <- 1
    dag_confounders <- cbind(dag_confounders, col_to_append)
    dag_confounders <- rbind(dag_confounders, row_to_append)
  }

  dag_confounders <- unname(dag_confounders)
  pos_confounders <-  (p + 1):NROW(dag_confounders)

  # Return data
  return(list (dag_confounders = dag_confounders,
               pos_confounders = pos_confounders))
}


#' Simulate noise observations
#'
#' Sample \code{n} observations for \code{p} independent noise variables
#' from the distribution \code{distr}. The distribution is one of:
#' \itemize{
#' \item \code{student_t}, in this case the user has to specify the
#' \code{tail_index}, i.e., the degrees of freedom,
#' \item \code{gaussian},
#' \item \code {log_normal}.
#' }
#' @param n Positive integer. The number of observations.
#' @param p Positive integer. The number of variables.
#' @param distr Character. The distribution of the noise.
#' @param tail_index Positive numeric. The tail index, i.e., degrees
#' of freedom, of the noise.
#' @return Numeric matrix. Dataset matrix with \code{n}
#' rows (observations) and \code{p} columns (variables).
simulate_noise <- function(n, p, distr = c("student_t", "gaussian",
                                          "log_normal")[1], tail_index){
  switch(distr,
         "student_t" = {

           noise <- array(rt(n * p, df = tail_index), dim = c(n, p))

         },
         "gaussian" = {

           noise <- array(rnorm(n * p), dim = c(n, p))
         },
         "log_normal" = {

           noise <- array(rlnorm(n * p), dim = c(n, p))

         },
         stop("Wrong distribution. Enter one of 'student_t',
              'gaussian', 'log_normal'."))

  return(noise)
}

#' Generate data from non-linear Structural Equation Model
#'
#' Generate data from non-linear Structural Causal Model as shown in the
#' paper "Causality in heavy-tailed models".
#'
#' @inheritParams get_all_paths
#' @param noise Numeric matrix. Dataset matrix with \code{n}
#' rows (observations) and \code{p} columns (variables).
#' @return Numeric matrix. Dataset matrix with \code{n}
#' rows (observations) and \code{p} columns (variables).
nonlinear_scm <- function(adj_mat, noise){

  n <- NROW(noise)
  p <- NROW(adj_mat)
  dag <- (adj_mat != 0) * 1
  nodes <- compute_caus_order(dag)

  dataset <- matrix(0, nrow = n, ncol = p)
  dataset_transf <- matrix(0, nrow = n, ncol = p)

  for (i in nodes){
    betas <- adj_mat[, i]
    eps <- noise[, i]
    dataset[, i] <- dataset_transf %*% betas + eps
    dataset_transf[, i] <- broken_hockeystick(dataset[, i])
  }

  return(dataset)
}


#' Apply shifted-hockeystick function
#'
#' Apply shifted-hockeystick function to a numeric vector \code{v}.
#' The shifted-hockeystick function sets to zero all entries of
#' \code{v} greater or equal to the \code{q_low}-quantile and
#' strictly less than the \code{q_high}-quantile.
#'
#' @param v Numeric vector.
#' @param q_low Numeric, between 0 and 1.
#' @param q_high Numeric, between 0 and 1.
#' @return Numeric vector.
broken_hockeystick <- function(v, q_low = 0.01, q_high = 0.99){
  n <- length(v)
  r <- rank(v, ties.method = "first")
  ind <- which(r > floor(n * q_low) & r <= ceiling(n * q_high))
  v[ind] <- 0
  return(v)
}


#' Transform to uniform margins
#'
#' Transform the numeric vector \code{x} to uniform margin
#' between 0 and 1.
#'
#' @inheritParams broken_hockeystick
#' @return Numeric vector.
uniform_margin <- function(v){
  n <- length(v)
  rank(v) / n
}
