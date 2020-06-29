#' Causal order of a DAG
#'
#' Produces one causal order of the given DAG \code{dag}.
#' Copyright (c) 2013 Jonas Peters \email{peters@@math.ku.dk}.
#' All rights reserved.
#'
#' @param dag Square binary matrix. A matrix representing a DAG.
#' @return Numeric vector. The causal order of the DAG \code{dag}.
#' @noRd
compute_caus_order <- function(dag) {

  # check if dag is a DAG (i.e., a binary matrix)
  if (!all(dag %in% c(0, 1))){
    stop("The entries of dag must be either 0 or 1.")
  }

  p <- dim(dag)[2]
  remaining <- 1:p
  caus_order <- rep(NA, p)

  for (i in 1:(p - 1)) {
    root <- min(which(colSums(dag) == 0))
    caus_order[i] <- remaining[root]
    remaining <- remaining[-root]
    dag <- dag[-root, -root]
  }
  caus_order[p] <- remaining[1]
  return(caus_order)
}


#' Check the causal order against a DAG
#'
#' Produces \code{TRUE} if \code{caus_order} is a causal order of the DAG
#' \code{dag}, \code{FALSE} otherwise.
#'
#' @param caus_order Numeric vector.
#' Represents a causal order.
#' @inheritParams compute_caus_order
#' @return Boolean.
#' @noRd
check_caus_order <- function(caus_order, dag) {

  # check if dag is a DAG (i.e., a binary matrix)
  if (!all(dag %in% c(0, 1))){
    stop("The entries of dag must be either 0 or 1.")
  }

  # check if caus_order has some NA
  if (any(is.na(caus_order))){
    stop("The causal order caus_order cannot contain NA.")
  }

  p <- NROW(dag)

  for (i in 1:p) {
    current_node <- caus_order[i]
    indegree <- sum(dag[, current_node])

    if (indegree > 0) {
      return(FALSE)
    }

    dag[current_node, ] <- 0
  }

  return(TRUE)
}


#' Get the ancestors of a DAG
#'
#' Produces a matrix where the entry \eqn{(i, j) = 1}
#' if \eqn{i} is an ancestor of \eqn{j} in the DAG \code{dag},
#' and \eqn{0} otherwise.
#'
#' @inheritParams compute_caus_order
#' @return Square binary matrix. The ancestral matrix of the DAG \code{dag}.
#' @noRd
get_ancestors <- function(dag){

  # check if dag is a DAG (i.e., a binary matrix)
  if (!all(dag %in% c(0, 1))){
    stop("The entries of dag must be either 0 or 1.")
  }

  p <- NROW(dag)

  # compute ancestors
  ident_mat <- diag(1, p)
  ancestors <- ident_mat + dag

  for (i in 1:p){
    ancestors <- ident_mat + dag %*% ancestors
  }

  return( (ancestors != 0) * 1)
}


#' Get the descendants of a DAG
#'
#' Produces a matrix where the entry \eqn{(i, j) = 1}
#' if \eqn{i} is a descendant of \eqn{j} in the DAG \code{dag},
#' and \eqn{0} otherwise.
#'
#' @inheritParams compute_caus_order
#' @return Square binary matrix. The descendant matrix of
#' the DAG \code{dag}.
#' @noRd
get_descendants <- function(dag){

  # check if dag is a DAG (i.e., a binary matrix)
  if (!all(dag %in% c(0, 1))){
    stop("The entries of dag must be either 0 or 1.")
  }

  # return descendants
  descendants <- t(get_ancestors(dag))
  (descendants != 0) * 1
}


#' Get the parents of a DAG
#'
#' Produces a matrix where the entry \eqn{(i, j) = 1} if \eqn{i} is a parent
#' of \eqn{j} in the DAG \code{dag}, and \eqn{0} otherwise.
#'
#' @inheritParams compute_caus_order
#' @return Square binary matrix. The parental matrix of
#' the DAG \code{dag}, i.e., the DAG itself.
#' @noRd
get_parents <- function(dag){

  # check if dag is a DAG (i.e., a binary matrix)
  if (!all(dag %in% c(0, 1))){
    stop("The entries of dag must be either 0 or 1.")
  }

  # return parents
  dag
}


#' Get the children of a DAG
#'
#' Produces a matrix where the entry \eqn{(i, j) = 1} if \eqn{i} is a child
#' of \eqn{j} in the DAG \code{dag}, and \eqn{0} otherwise.
#'
#' @inheritParams compute_caus_order
#' @return Square binary matrix. The children matrix of
#' the DAG \code{dag}.
#' @noRd
get_children <- function(dag){

  # check if dag is a DAG (i.e., a binary matrix)
  if (!all(dag %in% c(0, 1))){
    stop("The entries of dag must be either 0 or 1.")
  }

  # return children
  t(dag)
}


#' Get all the directed paths of a DAG
#'
#' Given the adjacency matrix \code{adj_mat}, produces a matrix
#' where the entry \eqn{(i, j)} measures the total weight of the paths
#' from node \eqn{i} to node \eqn{j}.
#' \strong{Note}: if \eqn{i = j}, the entry \eqn{(i, j) = 1}.
#'
#' @param adj_mat Square numeric matrix. The adjacency matrix of a DAG.
#' @return Square numeric matrix. The path matrix associated to \code{adj_mat}.
#' @noRd
get_all_paths <- function(adj_mat){

  p <- NROW(adj_mat)

  # compute paths
  ident_mat <- diag(1, p)
  paths <- ident_mat + adj_mat

  for (i in 1:p){
    paths <- ident_mat + adj_mat %*% paths
  }

  return(paths)
}


#' Convert a causal order into a fully connected DAG
#'
#' Convert the given causal order \code{caus_order} into a fully connected
#' DAG.
#'
#' @inheritParams check_caus_order
#' @return Square binary matrix. A fully connected DAG that agrees
#' with the given order.
#' @noRd
caus_order_to_dag <- function(caus_order){

  # check if caus_order has some NA
  if (any(is.na(caus_order))){
    stop("The causal order caus_order cannot contain NA.")
  }

  p <- length(caus_order)
  dag <- upper.tri(x = matrix(0, nrow = p, ncol = p)) * 1
  inv_caus_order <- order(caus_order)
  return(dag[inv_caus_order, inv_caus_order])
}


#' Convert a DAG into a CPDAG
#'
#' Convert the given DAG \code{dag} into a CPDAG
#' DAG.
#'
#' @inheritParams compute_caus_order
#' @return Square binary matrix. A CPDAG that agrees
#' with the given DAG.
#' @noRd
dag_to_cpdag <- function(dag){

  # check if dag is a DAG (i.e., a binary matrix)
  if (!all(dag %in% c(0, 1))){
    stop("The entries of dag must be either 0 or 1.")
  }

  # return cpdag
  (pcalg::dag2cpdag(dag)) * 1
}


#' Convert a causal order into a CPDAG
#'
#' Convert the given causal order \code{caus_order} into a CPDAG
#' DAG.
#'
#' @inheritParams check_caus_order
#' @return Square binary matrix. A CPDAG that agrees
#' with the given DAG.
#' @noRd
caus_order_to_cpdag <- function(caus_order){

  # check if caus_order has some NA
  if (any(is.na(caus_order))){
    stop("The causal order caus_order cannot contain NA.")
  }

  # return cpdag
  dag <- caus_order_to_dag(caus_order)
  (pcalg::dag2cpdag(dag)) * 1
}


#' Compute ancestral distance
#'
#' Computes the ancestral distance between the DAG \code{dag} and
#' the causal order \code{caus_order}. The ancestral distance is defined
#' as the number of pairs \eqn{(i, j)} where \eqn{i} is an ancestor
#' of \eqn{j} and \eqn{i} is placed after \eqn{j} in the causal order
#' \code{caus_order}.
#'
#' @inheritParams check_caus_order
#' @return Numeric --- between 0 and 1. The ancestral distance between the
#' DAG \code{dag} and the causal order \code{caus_order}.
#' @noRd
compute_ancestral_distance <- function(dag, caus_order){

  # check if dag is a DAG (i.e., a binary matrix)
  if (!all(dag %in% c(0, 1))){
    stop("The entries of dag must be either 0 or 1.")
  }

  # check if caus_order contains some NA
  if (any(is.na(caus_order))){
    stop("The causal order caus_order cannot contain NA.")
  }

  p <- NROW(dag)
  max_inversion <- p * (p - 1) / 2
  ancestors <- get_ancestors(dag)
  ancestors_ord <- ancestors[caus_order, caus_order]
  return(sum(ancestors_ord[lower.tri(ancestors_ord)]) / max_inversion)
}


#' Compute structural intervention distance
#'
#' Compute the structural intervention distance between the DAG \code{dag}
#' and the estimated DAG (or CPDAG) \code{est_g}.
#' The structural intervention distance is defined as in the paper from
#' Peters J., and Bühlmann P.,
#' \url{https://www.mitpressjournals.org/doi/full/10.1162/NECO_a_00708}.
#' In general, the structural intervention distance is
#' not symmetric, i.e., \code{compute_str_int_distance(dag, est_g) !=}
#' \code{compute_str_int_distance(est_g, dag)}.
#'
#' @inheritParams compute_caus_order
#' @param est_g Square binary matrix. The estimated DAG or CPDAG.
#' @return  Numeric --- between 0 and 1. The structural intervention
#' distance between the DAG \code{dag} and the DAG (or CPDAG) \code{est_g}.
#' @noRd
compute_str_int_distance <- function(dag, est_g){

  # check if dag is a DAG (i.e., a binary matrix)
  if (!all(dag %in% c(0, 1))){
    stop("The entries of dag must be either 0 or 1.")
  }

  # check if est_g is a DAG or CPDAG (i.e., a binary matrix)
  if (!all(est_g %in% c(0, 1))){
    stop("The entries of est_g must be either 0 or 1.")
  }

  p <- NROW(dag)
  s <- SID::structIntervDist(dag, est_g)
  return(s$sidLowerBound / (p * (p - 1)))
}


#' Compute structural Hamming distance
#'
#' Compute the structural Hamming distance between the true CPDAG
#' \code{cpdag} and the estimated CPDAG \code{est_cpdag}.
#' The structural Hamming distance is defined as in the paper from
#' Tsamardinos I., Brown L.E., and Aliferis C.F.,
#' \url{https://link.springer.com/article/10.1007/s10994-006-6889-7}.
#'
#' @param cpdag Square binary matrix. A matrix representing a CPDAG.
#' @param est_cpdag Square binary matrix. The estimated CPDAG.
#' @return Numeric --- between 0 and 1. The structural Hamming
#' distance between a true CPDAG \code{cpdag} and the
#' estimated CPDAG \code{est_cpdag}.
#' @noRd
compute_str_ham_distance <- function(cpdag, est_cpdag){

  # check if cpdag is a CPDAG (i.e., a binary matrix)
  if (!all(cpdag %in% c(0, 1))){
    stop("The entries of cpdag must be either 0 or 1.")
  }

  # check if est_cpdag is a CPDAG (i.e., a binary matrix)
  if (!all(est_cpdag %in% c(0, 1))){
    stop("The entries of est_cpdag must be either 0 or 1.")
  }

  p <- NROW(cpdag)
  s <- SID::hammingDist(cpdag, est_cpdag)
  return(s / (p * (p - 1) / 2))
}
