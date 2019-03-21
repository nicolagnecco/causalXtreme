#' Causal order of a DAG
#'
#' Produces one causal order of the DAG, represented as a non-weighted
#' adjacency matrix \code{g}.
#' Copyright (c) 2013 Jonas Peters \email{peters@@math.ku.dk}.
#' All rights reserved.
#'
#' @param g Numeric matrix. The adjacency matrix of a DAG.
#' @return Numeric vector. The causal order of the DAG \code{g}.
compute_caus_order <- function(g) {

  # check if g is a (non-weighted) adjacency matrix
  if (!all(g %in% c(0, 1))){
    stop("The entries of g must be either 0 or 1.")
  }

  p <- dim(g)[2]
  remaining <- 1:p
  caus_order <- rep(NA, p)

  for (i in 1:(p - 1)) {
    root <- min(which(colSums(g) == 0))
    caus_order[i] <- remaining[root]
    remaining <- remaining[-root]
    g <- g[-root, -root]
  }
  caus_order[p] <- remaining[1]
  return(caus_order)
}


#' Check the causal order against a DAG
#'
#' Produces \code{TRUE} if \code{caus_order} is a causal order of the DAG
#' \code{g}, \code{FALSE} otherwise.
#'
#' @param caus_order Numeric vector.
#' Represents a causal order.
#' @inheritParams compute_caus_order
#' @return Logical.
check_caus_order <- function(caus_order, g) {

  # check if g is a (non-weighted) adjacency matrix
  if (!all(g %in% c(0, 1))){
    stop("The entries of g must be either 0 or 1.")
  }

  # check if caus_order has some NA
  if (any(is.na(caus_order))){
    stop("The causal order caus_order cannot contain NA.")
  }

  p <- NROW(g)

  for (i in 1:p) {
    current_node <- caus_order[i]
    indegree <- sum(g[, current_node])

    if (indegree > 0) {
      return(FALSE)
    }

    g[current_node, ] <- 0
  }

  return(TRUE)
}


#' Get the ancestors of a DAG
#'
#' Produces a matrix where the entry \eqn{(i, j) = 1}
#' if \eqn{i} is an ancestor of \eqn{j} in the DAG \code{g},
#' and \eqn{0} otherwise.
#'
#' @inheritParams compute_caus_order
#' @return Numeric matrix. The ancestral matrix of the DAG \code{g}.
get_ancestors <- function(g){

  # check if g is a (non-weighted) adjacency matrix
  if (!all(g %in% c(0, 1))){
    stop("The entries of g must be either 0 or 1.")
  }

  p <- NROW(g)

  # return ancestors
  ancestors <- solve(diag(p) - g)
  (ancestors != 0) * 1
}


#' Get the descendants of a DAG
#'
#' Produces a matrix where the entry \eqn{(i, j) = 1}
#' if \eqn{i} is a descendant of \eqn{j} in the DAG \code{g},
#' and \eqn{0} otherwise.
#'
#' @inheritParams compute_caus_order
#' @return Numeric matrix. The descendant matrix of
#' the DAG \code{g}.
get_descendants <- function(g){

  # check if g is a (non-weighted) adjacency matrix
  if (!all(g %in% c(0, 1))){
    stop("The entries of g must be either 0 or 1.")
  }

  p <- NROW(g)

  # return descendants
  descendants <- t(solve(diag(p) - g))
  (descendants != 0) * 1
}


#' Get the parents of a DAG
#'
#' Produces a matrix where the entry \eqn{(i, j) = 1} if \eqn{i} is a parent
#' of \eqn{j} in the DAG \code{g}, and \eqn{0} otherwise.
#'
#' @inheritParams compute_caus_order
#' @return Numeric matrix. The parental (adjacency) matrix of
#' the DAG \code{g}.
get_parents <- function(g){

  # check if g is a (non-weighted) adjacency matrix
  if (!all(g %in% c(0, 1))){
    stop("The entries of g must be either 0 or 1.")
  }

  # return parents
  g
}


#' Get the children of a DAG
#'
#' Produces a matrix where the entry \eqn{(i, j) = 1} if \eqn{i} is a child
#' of \eqn{j} in the DAG \code{g}, and \eqn{0} otherwise.
#'
#' @inheritParams compute_caus_order
#' @return Numeric matrix. The children matrix of
#' the DAG \code{g}.
get_children <- function(g){

  # check if g is a (non-weighted) adjacency matrix
  if (!all(g %in% c(0, 1))){
    stop("The entries of g must be either 0 or 1.")
  }

  # return children
  t(g)
}


#' Get all the directed paths of a DAG
#'
#' Given the DAG \code{g}, produces a matrix where the entry \eqn{(i, j)}
#' \itemize{
#' \item counts the number of directed paths from node \eqn{i} to node \eqn{j},
#' if \code{type == "count"},
#' \item measures the total weight of the paths from node \eqn{i}
#' to node \eqn{j}, if \code{type == "weighted"}.
#' }
#' \strong{Note}: if \eqn{i = j}, the entry \eqn{(i, j) = 1}.
#'
#' @param wg Numeric matrix. The (possibly weighted) adjacency matrix of a DAG.
#' @param type String. Is one of:
#' \itemize{
#' \item \code{"count"} (default),
#' \item \code{"weighted"}.
#' }
#' @return Numeric matrix. The path matrix of the DAG \code{g}.
get_all_paths <- function(wg, type = c("count", "weighted")[1]){
  switch(type,
         "count" = {
           wg <- (wg != 0) * 1
         },
         "weighted" = {
           wg <- wg
         },
         stop("Wrong type. Please enter one of the following:
              'count', 'weighted'."))

  p <- NROW(wg)
  solve(diag(p) - wg)
}


#' Compute ancestral distance
#'
#' Computes the ancestral distance between the DAG \code{g} and
#' the causal order \code{caus_order}. The ancestral distance is defined
#' as the number of pairs \eqn{(i, j)} where \eqn{i} is an ancestor
#' of \eqn{j} and \eqn{i} is placed after \eqn{j} in the causal order
#' \code{caus_order}.
#'
#' @inheritParams check_caus_order
#' @return Numeric --- between 0 and 1. The ancestral distance between the
#' DAG \code{g} and the causal order \code{caus_order}.
compute_ancestral_distance <- function(g, caus_order){

  # check if g is a (non-weighted) adjacency matrix
  if (!all(g %in% c(0, 1))){
    stop("The entries of g must be either 0 or 1.")
  }

  # check if caus_order contains some NA
  if (any(is.na(caus_order))){
    stop("The causal order caus_order cannot contain NA.")
  }

  p <- NROW(g)
  max_inversion <- p * (p - 1) / 2
  ancestors <- get_ancestors(g)
  ancestors_ord <- ancestors[caus_order, caus_order]
  return(sum(ancestors_ord[lower.tri(ancestors_ord)]) / max_inversion)
}


#' Compute structural intervention distance
#'
#' Compute the structural intervention distance between the DAG \code{g} and
#' the estimated DAG \code{est_g}. The structural intervention distance is
#' defined as in the paper from Peters J. and Bühlmann P.,
#' \url{https://www.mitpressjournals.org/doi/full/10.1162/NECO_a_00708}.
#' In general, the structural intervention distance is
#' not symmetric, i.e., \code{compute_str_int_distance(g, est_g) !=}
#' \code{compute_str_int_distance(est_g, g)}. Also, \code{est_g} can be
#' a CPDAG.
#'
#' @inheritParams compute_caus_order
#' @param est_g Numeric matrix. The (estimated) adjacency matrix of a
#' DAG or CPDAG.
#' @return  Numeric --- between 0 and 1. The structural intervention
#' distance between the DAG \code{g} and the DAG \code{est_g}.
compute_str_int_distance <- function(g, est_g){

  # check if g is a (non-weighted) adjacency matrix
  if (!all(g %in% c(0, 1))){
    stop("The entries of g must be either 0 or 1.")
  }

  # check if est_g is a (non-weighted) adjacency matrix
  if (!all(est_g %in% c(0, 1))){
    stop("The entries of est_g must be either 0 or 1.")
  }

  p <- NROW(g)
  s <- SID::structIntervDist(g, est_g)
  return(s$sidLowerBound / (p * (p - 1)))
}


#' Convert a causal order into a fully connected DAG
#'
#' Convert the given causal order \code{caus_order} into a fully connected
#' DAG.
#'
#' @inheritParams check_caus_order
#' @return Numeric matrix. The adjacency matrix of the fully connected DAG.
#'
caus_order_to_adjmat <- function(caus_order){

  # check if caus_order has some NA
  if (any(is.na(caus_order))){
    stop("The causal order caus_order cannot contain NA.")
  }

  p <- length(caus_order)
  adj_mat <- upper.tri(x = matrix(0, nrow = p, ncol = p)) * 1
  inv_caus_order <- order(caus_order)
  return(adj_mat[inv_caus_order, inv_caus_order])
}
