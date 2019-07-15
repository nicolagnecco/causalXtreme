#' Perform causal discovery
#'
#' Perform causal discovery using a particular \code{method}.
#' For more information see 'Details'.
#'
#' This method is a wrapper around the individual causal
#' search functions. For each causal method, it returns a list
#' that can be passed directly to the function \code{\link{causal_metrics}}
#' for evaluation. In particular, the first entry of the returned list
#' is a DAG (or CPDAG) and it is used to compute the structural intervention
#' distance (see \code{\link{compute_str_int_distance}}).
#' The second entry is a CPDAG and it is used to compute the
#' structural Hamming distance (see \code{\link{compute_str_ham_distance}}).
#'
#' @inheritParams causal_tail_matrix
#' @param method String. Is one of:
#' \itemize{
#' \item \code{"greedy"}, see \code{\link{greedy_ancestral_search}}.
#' \item \code{"lingam"}, see \code{\link{lingam_search}}.
#' \item \code{"pc"}, see \code{\link{pc_search}}.
#' \item \code{"pc_rank"}, see \code{\link{pc_rank_search}}.
#' \item \code{"random"}, see \code{\link{random_search}}.
#' }
#' @param ... Parameters to be passed to the respective method's function.
#' @return List. The list is made of:
#' \itemize{
#' \item \code{est_g} --- Square binary matrix
#' (or \code{NA} in case of error).
#' The estimated DAG (or CPDAG when the method is \code{pc} or
#' \code{pc_rank}).
#' \item \code{est_cpdag} --- Square binary matrix
#' (or \code{NA} in case of error).
#' The estimated CPDAG.
#' }
causal_discovery <- function(dat, method = c("greedy", "lingam", "pc",
                                        "pc_rank", "random")[1],
                             ...){

  # Collect inputs
  argms <- list(...)
  # !!! continue from here

  # set up list
  out <- list(est_g = NA, order = NA)

  # compute causal order
  if (method %in% c("fast", "full", "greedy", "maxmin")){

    if (is.null(argms$delta)){
      stop(paste("Please, provide the argument",
                 "delta = <delta_matrix> for", method, "method."))
    }

    switch(method,
           "fast" = {
             out$order <- fast_perm_search(argms$delta, silent = TRUE,
                                           mode = "sum")$order
             out$est_g <- caus_order_to_dag(out$order)

           },
           "full" = {
             p <- NROW(argms$delta)
             if (p < 10){
               out$order <- full_perm_search(argms$delta, silent = TRUE)$order
               out$est_g <- caus_order_to_dag(out$order)
             }
           },
           "greedy" = {
             out$order <- greedy_perm_search(argms$delta, silent = TRUE)$order
             out$est_g <- caus_order_to_dag(out$order)
           },
           "maxmin" = {
             out$order <- fast_perm_search(argms$delta, silent = TRUE,
                                           mode = "maxmin")$order
             out$est_g <- caus_order_to_dag(out$order)

           })

  } else if (method %in% c("lingam", "pc")){

    if (is.null(argms$mat)){
      stop(paste("Please, provide the argument",
                 "mat = <observations_matrix> for", method, "method."))
    }

    switch(method,
           "lingam" = {
             out$est_g <- lingam_search(argms$mat) # DAG
             out$order <- if (any(is.na(out$est_g))){
               NA
             } else {
               compute_caus_order(out$est_g)
             }

           },
           "pc" = {
             out$est_g <- pc_search(argms$mat) # CPDAG
             dag <- out$est_g * (t(out$est_g) == 0)
             out$order <- if (any(is.na(out$est_g))){
               NA
             } else {
               compute_caus_order(dag)
             }

           })

  } else if (method %in% c("oracle", "random")){

    if (is.null(argms$g)){
      stop(paste("Please, provide the argument",
                 "g = <adjacency_matrix> for", method, "method."))
    }

    switch(method,
           "oracle" = {
             out$order <- oracle_search(argms$g)
             out$est_g <- caus_order_to_dag(out$order)

           },
           "random" = {
             out$order <- random_perm_search(argms$g)
             out$est_g <- caus_order_to_dag(out$order)

           })

  } else if (method %in% c("minimax")){

    if (is.null(argms$gamma)){
      stop(paste("Please, provide the argument",
                 "gamma = <gamma_matrix> for", method, "method."))
    }

    switch(method,
           "minimax" = {
             out$order <- minimax_search(argms$gamma)
             out$est_g <- caus_order_to_dag(out$order)

           })

  } else {
    stop("Wrong method. Enter one of 'fast', 'full', 'greedy',
              'lingam', 'maxmin', 'minimax', 'oracle', 'random'.")
  }

  # return list
  return(out)
}

#' Causal evaluation metrics
#'
#' Evaluate the output of some causal inference model,
#' see \code{\link{causal_discovery}}, with different metrics.
#'
#' @inheritParams compute_str_int_distance
#' @inheritParams compute_caus_order
#' @inheritParams check_caus_order
#'
#' @references List. The list is made of:
#' \itemize{
#' \item \code{str_int_dist} Numeric --- between 0 and 1 --- (or \code{NA} if
#' \code{est_g} is \code{NA}). The structural intervention distance between
#' the DAG \code{g} and the estimated adjacency matrix \code{est_g}. See also
#' \code{\link{compute_str_int_distance}}.
#' \item \code{ancestral_dist} Numeric --- between 0 and 1 --- (or \code{NA} if
#' \code{order} is \code{NA}). The ancestral distance between the DAG \code{g}
#' and the estimated order \code{order}.
#' See also \code{\link{compute_ancestral_distance}}.
#' }
causal_metrics <- function(g, est_g, caus_order){

  # set up list
  out <- list(str_int_dist = NA, ancestral_dist = NA)

  # compute structural intervention distance (if caus_order is not NA)
  if (!any(is.na(est_g))){
    out$str_int_dist <- compute_str_int_distance(g, est_g)
  } else {
    out$str_int_dist <- NA
  }

  # compute ancestral distance (if est_g is not NA)
  if (!any(is.na(caus_order))){
    out$ancestral_dist <- compute_ancestral_distance(g, caus_order)
  }else {
    out$ancestral_dist <- NA
  }

  # return list
  return(out)
}
