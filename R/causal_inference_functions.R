#' Perform causal discovery
#'
#' Perform causal discovery using a particular \code{method}.
#'
#' @param method String. Is one of:
#' \itemize{
#' \item \code{"fast"}, see \code{\link{fast_perm_search}}.
#' \item \code{"full"}, see \code{\link{full_perm_search}}.
#' \item \code{"greedy"}, see \code{\link{greedy_perm_search}}.
#' \item \code{"lingam"}, see \code{\link{lingam_search}}.
#' \item \code{"maxmin"}, see \code{\link{fast_perm_search}}.
#' \item \code{"minimax"}, see \code{\link{minimax_search}}.
#' \item \code{"oracle"}, see \code{\link{oracle_search}}.
#' \item \code{"pc"}, see \code{\link{pc_search}}.
#' \item \code{"random"}, see \code{\link{random_perm_search}}.
#' }
#' @param ... The argument for the given \code{method}. Provide:
#' \itemize{
#' \item \code{delta = <delta_matrix>},
#'  when \code{method} is \code{"fast", "full", "greedy", "maxmin"}.
#' \item \code{gamma = <gamma_matrix>},
#'  when \code{method} is \code{"minimax"}.
#' \item \code{mat = <observations_matrix>},
#'  when \code{method} is \code{"lingam", "pc"}.
#' \item \code{g = <adjacency_matrix>},
#'  when \code{method} is \code{"oracle", "random"}.
#' }
#' @return List. The list is made of:
#' \itemize{
#' \item \code{est_g} --- Numeric matrix (or \code{NA} in case of error).
#' The estimated adjacency matrix of a DAG or CPDAG.
#' \item \code{order} --- Numeric vector (or \code{NA} in case of error).
#' The causal order of the DAG \code{est_g}. If \code{est_g} is a CPDAG,
#' it is first converted to a DAG in order to obtain the causal order.
#' }
causal_discovery <- function(method = c("fast", "full", "greedy",
                                        "lingam", "maxmin", "minimax",
                                        "oracle", "pc", "random")[6],
                             ...){

  # Collect inputs
  argms <- list(...)

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
             out$est_g <- caus_order_to_adjmat(out$order)

           },
           "full" = {
             p <- NROW(argms$delta)
             if (p < 10){
               out$order <- full_perm_search(argms$delta, silent = TRUE)$order
               out$est_g <- caus_order_to_adjmat(out$order)
             }
           },
           "greedy" = {
             out$order <- greedy_perm_search(argms$delta, silent = TRUE)$order
             out$est_g <- caus_order_to_adjmat(out$order)
           },
           "maxmin" = {
             out$order <- fast_perm_search(argms$delta, silent = TRUE,
                                           mode = "maxmin")$order
             out$est_g <- caus_order_to_adjmat(out$order)

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
             out$est_g <- caus_order_to_adjmat(out$order)

           },
           "random" = {
             out$order <- random_perm_search(argms$g)
             out$est_g <- caus_order_to_adjmat(out$order)

           })

  } else if (method %in% c("minimax")){

    if (is.null(argms$gamma)){
      stop(paste("Please, provide the argument",
                 "gamma = <gamma_matrix> for", method, "method."))
    }

    switch(method,
           "minimax" = {
             out$order <- minimax_search(argms$gamma)
             out$est_g <- caus_order_to_adjmat(out$order)

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
