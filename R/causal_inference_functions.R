#' Perform causal discovery
#'
#' Perform causal discovery using a particular \code{method}.
#' For more information see 'Details'.
#'
#' This method is a wrapper around the individual causal
#' search functions. For each causal method, it returns a list
#' that can be passed directly to the function \code{\link{causal_metrics}}
#' for evaluation. In particular, the first element of the returned list
#' is a DAG (or CPDAG) and it is used to compute the structural intervention
#' distance (see \code{\link{compute_str_int_distance}}).
#' The second element of the list is a CPDAG and it is used to compute the
#' structural Hamming distance (see \code{\link{compute_str_ham_distance}}).
#'
#' @inheritParams causal_tail_matrix
#' @param method String. Is one of:
#' \itemize{
#' \item \code{"ease"} (the default choice), see \code{\link{ease}}.
#' \item \code{"lingam"}, see \code{\link{lingam_search}}.
#' \item \code{"pc"}, see \code{\link{pc_search}}.
#' \item \code{"pc_rank"}, see \code{\link{pc_rank_search}}.
#' \item \code{"random"}, see \code{\link{random_search}}.
#' }
#' @param set_args List containing the \emph{named} arguments
#' to be passed to the respective method's function. If an argument
#' is missing, then it is set to its default value. If \code{set_args}
#' is not specified, or if it is passed as an empty list,
#' then all arguments of the underlying method are set to their default values.
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
#' @export
causal_discovery <- function(dat, method = c("ease", "lingam", "order_lingam",
                                             "pc", "pc_rank", "random"),
                             set_args = list()){

  # check method
  method <- match.arg(method)

  # set up output list
  out <- list(est_g = NA, est_cpdag = NA)

  # run causal search
  switch(method,
         "ease" = {

           # check arguments
           if (length(set_args) == 0){

             caus_order <- ease(dat)

           } else if (length(set_args) <= 2){

             if (all(names(set_args) %in% c("k", "both_tails"))){

               caus_order <- purrr::pmap(set_args, ease, dat = dat)[[1]]

             } else{

               stop(paste("Arguments for", toupper(method),
                          "must be 'k' and 'both_tails'."))
             }
           } else{

             stop(paste(toupper(method), "accepts at most 2 arguments."))
           }

           # compute DAG/CPDAG and CPDAG
           out$est_g <- caus_order_to_dag(caus_order)
           out$est_cpdag <- dag_to_cpdag(out$est_g)

         },
         "lingam" = {

           # check arguments
           if (length(set_args) == 0){

             dag <- lingam_search(dat)

           } else if (length(set_args) <= 1){

             if (all(names(set_args) %in% c("contrast_fun"))){

               dag <- purrr::pmap(set_args, lingam_search, dat = dat)[[1]]

             } else{

               stop(paste("Arguments for", toupper(method),
                          "must be 'contrast_fun'."))
             }
           } else{

             stop(paste(toupper(method), "accepts at most 1 argument."))
           }

           # compute DAG/CPDAG and CPDAG
           out$est_g <- dag
           out$est_cpdag <- if (all(is.na(out$est_g))) {
             NA
           } else {
             dag_to_cpdag(out$est_g)
           }

         },
         "order_lingam" = {

           # check arguments
           if (length(set_args) == 0){

             caus_order <- order_lingam_search(dat)

           } else if (length(set_args) <= 1){

             if (all(names(set_args) %in% c("contrast_fun"))){

               caus_order <- purrr::pmap(set_args,
                                         order_lingam_search, dat = dat)[[1]]

             } else{

               stop(paste("Arguments for", toupper(method),
                          "must be 'contrast_fun'."))
             }
           } else{

             stop(paste(toupper(method), "accepts at most 1 argument."))
           }

           # compute DAG/CPDAG and CPDAG
           out$est_g <- if(all(is.na(caus_order))){
             NA
           } else {
             caus_order_to_dag(caus_order)
           }

           out$est_cpdag <- if (all(is.na(out$est_g))) {
             NA
           } else {
             dag_to_cpdag(out$est_g)
           }

         },
         "pc" = {

           # check arguments
           if (length(set_args) == 0){

             cpdag <- pc_search(dat)

           } else if (length(set_args) <= 1){

             if (all(names(set_args) %in% c("alpha"))){

               cpdag <- purrr::pmap(set_args, pc_search, dat = dat)[[1]]

             } else{

               stop(paste("Arguments for", toupper(method),
                          "must be 'alpha'."))
             }
           } else{

             stop(paste(toupper(method), "accepts at most 1 argument."))
           }

           # compute DAG/CPDAG and CPDAG
           out$est_g <- cpdag
           out$est_cpdag <- if (all(is.na(out$est_g))) {
             NA
           } else {
             out$est_g
           }
         },
         "pc_rank" = {

           # check arguments
           if (length(set_args) == 0){

             cpdag <- pc_rank_search(dat)

           } else if (length(set_args) <= 1){

             if (all(names(set_args) %in% c("alpha"))){

               cpdag <- purrr::pmap(set_args, pc_rank_search, dat = dat)[[1]]

             } else{

               stop(paste("Arguments for", toupper(method),
                          "must be 'alpha'."))
             }
           } else{

             stop(paste(toupper(method), "accepts at most 1 argument."))
           }

           # compute DAG/CPDAG and CPDAG
           out$est_g <- cpdag
           out$est_cpdag <- if (all(is.na(out$est_g))) {
             NA
           } else {
             out$est_g
           }
         },
         "random" = {

           # check arguments
           if (length(set_args) == 0){

             dag <- random_search(dat)

           } else {

             stop(paste(toupper(method), "accepts no arguments."))
           }

           # compute DAG/CPDAG and CPDAG
           out$est_g <- dag
           out$est_cpdag <- dag_to_cpdag(out$est_g)
         })

  # return list
  return(out)
}


#' Causal evaluation metrics
#'
#' Evaluate the output of the causal inference method called by
#' \code{\link{causal_discovery}}.
#'
#' The evaluation is done with respect to
#' the structural intervention
#' distance (see \code{\link{compute_str_int_distance}}).
#' and the
#' structural Hamming distance (see \code{\link{compute_str_ham_distance}}).
#'
#' @param simulated_data List returned by \code{\link{simulate_data}}.
#' The list is made of:
#' \itemize{
#' \item \code{dataset} --- Numeric matrix. Dataset of simulated data with
#' \code{n} rows and \code{p} columns (note that the hidden variables are not
#' included in this matrix).
#' \item \code{dag} --- Square binary matrix. The generated DAG, including
#' both the observed variables and the confounders,
#' if the argument \code{has_confounder = TRUE} when calling
#' \code{\link{simulate_data}}.
#' \item \code{pos_confounders} --- Integer vector. Represents the position
#' of confounders (rows and columns) in \code{dag}.
#' If the argument \code{has_confounder = FALSE} when calling
#' \code{\link{simulate_data}}, then \code{pos_confounders = integer(0)}.
#' }
#' @param estimated_graphs List returned by \code{\link{causal_discovery}}.
#' The list is made of:
#' \itemize{
#' \item \code{est_g} --- Square binary matrix
#' (or \code{NA} in case of error).
#' The estimated DAG (or CPDAG when the method is \code{pc} or
#' \code{pc_rank}).
#' \item \code{est_cpdag} --- Square binary matrix
#' (or \code{NA} in case of error).
#' The estimated CPDAG.
#' }
#'
#' @return List. The list is made of:
#' \itemize{
#' \item \code{sid} Numeric --- between 0 and 1 --- (or \code{NA} if
#' \code{est_g} is \code{NA}). The structural intervention distance between
#' the true DAG \code{dag} and the estimated DAG (or CPDAG) \code{est_g}.
#' See also \code{\link{compute_str_int_distance}}.
#' \item \code{shd} Numeric --- between 0 and 1 --- (or \code{NA} if
#' \code{est_cpdag} is \code{NA}). The structural Hamming distance
#' between the true CPDAG (\code{dag_to_cpdag(dag)})
#' and the estimated CPDAG \code{est_cpdag}.
#' See also \code{\link{compute_str_ham_distance}}.
#' }
#' @export
causal_metrics <- function(simulated_data, estimated_graphs){

  # collect variables
  true_dag <- simulated_data$dag
  pos_confounders <- simulated_data$pos_confounders
  est_g <- estimated_graphs$est_g
  est_cpdag <- estimated_graphs$est_cpdag

  # set up list
  out <- list(sid = NA, shd = NA)

  # If there are any NA
  if (any(is.na(est_g)) | any(is.na(est_cpdag))){
    out$sid <- NA
    out$shd <- NA

    # return list
    return(out)
  }

  # If there are no confounders
  if (length(pos_confounders) == 0) {
    # SID: compute SID between true DAG and estimated DAG/CPDAG
    out$sid <- compute_str_int_distance(true_dag, est_g)

    # SHD: cast true DAG into CPDAG
    true_cpdag <- dag_to_cpdag(true_dag)
    out$shd <- compute_str_ham_distance(true_cpdag, est_cpdag)
  }

  # If there are confounders
  if (length(pos_confounders) > 0) {
    # SID: extend estimated DAG/CPDAG so that they have the confounders
    extended_est_g <- true_dag
    extended_est_g[-pos_confounders, -pos_confounders] <- est_g
    sid <- compute_str_int_distance(true_dag, extended_est_g)
    p_obs <- NCOL(est_g)
    p <- NCOL(true_dag)
    # adjust for the number of observed variables (not for the total number of
    # variables p = p_obs + p_confounded)
    out$sid <- sid * (p * (p - 1)) / (p_obs * (p_obs - 1))


    # SHD: remove confounders from the true DAG and cast it into CPDAG
    reduced_true_dag <- true_dag[-pos_confounders, -pos_confounders]
    reduced_true_cpdag <- dag_to_cpdag(reduced_true_dag)
    out$shd <- compute_str_ham_distance(reduced_true_cpdag, est_cpdag)
  }

  # return list
  return(out)
}
