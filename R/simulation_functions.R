#' Generate data from a random DAG
#'
#' Generate data from a random DAG.
#' @inheritParams simulate_noise
#' @inheritParams random_dag
#' @param has_confounder Boolean. Are there confounders in the system?
#' @param is_nonlinear Boolean. Is the data generated non linear?
#' @param has_uniform_margins Boolean. Are the variables rescaled uniformly
#' between 0 and 1?
#' @return List. The list is made of:
#' \itemize{
#' \item \code{dataset} --- Numeric matrix. Dataset of simulated data with
#' \code{n} rows and \code{p} columns (note that the hidden variables are not
#' included in this matrix).
#' \item \code{dag} --- Square binary matrix. The generated DAG, including
#' both the observed variables and the confounders,
#' if \code{has_confounder = TRUE}.
#' \item \code{pos_confounders} --- Integer vector. Represents the position
#' of confounders (rows and columns) in \code{dag}.
#' If \code{has_confounder = FALSE}, then \code{pos_confounders = integer(0)}.
#' }
#' @export
simulate_data <- function(n, p, prob_connect,
                         distr = c("student_t", "log_norm", "norm")[1],
                         tail_index = 1.5, has_confounder = FALSE,
                         is_nonlinear = FALSE, has_uniform_margins = FALSE){

  if (p <= 1 | n <= 1){
    stop("n and p must be larger than 1!")
  }

  # Simulate random DAG
  dag <- random_dag(p = p, prob_connect = prob_connect)

  # Add confounders to DAG
  if (has_confounder){
    ll <- add_random_confounders(dag, prob_confound = 2 / (3 * (p - 1)))
    dag <- ll$dag_confounders
    pos_confounders <- ll$pos_confounders
    p <- p + length(pos_confounders)
  } else{
    dag <- dag
    pos_confounders <- integer(0)
  }

  # Create random adjacency matrix
  adj_mat <- random_coeff(dag)

  # Simulate the noise variables
  noise <- simulate_noise(n, p, distr, tail_index)

  # Simulate data
  if (is_nonlinear){
    dataset <- nonlinear_scm(adj_mat, noise)
  } else {
    dataset <- t(solve(diag(p) - t(adj_mat), t(noise)))
  }

  # Marginally transform each variable?
  if (has_uniform_margins){
    dataset <- apply(dataset, 2, uniform_margin)
  }

  # Remove confounders
  if (has_confounder){
    if (length(pos_confounders) > 0){
      dataset <- dataset[, -pos_confounders]
    }
  }

  # Return list
  ll <- list()
  ll$dataset <- dataset
  ll$dag <- dag
  ll$pos_confounders <- pos_confounders
  return(ll)
}
