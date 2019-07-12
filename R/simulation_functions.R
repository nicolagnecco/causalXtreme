simulateData <- function(n, p, prob_connect, distr = c("student_t", "log_norm", "norm")[1],
                         tail_index = 1.5, has_confounder = FALSE,
                         is_nonlinear = FALSE, has_uniform_margins = FALSE){
  ## integer integer character numeric_matrix boolean -> list
  ## produces a list given n observations, p variables, a certain distribution,
  ## its degrees of freedom a weighted adjacency matrix and whether to include
  ## some confounders or not.
  ## If adj.mat is not specified a random one will be created.
  ## The produced list is made of:
  ## - data: simulated data with noise following distr
  ## - adj.mat: weighted adjacency matrix of the simulated data
  ## ASSUME:
  ## 1. distr is one of:
  ##      - 'student.t': t-Student
  ##      - 'log.norm  : log-Normal
  ##      - 'norm'     : Normal
  ##      - 'unif'     : Uniform
  ## 2. adj.matr has non-negative weights

  dag <- random_dag(p = p, prob_connect = prob_connect)

  if (has_confounder){
    ll <- add_random_confounders(dag, prob_confound = 2/(3*(p-1)))
    dag_confounders <- ll$dag_confounders
    pos_confounders <- ll$pos_confounders
    p <- p + length(pos_confounders)
  }

  adj_mat <- random_coeff(dag_confounders)

  noise <- simulate_noise(n, p, distr, tail_index)

  # Simulate data
  if(is_nonlinear){
    dataset <- nonlinear_scm(adj_mat, noise)
  } else {
    dataset <- t(solve(diag(p) - t(adj_mat), t(noise)))
  }

  # Marginally transform each variable?
  if(has_uniform_margins){
    dataset <- apply(dataset, 2, unif)
  }

  # Return list
  ll <- list()
  ll$dataset <- dataset
  ll$dag <- dag_confounders
  ll$pos_confounders <- pos_confounders
  return(ll)

}

