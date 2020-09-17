## Project: Causal discovery in heavy-tailed models
## Title: Simulation 4 --- Sensitivity analysis of Rank PC to
##        significance level.
## Authors: Nicola Gnecco, Nicolai Meinshausen, Jonas Peters, Sebastian Engelke

simulation_4 <- function(log_file,
                         result_file,
                         is_demo = TRUE){

  source("simulation_functions.R")

  ## BASIC SIMULATIONS ####
  # Simulation arguments
  exp_ids <- (simulation_settings() %>%
                filter(has_confounder == FALSE,
                       is_nonlinear == FALSE,
                       has_uniform_margins == FALSE,
                       p %in% c(4, 10, 30, 50),
                       n == 1e3,
                       tail_index == 3.5))$id
  settings <- set_simulations(seed = 1321, # Dante's Paradise
                              experiment_ids = exp_ids,
                              method_nms = c("pc_rank", "random"))

  if (is_demo) {
    settings <- set_simulations(seed = 1321,
                                method_nms = c("pc_rank", "random"))
    # Dante's Paradise
  }

  method_args <- settings$method_arguments %>%
    mutate(set_args = list(list(list(alpha = 5e-4),
                                list(alpha = 5e-3),
                                list(alpha = 5e-2),
                                list(alpha = 5e-1)),
                                list(list()))) %>%
    unnest(cols = set_args)

  sims_args <- settings$simulation_arguments


  if (is_demo){
    set.seed(42)
    sims_args <- sims_args %>%
      ungroup() %>%
      filter(p %in% c(4, 10)) %>%
      sample_n(20)
  }

  rm(settings)

  m <- NROW(sims_args)
  k <- NROW(method_args)

  # Loop through all simulations
  inds <- expand.grid(j = 1:k, i = 1:m)
  sink(file = log_file)
  cat("**** Simulation 4 **** \n")
  ll <- map2_dfr(inds$i, inds$j, wrapper_sim, sims_args, method_args,
                 trimdata = FALSE, meth_args = TRUE)
  sink()
  closeAllConnections()

  # Collect results
  ll <- ll %>%
    ungroup() %>%
    left_join(sims_args, by = "id")

  # save results
  saveRDS(ll, file = result_file)
}
