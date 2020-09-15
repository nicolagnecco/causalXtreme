## Project: Causal discovery in heavy-tailed models
## Title: Simulation 1 --- Algorithms performance
## Authors: Nicola Gnecco, Nicolai Meinshausen, Jonas Peters, Sebastian Engelke

simulation_1 <- function(log_file,
                         result_file,
                         is_demo = TRUE){

  source("simulation_functions.R")

  ## BASIC SIMULATIONS ####
  # Simulation arguments
  settings <- set_simulations(seed = 1321) # Dante's Paradise

  method_args <- settings$method_arguments
  sims_args <- settings$simulation_arguments

  if (is_demo){
    set.seed(42)
    sims_args <- sims_args %>%
      ungroup() %>%
      filter(n %in% c(500, 100),
             p %in% c(4, 10)) %>%
      sample_n(20)
  }
  rm(settings)

  m <- NROW(sims_args)
  k <- NROW(method_args)

  # Loop through all simulations
  inds <- expand.grid(j = 1:k, i = 1:m)
  sink(file = log_file)
  cat("**** Simulation 1 **** \n")
  ll <- map2_dfr(inds$i, inds$j, wrapper_sim, sims_args, method_args)
  sink()
  closeAllConnections()

  # Collect results
  ll <- ll %>%
    ungroup() %>%
    left_join(sims_args, by = "id")

  # save results
  saveRDS(ll, file = result_file)
}
