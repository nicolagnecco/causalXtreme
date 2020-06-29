## Project: Causal discovery in heavy-tailed models
## Title: Simulation 3 --- ICA-LiNGAM and Pairwise LiNGAM performance on
##        Setting 3, where we discard 50% of the data in the bulk.
## Authors: Nicola Gnecco, Nicolai Meinshausen, Jonas Peters, Sebastian Engelke

simulation_3 <- function(log_file,
                         result_file,
                         is_demo = TRUE){

  source("simulation_functions.R")

  ## BASIC SIMULATIONS ####
  # Simulation arguments
  exp_ids <- (simulation_settings() %>%
                filter(is_nonlinear == TRUE,
                       p %in% c(10, 20, 30, 50),
                       n == 1e4))$id
  settings <- set_simulations(seed = 1321, # Dante's Paradise
                              experiment_ids = exp_ids,
                              method_nms = c("ica_lingam", "direct_lingam"))

  if (is_demo) {
    settings <- set_simulations(seed = 1321,
                                method_nms = c("ica_lingam", "direct_lingam"))
    # Dante's Paradise
  }

  method_args <- settings$method_arguments
  sims_args <- settings$simulation_arguments

  if (is_demo){
    set.seed(42)
    sims_args <- sims_args %>%
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
  cat("**** Simulation 3 **** \n")
  ll <- map2_dfr(inds$i, inds$j, wrapper_sim, sims_args, method_args, TRUE)
  sink()
  closeAllConnections()

  # Collect results
  ll <- ll %>%
    ungroup() %>%
    left_join(sims_args, by = "id")

  # save results
  saveRDS(ll, file = result_file)
}
