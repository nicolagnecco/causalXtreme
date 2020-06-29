## Project: Causal discovery in heavy-tailed models
## Title: Simulation 2 --- Computing times
## Authors: Nicola Gnecco, Nicolai Meinshausen, Jonas Peters, Sebastian Engelke


simulation_2 <- function(log_file,
                         result_file,
                         is_demo = TRUE){

  source("simulation_functions.R")

  ## BASIC SIMULATIONS ####
  # Simulation arguments
  set.seed(42)
  exp_ids <- (simulation_settings() %>%
                filter(p %in% c(4, 10, 20, 50)) %>%
                ungroup() %>%
                group_by(n, p) %>%
                sample_n(10))$id
  settings <- set_simulations(seed = 42, experiment_ids = exp_ids,
                              method_nms = c("ease", "direct_lingam",
                                             "ica_lingam", "pc_rank"))
  method_args <- settings$method_arguments

  if (is_demo) {
    method_args <- method_args %>% filter(method == "ease")
  }

  sims_args <- settings$simulation_arguments
  rm(settings)

  m <- NROW(sims_args)
  k <- NROW(method_args)

  # Loop through all simulations
  inds <- expand.grid(j = 1:k, i = 1:m)
  sink(file = log_file)
  cat("**** Simulation 2 **** \n")
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
