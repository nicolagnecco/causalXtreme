## Project: Causal discovery in heavy-tailed models
## Title: Simulation 1 --- Algorithms performance
## Authors: Nicola Gnecco, Nicolai Meinshausen, Jonas Peters, Sebastian Engelke
## NOTE: this code runs in parallel


simulation_1_par <- function(log_file,
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

  m <- NROW(sims_args); m <- 5
  k <- NROW(method_args)

  # create clusters
  cores <-  if(Sys.info()["user"] == "elvis"){
    detectCores()
  } else {
    detectCores() - 1
  }
  cl <- makeCluster(cores, type = "PSOCK", outfile = "")
  registerDoParallel(cl)

  # update each worker's environment
  clusterEvalQ(cl, {
    library(causalXtreme)
    library(doParallel)
    library(tidyverse)
    source("simulation_functions.R")
  })

  clusterExport(cl, c("sims_args", "method_args"), envir = environment())


  # Loop through all simulations
  tic()
  sink(file = log_file)
  cat("**** Simulation 1 **** \n")

  ll <- foreach(i = 1:m, .combine = bind_rows)  %:%
    foreach(j = 1:k, .combine = bind_rows) %dopar%{

      cat("Simulation", i, "out of", m, "--- Inner iteration", j,
          "out of", k, "\n", file = log_file, append = TRUE)
      wrapper_sim(i, j, sims_args, method_args)

    }
  sink()
  toc()
  stopCluster(cl)
  closeAllConnections()

  # Collect results
  ll <- ll %>%
    ungroup() %>%
    left_join(sims_args, by = "id")

  # save results
  saveRDS(ll, file = result_file)
}
