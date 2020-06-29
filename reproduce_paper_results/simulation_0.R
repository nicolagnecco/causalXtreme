## Project: Causal discovery in heavy-tailed models
## Title: Simulation 0 --- Robustness of hyperparameter k
## Authors: Nicola Gnecco, Nicolai Meinshausen, Jonas Peters, Sebastian Engelke
## NOTE: this code runs in parallel

simulation_0 <- function(log_file,
                         result_file,
                         is_demo = TRUE){


  ## ROBUSTNESS OF K ####
  set.seed(1436) # Gutenberg's press

  # simulation settings
  nexp <- 10
  my_args <- list(
    experiment = 1:nexp,
    n = c(5e2, 1e3, 1e4),
    p = c(4, 7, 10, 15, 20, 30, 50),
    distr = c('student_t'),
    tail_index = c(1.5, 2.5, 3.5),
    has_confounder = c(F),
    is_nonlinear = c(F),
    has_uniform_margins = c(F))

  n <- c("500" = 5e2, "1000" = 1e3, "10000" = 1e4)
  root <- seq(.2, .69, by = 5e-2)
  n_root <- as_tibble(cbind(sapply(n, function(n) floor(n ** root)), root)) %>%
    gather('500', '1000', '10000', key = 'n', value = 'k') %>%
    mutate(n = as.numeric(n))


  my_args <- expand.grid(my_args, stringsAsFactors = F) %>%
    as_tibble() %>%
    left_join(n_root, by = 'n') %>%
    rowwise() %>%
    mutate(prob_connect = min(5/(p - 1), 1/2)) %>%
    filter(has_confounder + is_nonlinear + has_uniform_margins <= 1) %>%
    mutate(id = group_indices())

  niter <- NROW(my_args)

  if (is_demo){
    niter <- 100
  }

  # create clusters
  cores <-  if(Sys.info()["user"] == "elvis"){
    detectCores()
  } else {
    detectCores() - 1
  }
  cl <- makeCluster(cores, outfile = "")
  registerDoParallel(cl)

  # update each worker's environment
  clusterEvalQ(cl, {
    library(causalXtreme)
    library(doParallel)
    library(tidyverse)
    source("simulation_functions.R")
  })


  # Loop through all simulations
  tic()
  sink(file = log_file)
  cat("**** Simulation 0 **** \n")

  ll <- foreach(i = 1:niter, .combine = rbind) %dorng% {
    cat("Experiment number", i, "out of", niter, "\n",
        file = log_file, append = TRUE)

    # Generate data
    current_exp <- my_args[i, ]
    args_simulate <- current_exp %>% select(-experiment, -id, -root, -k)
    X <- do.call(simulate_data, args_simulate)

    # Methods arguments
    args_methods <- tibble(method = c("ease"),
                           set_args = list(
                             list(k = current_exp$k)))

    # Run algorithms
    out <- pmap(args_methods, causal_discovery_wrapper, dat = X$dataset)
    time_elapsed <- map(out, "time") %>% unlist()
    algo_results <- map(out, magrittr::extract, c("est_g", "est_cpdag"))

    # Evaluate algorithms
    algo_evaluations <- algo_results %>%
      map_dfr(causal_metrics, simulated_data = X)

    # Collect results
    res <- algo_evaluations %>%
      mutate(time = time_elapsed) %>%
      mutate(method = c("ease")) %>%
      mutate(id = current_exp$id) %>%
      left_join(current_exp, by = "id") %>%
      group_by_at(vars(-method, -time, -sid, -shd)) %>%
      nest()



  }
  sink()
  toc()
  stopCluster(cl)
  closeAllConnections()

  # save results
  saveRDS(ll, file = result_file)

}
