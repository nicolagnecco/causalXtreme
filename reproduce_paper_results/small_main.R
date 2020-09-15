# Simulations ----
rm(list = ls())
# devtools::install_github("nicolagnecco/causalXtreme")
library(causalXtreme)
library(tidyverse)
library(doParallel)
library(doRNG)
library(rngtools)
library(tictoc)

# set is_demo = FALSE if you want to reproduce full simulations
# (it can take up to 10 hours)
is_demo <- FALSE

# set is_parallel = TRUE if you want parallel implementation of Simulation 1
# (faster for all algorithms except direct_lingam when n and p are high)
is_parallel <- FALSE




# Simulation 4 --- sensitivity analysis of Rank PC to significance level
source("simulation_4.R")
file.name <- paste("output/simulations_rankpc", ".rds", sep = "")
simulation_4(log_file = "output/sims_4.txt", result_file = file.name,
             is_demo = is_demo)
