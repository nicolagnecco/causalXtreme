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
is_demo <- TRUE

# set is_parallel = TRUE if you want parallel implementation of Simulation 1
# (faster for all algorithms except direct_lingam when n and p are high)
is_parallel <- FALSE


# Simulation 0 --- EASE robustness wrt k
source("simulation_0.R")
file.name <- paste("output/k_robustness", ".rds", sep = "")

simulation_0(log_file = "output/sims_0.txt", result_file = file.name,
             is_demo = TRUE)


# Simulation 1 --- main simulations
file.name <- paste("output/simulations", ".rds", sep = "")

if (is_parallel){
  source("simulation_1_parallel.R")
  simulation_1_par(log_file = "output/sims_1.txt", result_file = file.name,
                   is_demo = is_demo)
} else{
  source("simulation_1.R")
  simulation_1(log_file = "output/sims_1.txt", result_file = file.name,
               is_demo = is_demo)
}


# Simulation 2 --- computing times
source("simulation_2.R")
file.name <- paste("output/time", ".rds", sep = "")
simulation_2(log_file = "output/sims_2.txt", result_file = file.name,
             is_demo = is_demo)


# Simulation 3 --- remove data from bulk from ICA-LiNGAM and Pairwise LiNGAM
source("simulation_3.R")
file.name <- paste("output/simulations_lingam",".rds", sep = "")
simulation_3(log_file = "output/sims_3.txt", result_file = file.name,
             is_demo = is_demo)

# Simulation 4 --- sensitivity analysis of Rank PC to significance level
source("simulation_4.R")
file.name <- paste("output/simulations_rankpc", ".rds", sep = "")
simulation_4(log_file = "output/sims_4.txt", result_file = file.name,
             is_demo = is_demo)


# Plot results
library(latex2exp)
library(tools)
library(kableExtra)
source("produce_charts.R")
produce_charts(sim0_file = "original_output/k_robustness.rds",
               sim1_file = "original_output/simulations.rds",
               sim2_file = "original_output/time.rds",
               sim3_file = "original_output/simulations_lingam.rds",
               sim4_file = "original_output/simulations_rankpc.rds")


# Financial application ----
source("financial_data.R")


# River data ----
source("river_data.R")
source("river_map.R")
