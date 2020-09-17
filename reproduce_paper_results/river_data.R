## Project: Causal discovery in heavy-tailed data
## Descrip: Do analysis on the upper Danube basin discharges
## Authors: Nicola Gnecco, Nicolai Meinshausen, Jonas Peters, Sebastian Engelke

## Load libraries ####
rm(list = ls())
library(causalXtreme)
library(doParallel)
library(doRNG)
library(Hmisc)
library(latex2exp)
require(reshape2)
library(tidyquant)
library(timetk)
library(tsibble)
library(ismev)
library(pracma)

theme_set(theme_bw() +
            theme(plot.background = element_blank(),
                  legend.background = element_blank()))

## Define constants ####
OUTPUT_FILE <- "output/river_results.txt"
tolPalette <- c(tolBlack = "#000000",
                tolBlue = "#4477AA",
                tolRed = "#EE6677",
                tolGreen = "#228833",
                tolYellow = "#CCBB44",
                tolCyan = "#66CCEE",
                tolPurple = "#AA3377",
                tolGrey = "#BBBBBB") %>%
  unname()



## Define functions ####
dfr2tibble <- function(dfr){
  ## dataframe -> tibble
  ## convert dataframe to tibble and rename column "Date" to date

  dfr %>%
    mutate(Date = as_date(Date)) %>%
    rename(date = Date) %>%
    as_tibble()
}

PP.lik.linear <- function(Params, Data, u, NoYears, CovarMu, CovarSc, CovarXi){
  ## Function from the paper "EXTREMES ON RIVER NETWORKS"
  ## By PEIMAN ASADI, ANTHONY C. DAVISON and SEBASTIAN ENGELKE

  NoSt <- length(Data)
  NoParMu <- ncol(CovarMu)
  NoParSc <- ncol(CovarSc)
  NoParXi <- ncol(CovarXi)
  NoPar <- NoParMu + NoParSc + NoParXi
  Out <- 0

  for (i in 1:NoSt) {
    mu <- sum(Params[1:NoParMu]*(CovarMu[i,]))
    sc <- sum(Params[(NoParMu+1):(NoParMu+NoParSc)]*CovarSc[i,])
    xi <- sum(Params[(NoParMu+NoParSc+1):NoPar]*CovarXi[i,])

    y1 <- 1 + xi * ((u[i] - mu)/sc)
    y2 <- 1 + xi * ((Data[[i]]- mu)/sc)
    InterceptSc <- Params[NoParMu+1]
    InterceptMu <- Params[1]

    if ((sc <= 0) | (min(y1) <= 0) | (min(y2) <= 0) |
        InterceptSc <= 0 | InterceptMu < 0 ) {
      l <- 10^6
    } else {
      l <- NoYears[i] * (y1 ^ (-1 / xi)) +
        length(Data[[i]]) * log(sc) + (1 / xi + 1) * sum(log(y2))
    }
    Out <- Out + l
  }
  return(Out)
}

PPFit <- function(Data, u, NoYears, Init, CovarMu, CovarSc, CovarXi,
                  method = method, control = control) {
  ## Function from the paper "EXTREMES ON RIVER NETWORKS"
  ## By PEIMAN ASADI, ANTHONY C. DAVISON and SEBASTIAN ENGELKE

  x <- optim(Init, PP.lik.linear,
             Data = Data, u = u, NoYears = NoYears, CovarMu = CovarMu,
             CovarSc = CovarSc, CovarXi = CovarXi,
             method = method, control = control, hessian = TRUE)
  output <- x
  return(output)

}

split_stations <- function(dat, year = NULL, stations = NULL){
  ## tibble integer character_vector -> list
  ## produce a named list with discharges for the stations in given year
  if (!is.null(year) & !is.null(stations)){
    dat_tmp <- dat %>%
      filter(year(date) %in% year, station %in% stations)

    if (nrow(dat_tmp) == 0){
      return(list())
    }
  } else {
    dat_tmp <- dat
  }

  unique_names <- unique(dat_tmp$station)

  dat_tmp %>%
    select(station, disc) %>%
  group_split() %>%
    set_names(unique_names) %>%
    lapply(function(tbl){tbl$disc})

}

## Import dataset ####
load("data/river_data/StsTSs.RData")
load("data/river_data/StsInfo.RData")


# Select the 31 stations that are in the AOAS paper plus station 19 (renamed 32)
StsChos <- c(c(1:47)[-c(16,30,31,34,42,43,44,45,46,47,3,1,2,29,18,19)], 19)

NoSt <- length(StsChos)

StsTSsChos <- StsTSs[StsChos]
StsInfoChos <- StsInfo[StsChos,] %>%
  mutate(id_old = StsChos,
         id = 1:NoSt)

# clean data and find common dates
river_dat <- map(.x = StsTSsChos, .f = dfr2tibble) %>%
  reduce(.f = inner_join, by = "date") %>%
  filter(month(date) %in% c(6, 7, 8),
         year(date) < 2010)

# rename columns (according to AOAS paper)
colnames(river_dat)[2:(NoSt + 1)] <- sprintf("station_%02d", 1:NoSt)

# select stations
station_names <- c(11, 9, 21, 7, 19, 14, 26, 23, 28, 1, 13, 32)
river_dat <- river_dat[, c(1, station_names + 1)]
NoSt <- length(station_names)

station_info <- StsInfoChos %>%
  filter(id %in% station_names) %>%
  rename(name = RivNames,
         lat = Lat,
         lon = Long,
         ave_vol = AveVol)

# clean workspace
rm(StsInfo, StsInfoChos, StsTSs, StsTSsChos, StsChos)

## Check model assumptions ####
method <- "BFGS"
control <- list(maxit = 5000, reltol=10^(-30), abstol=0.0001, trace=0)

# Take maxima
SummerMaxima <- river_dat %>%
  mutate(year = year(date)) %>%
  gather(key = "station", value = "disc", -date, -year) %>%
  group_by(year, station) %>%
  summarise(maxDisc = max(disc)) %>%
  spread(key = station, value = maxDisc) %>%
  select(year, sprintf("station_%02d", station_names))

# GEV analysis
GevPars <- tibble(station = character(), location = double(),
                  scale = double(), shape = double())

for (i in 1:NoSt){

  vec <- SummerMaxima[[i + 1]]
  station_nm <- colnames(SummerMaxima)[i + 1]
  tmpres <- gev.fit(vec, show = FALSE)
  GevPars <- bind_rows(GevPars,
                       tibble(station = station_nm,
                              location = tmpres$mle[1],
                              scale = tmpres$mle[2],
                              shape = tmpres$mle[3]))
}


# Compute threshold data
q <- .9
threshold_dat <- river_dat %>%
  gather(key = "station", value = "disc", -date) %>%
  mutate(station = factor(station,
                          levels = sprintf("station_%02d", station_names))) %>%
  group_by(station) %>%
  mutate(thres = quantile(disc, q)) %>%
  filter(disc > thres) %>%
  mutate(noyears = length(unique(year(date))))


AllEvents <- split_stations(threshold_dat)
Threshold <- threshold_dat %>% distinct(station, thres) %>% deframe()
NoOfYears <- threshold_dat %>% distinct(station, noyears) %>% deframe()


# Fitting a GEVD to each station (by maximizing the joint Poisson
# process likelihood; cf. formula (21) in the paper of
# Asadi, Davison, Engelke (2015)
ParsPP <- tibble(station = character(), location = double(),
                 scale = double(), shape = double())

for (i in 1:NoSt) {
  Init <- GevPars %>% select(-station) %>% slice(i) %>% as.double()
  station_nm <- GevPars %>% select(station) %>% slice(i) %>% as.character()

  tmpres <- PPFit(Data=AllEvents[i],
                  u=Threshold[i],
                  NoYears=NoOfYears[i],
                  Init=Init,
                  CovarMu=matrix(1,1,1),CovarSc=matrix(1,1,1),
                  CovarXi=matrix(1,1,1),method=method,control=control)
  ParsPP <- bind_rows(ParsPP,
                      tibble(station = station_nm,
                             location = tmpres$par[1],
                             scale = tmpres$par[2],
                             shape = tmpres$par[3]))
}



# Regionalized model from paper of Asadi, Davison, Engelke (2015)
Grp1 <- sprintf("station_%02d", c(11, 19, 21))
Grp2 <- sprintf("station_%02d", c(13, 28, 32))
Grp3 <- sprintf("station_%02d", c(1, 7, 9, 14))
Grp4 <- sprintf("station_%02d", c(23, 26))
Grp <- list(Grp1, Grp2, Grp3, Grp4)

# Estimate parameters: shape, location, scale
ParsCovModel <- tibble(station = character(), region = double(),
                       location = double(), scale = double(),
                       shape = double(), shape_se = double())

hessian_list <- list()

for (i in 1:length(Grp)) {
  GrSts <- Grp[[i]]

  n_st <- length(GrSts)

  CovarMuReg <- diag(n_st)
  CovarScReg <- diag(n_st)
  CovarXi <- matrix(1, nrow = n_st)

  station_nms <- ParsPP %>%
    filter(station %in% GrSts) %>%
    select(station) %>%
    deframe()

  InitMu <- ParsPP %>% filter(station %in% GrSts) %>% select(location) %>%
    deframe()
  InitSc <- ParsPP %>% filter(station %in% GrSts) %>% select(scale) %>%
    deframe()
  InitXi <- ParsPP %>% filter(station %in% GrSts) %>% select(shape) %>%
    summarise(mean(shape)) %>% deframe()
  Init <- c(InitMu, InitSc, InitXi)

  tmpres <- PPFit(Data=AllEvents[names(AllEvents) %in% GrSts],
                             u=Threshold[names(Threshold) %in% GrSts],
                             NoYears=NoOfYears[names(NoOfYears) %in% GrSts],
                             Init=Init,
                             CovarMu=CovarMuReg,
                             CovarSc=CovarScReg,
                             CovarXi=CovarXi,
                             method=method, control=control)


  MuCov <- tmpres$par[1:n_st]
  ScCov <- tmpres$par[1:n_st + n_st]
  XiCov <- tmpres$par[2 * n_st + 1]

  nms <- c(paste("loc_", station_nms, sep = ""),
           paste("scale_", station_nms, sep = ""),
           "shape")
  H <- tmpres$hessian
  dimnames(H) <- list(nms, nms)
  hessian_list[[i]] <- H


  sderror <- sqrt(solve(H)[2 * n_st + 1, 2 * n_st + 1])

  ParsCovModel <- bind_rows(ParsCovModel,
                            tibble(station = station_nms,
                                   region = i, location = MuCov, scale = ScCov,
                                   shape = XiCov, shape_se = sderror))
}



# Compute adjusted se taking into accout time dependence
# (see Fawcett and Walshaw, 2007 and 2012)
se_adj <- tibble(region = double(), shape_se_adj = double())

for (i in 1:length(Grp)) {
  GrSts <- Grp[[i]]

  my_grad <- NULL
  yrs <- unique(year(threshold_dat$date))

  for (j in seq_along(yrs)){
    y <- unique(year(threshold_dat$date))[j]
    AllEvents_year <- split_stations(threshold_dat, year = y, stations = GrSts)

    GrSts_year <- names(AllEvents_year)

    if(length(GrSts) != length(GrSts_year)){
      next()
    }

    cat(GrSts_year, "--- year =", y, "\n")


    n_st <- length(GrSts_year)

    # Initialize parameters
    CovarMuReg <- diag(n_st)
    CovarScReg <- diag(n_st)
    CovarXi <- matrix(1, nrow = n_st)

    station_nms <- ParsPP %>%
      filter(station %in% GrSts_year) %>%
      select(station) %>%
      deframe()

    InitMu <- ParsCovModel %>%
      filter(station %in% GrSts_year) %>%
      select(location) %>%
      deframe()
    InitSc <- ParsCovModel %>%
      filter(station %in% GrSts_year) %>%
      select(scale) %>%
      deframe()
    InitXi <- ParsCovModel %>%
      filter(station %in% GrSts_year) %>%
      select(shape) %>%
      summarise(mean(shape)) %>%
      deframe()

    Init <- c(InitMu, InitSc, InitXi)

    # Compute gradient at optimal paramter value
    tmpres <- pracma::grad(f = PP.lik.linear,
                           x0 = Init, Data=AllEvents_year[GrSts_year],
                           u=Threshold[GrSts_year],
                           NoYears=NoOfYears[GrSts_year],
                           CovarMu=CovarMuReg,CovarSc=CovarScReg,
                           CovarXi=CovarXi)

    my_grad <- rbind(my_grad, tmpres)

  }

  H <- hessian_list[[i]]
  V <- max(NoOfYears[GrSts_year]) * cov(my_grad)

  d <- ncol(my_grad)
  cov_theta <- solve(H) %*% V %*% solve(H)

  se_adj <- bind_rows(se_adj,
                      tibble(region = i,
                             shape_se_adj = sqrt(cov_theta[d, d])))
}

# Paper results
sink(OUTPUT_FILE)
cat("--------------------------------------------------------------------------",
    "\n")
cat("MLE estimates of shape parameter xi and CI (lower bound = lb, upper bound = up)",
    "\n")
cat("--------------------------------------------------------------------------",
    "\n")
ParsCovModel  %>% left_join(se_adj) %>%
  mutate(lb = round(shape - 2 * shape_se_adj, 3),
         ub = round(shape + 2 * shape_se_adj, 3)) %>%
  group_by(region) %>%
  summarise(shape = mean(shape),
            lb = mean(lb),
            ub = mean(ub)) %>%
  print()
cat("\n\n\n")
sink()


## Spatial structure ####
# true DAG
true_dag <- rbind(c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                  c(0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0),
                  c(0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0),
                  c(rep(0, 9), 1, 0, 0),
                  c(0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0),
                  c(rep(0, 9), 1, 0, 0),
                  c(rep(0, 9), 1, 0, 0),
                  c(rep(0, 9), 1, 0, 0),
                  c(rep(0, 9), 0, 1, 0),
                  rep(0, 12),
                  c(rep(0, 9), 1, 0, 0),
                  c(rep(0, 9), 0, 1, 0))
colnames(true_dag) <- rownames(true_dag) <- station_names

g <- igraph::graph_from_adjacency_matrix(true_dag)
igraph::V(g)$color <- "white"
# igraph::tkplot(g)

# data
mat <- river_dat %>%  select(-date) %>% as.matrix()

# ground truth
ll <- list(
  dataset = mat,
  dag = true_dag,
  pos_confounders = integer(0)
)

# run ease
sink(OUTPUT_FILE, append = TRUE)
cat("--------------------------------------------------------------------------",
    "\n")
cat("1. SPATIAL ANALYSIS", "\n")
cat("--------------------------------------------------------------------------",
    "\n")
cat("--------------------------------------------------------------------------",
    "\n")
cat("EASE \n")
cat("--------------------------------------------------------------------------",
    "\n")
ease <- causal_discovery(mat, "ease", set_args = list(both_tails = F))
cat("SID:", causal_metrics(ll, ease)$sid, "\n")
caus_ord <- ease(dat = mat, both_tails = F)
cat("Causal order:", colnames(mat)[caus_ord], "\n")
cat("\n\n")

# run ica_lingam
cat("--------------------------------------------------------------------------",
    "\n")
cat("ICA-LiNGAM \n")
cat("--------------------------------------------------------------------------",
    "\n")
ica_lingam <- causal_discovery(mat, "ica_lingam")
cat("SID:", causal_metrics(ll, ica_lingam)$sid, "\n")
cat("\n\n")

# run direct-lingam
cat("--------------------------------------------------------------------------",
    "\n")
cat("Pairwise LiNGAM \n")
cat("--------------------------------------------------------------------------",
    "\n")
direct_lingam <- causal_discovery(mat, "direct_lingam")
cat("SID:", causal_metrics(ll, direct_lingam)$sid, "\n")
cat("\n\n")


# run pc rank
cat("--------------------------------------------------------------------------",
    "\n")
cat("Rank PC \n")
cat("--------------------------------------------------------------------------",
    "\n")
pc_rank <- causal_discovery(mat, "pc_rank")
cat("SID:", causal_metrics(ll, pc_rank)$sid, "\n")
cat("\n\n\n")
sink()



## Time structure ####
k <- 6
river_dat_lagged <- river_dat %>%
  gather(key = "station", value = "disc", -date) %>%
  mutate(year = year(date)) %>%
  group_by(station, year) %>%
  tq_transmute(
    select = disc,
    mutate_fun = lag.xts,
    k = 0:k
  ) %>%
  ungroup()


dat_lagged <- river_dat_lagged %>%
  drop_na() %>%
  group_by(station, year) %>%
  mutate(counter = 1:n()) %>%
  filter(counter %in% seq(from = 1, to = n(), by = 1)) %>%
  ungroup() %>%
  select(-year, -counter)


true_time_str <- c(paste("disc.", k:1, sep = ""), "disc")
results <- matrix(nrow = length(station_names), ncol = 3)
colnames(results) <- c("ease", "ica_lingam", "direct_lingam")

for (j in 1:length(station_names)){
  current_station <- sprintf("station_%02d", station_names[j])
  mat_single_station <- dat_lagged %>%
    filter(station == current_station) %>%
    select(-station, -date) %>%
    as.matrix()

  # ease
  caus_ord <- ease(dat = mat_single_station, both_tails = F)
  cat("EASE causal order:", caus_ord, "\n")
  results[j, 1] <- all(colnames(mat_single_station)[caus_ord] == true_time_str)

  # ica lingam
  caus_ord <- ica_lingam_search(dat = mat_single_station)
  results[j, 2] <- all(colnames(mat_single_station)[caus_ord] == true_time_str)

  # direct lingam
  caus_ord <- direct_lingam_search(dat = mat_single_station)
  results[j, 3] <- all(colnames(mat_single_station)[caus_ord] == true_time_str)
}

rownames(results) <- colnames(mat)
colnames(results) <- c("EASE", "ICA-LiNGAM", "Pairwise LiNGAM")
sink(OUTPUT_FILE, append = TRUE)
cat("--------------------------------------------------------------------------",
    "\n")
cat("2. TIME ANALYSIS", "\n")
cat("--------------------------------------------------------------------------",
    "\n")
print(results)
cat("\n\n\n")
sink()


## Robustness wrt to k ####
set.seed(1436) # Gutenberg's press

# data
mat <- river_dat %>%
  select(-date) %>% as.matrix()

# varying k
len_root <- 51
n_subsample <- 50
n_iter <- len_root * n_subsample
root_vec <- seq(.2, .7, length.out = len_root)

# create clusters
cores <-  if(Sys.info()["user"] == "elvis"){detectCores()}else{detectCores() - 1}
cl <- makeCluster(cores, outfile = "")
registerDoParallel(cl)

# update each worker's environment
clusterEvalQ(cl, {
  library(causalXtreme)
  library(doParallel)
  library(tidyverse)
})


# Loop through all simulations
out_list <- foreach(i = 1:n_iter, .combine = rbind) %dorng% {
  cat("\r Simulation", i, "out of", n_iter)

  # Compute the fractional exponent
  idx <-  (i - 1) %% (len_root) + 1
  root <- root_vec[idx]

  # Bootstrap data
  idx_rows <- sample(x = 1:NROW(mat), replace = TRUE)
  current_mat <- mat[idx_rows, ]

  # Compute the number of observations in the tail
  k <- floor(NROW(current_mat) ** root)

  # Run EASE on bootstraped sample
  ease <- causal_discovery(current_mat, "ease",
                           set_args = list(both_tails = F, k  = k))
  sid <- causal_metrics(ll, ease)$sid

  # Return results
  res <- c(root = root,
           k = k,
           sid = sid)
}
stopCluster(cl)

# Plot results
dat <- as_tibble(out_list)

alpha <- .1

res <- dat %>%
  group_by(root) %>%
  summarise(mean_sid = mean(sid),
            se = sd(sid),
            lb = quantile(sid, probs = alpha / 2),
            ub = quantile(sid, probs = 1 - alpha / 2)) %>%
  ungroup()


g <- ggplot(res, aes(x = root, y = mean_sid)) +
  geom_ribbon(aes(ymin = lb, ymax = ub),
              alpha = 0.25) +
  geom_line(alpha = .9, size = 1,  color = tolPalette[2]) +
  scale_color_manual(values = unname(tolPalette[5])) +
  xlab(TeX("Fractional exponent of $k_n$")) +
  ylab("Structural Intervention Distance") +
  labs(color = TeX("Tail index $\\alpha$")) +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14)
  ) + ylim(0, 0.3); g

ggsave("output/river_k_robustness.pdf",
       g, width = 7.5, height = 5, units = c("in"))



## Granger causality ####
# Define matrices Y and X
Y <- river_dat %>%
  filter(!(month(date) == 6 & day(date) == 1)) %>%
  select(-date) %>%
  as.matrix(dimnames = colnames(river_dat))

X <- river_dat %>%
  filter(!(month(date) == 8 & day(date) == 31)) %>%
  select(-date) %>%
  as.matrix(dimnames = colnames(river_dat))

# compute robust pvalues
robust_pval <- function(fit){
  v <-  sandwich::sandwich(fit, meat. = sandwich::meatHAC)
  se <- sqrt(diag(v))
  df <- df.residual(fit)
  pval <- 2 * pt(abs(fit$coefficients/se), df = df, lower.tail = FALSE)
  (pval[-1])
}


# Run multivariate regressions
pval_mat <- matrix(nrow = NoSt, ncol = NoSt, dimnames = list(colnames(Y),
                                                             colnames(X)))
for (i in 1:NoSt){
  cat("Iteration", i, "out of", NoSt, "\n")
  y <- Y[, i, drop = FALSE]
  fit <- lm(y ~ X)
  pval <- robust_pval(fit)
  pval_mat[, i] <- pval
}

# take significant pairs (p-value < 0.05 / (# of tests))
significant_pairs <- pval_mat < 5e-2/(NoSt * (NoSt - 1))

# create a weighted adjacency matrix with significant pairs
w_graph <- matrix(Inf, nrow = NoSt, ncol = NoSt)
w_graph[significant_pairs] <- log(pval_mat[significant_pairs]) + 150
diag(w_graph) <- Inf
colnames(w_graph) <- rownames(w_graph) <- colnames(river_dat)[-1]



# build minimum spanning polytree
construct_polytree <- function(adj_mat){
  ## numeric_matrix -> numeric_matrix
  ## compute a polytree (or forest) adding the edges (i, j) in increasing order
  ## of weights, if nodes i and j are not connected.
  ## Stop if the polytree is weakly connected, or if all
  ## the edges have been considered
  p <- nrow(adj_mat)
  adj_mat0 <- adj_mat
  polytree <- matrix(0, nrow = p, ncol = p)

  for (i in seq_len(p)){
    for (j in seq_len(p)){

      # stop if the polytree is weakly connected
      g_tmp <- igraph::graph_from_adjacency_matrix(polytree, mode = "undirected")

      if (igraph::is.connected(g_tmp)){
        return(polytree)
      }

      # stop if all the edges have been considered
      if (min(adj_mat0) == Inf){
        return(polytree)
      }

      # take the edge (i, j) with minimum weight and update adj_mat0
      min_pair <- which(adj_mat0 == min(adj_mat0), arr.ind = T)
      adj_mat0[min_pair] <- Inf

      # if edge (i, j) are already connected skip, else add (i, j)
      ij_connected <- igraph::all_simple_paths(g_tmp,
                                               from = min_pair[1],
                                               to = min_pair[2])
      if (length(ij_connected) == 0){
        polytree[min_pair] <- 1
      } else {
        next()
      }

    }
  }
  return(polytree)
}
dag_est <- construct_polytree(w_graph)
colnames(dag_est) <- rownames(dag_est) <- colnames(river_dat)[-1]

# compute structural intervention distance
sink(OUTPUT_FILE, append = TRUE)
cat("--------------------------------------------------------------------------",
    "\n")
cat("3. GRANGER CAUSALITY ANALYSIS", "\n")
cat("--------------------------------------------------------------------------",
    "\n")
cat("SID polytree:",
    causalXtreme:::compute_str_int_distance(true_dag, dag_est), "\n")
sink()

# build minimum spanning dag
construct_dag <- function(adj_mat){
  ## numeric_matrix -> numeric_matrix
  ## compute a dag adding the edges (i, j) in increasing order of weights
  ## if edge (i, j) does not introduce cycle. Stop once all the edges have
  ## been considered

  # helper function
  find_cycles <- function(g){
    Cycles = NULL
    for(v1 in igraph::V(g)) {
      for(v2 in igraph::neighbors(g, v1, mode="out")) {
        Cycles = c(Cycles,
                   lapply(igraph::all_simple_paths(g, v2,v1, mode="out"),
                          function(p) c(v1,p)))
      }
    }
    return(Cycles)
  }

  # function's body
  p <- nrow(adj_mat)
  adj_mat0 <- adj_mat
  dag <- matrix(0, nrow = p, ncol = p)

  for (i in seq_len(p)){
    for (j in seq_len(p)){
      # stop if all the edges have been considered
      if (min(adj_mat0) == Inf){
        return(dag)
      }

      # take the edge (i, j) with minimum weight and update adj_mat0
      min_pair <- which(adj_mat0 == min(adj_mat0), arr.ind = T)
      adj_mat0[min_pair] <- Inf

      # if edge (i, j) introduce cycle skip, else add (i, j)
      dag_tmp <- dag
      dag_tmp[min_pair] <- 1
      g_tmp <- igraph::graph_from_adjacency_matrix(dag_tmp)
      # plot(g_tmp)
      cycles <- find_cycles(g_tmp)
      if (length(cycles) == 0){
        dag[min_pair] <- 1
      } else {
        next()
      }

    }
  }
}
dag_est <- construct_dag(w_graph)
colnames(dag_est) <- rownames(dag_est) <- colnames(river_dat)[-1]

# compute structural intervention distance
sink(OUTPUT_FILE, append = TRUE)
cat("SID DAG:",
    causalXtreme:::compute_str_int_distance(true_dag, dag_est), "\n")
sink()
