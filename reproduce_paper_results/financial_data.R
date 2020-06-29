## Project: Causal discovery in heavy-tailed data
## Descrip: Do analysis on financial data --- EURCHF exchange rate and major
##          Swiss stocks, i.e., Nestle, Novartis, Roche
## Authors: Nicola Gnecco, Nicolai Meinshausen, Jonas Peters, Sebastian Engelke

## Load libraries ####
rm(list = ls())
library(causalXtreme)
library(cowplot)
library(doParallel)
library(doRNG)
library(evd)
library(evir)
library(Hmisc)
library(latex2exp)
library(mthemer)
library(tidyquant)
library(timetk)
library(tsibble)



## Define constants ####
OUTPUT_FILE <- "output/financial_results.txt"
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
plot_rolling_causal_coeff <- function(dataset, cause, effect,
                                      ylab, span = 0.1){
  ## tibble character character character numeric -> plot
  ## plot the rolling causal coefficient through time

  dat_plot <- dataset %>%
    filter(str_detect(ticker, effect)) %>%
    mutate(ticker = factor(ticker,
                           levels = c(paste("psi_", cause, "_",
                                            effect, sep = ""),
                                      paste("psi_", effect, "_",
                                            cause, sep = ""))))

  cause_lab <- toupper(cause)
  effect_lab <- capitalize(effect)

  list_labels <- list(
    TeX(paste("$\\widehat{\\Psi}_{", cause_lab, "\\rightarrow",
              effect_lab, "}$")),
    TeX(paste("$\\widehat{\\Psi}_{", effect_lab, "\\rightarrow",
              cause_lab, "}$"))
  )

  dat_dummy <- dat_plot[1, ]
  dat_dummy$psi <- NaN

  plt <- ggplot() +
    geom_smooth(data = dat_plot,
                aes(x = date, y = psi, color = ticker),
                method = 'loess', span = span, se = FALSE,
                show.legend = FALSE) +
    geom_point(data = dat_dummy,
               aes(x = date, y = psi, color = ticker),
               shape = 15,
               size = 3) +
    mthemer() +
    scale_colour_manual(labels = list_labels,
                        values = tolPalette) +
    theme(legend.title=element_blank()) +
    xlab("Date") +
    ylab(ylab)

  return(plt)
}



plot_robustness_k <- function(dataset, cause, effect, ylab){
  ## tibble character character character -> plot
  ## plot the psi coefficient for varying k

  dat_plot <- dataset %>%
    filter(str_detect(ticker, effect)) %>%
    mutate(ticker = factor(ticker,
                           levels = c(paste("psi_", cause, "_",
                                            effect, sep = ""),
                                      paste("psi_", effect, "_",
                                            cause, sep = ""))))

  cause_lab <- toupper(cause)
  effect_lab <- capitalize(effect)

  list_labels <- list(
    TeX(paste("$\\widehat{\\Psi}_{", cause_lab, "\\rightarrow",
              effect_lab, "}$")),
    TeX(paste("$\\widehat{\\Psi}_{", effect_lab, "\\rightarrow",
              cause_lab, "}$"))
  )

  plt <- ggplot(dat_plot) +
    geom_ribbon(aes(x = root,
                    ymin = lb, ymax = ub,
                    fill = ticker),
                  alpha = .25) +
    geom_line(aes(x = root, y = psi_est, color = ticker),
              alpha = .9, size = 1) +
    geom_point(aes(x = root, y = psi_est, color = ticker),
               size = 3, shape = 21, fill = "white") +
    mthemer() +
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 16),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14)
    ) +
    scale_colour_manual(labels = list_labels,
                        values = tolPalette) +
    scale_fill_manual(labels = list_labels,
                        values = tolPalette) +
    theme(legend.title=element_blank()) +
    xlab(TeX("Fractional exponent of $k_n$")) +
    ylab(ylab) +
    guides(fill=FALSE)

  return(plt)
}




## Import dataset ####
dat_eurchf <- read_csv("data/financial_data/raw_dat_finance.csv",
                       col_types = cols(
                         date = col_date(format = "%Y-%m-%d")))




## Check model assumptions ####
mat <- dat_eurchf %>%  select(-date) %>% as.matrix()

shape_param <- matrix(nrow = NCOL(mat), ncol = 4)
colnames(shape_param) <- c("lt_shape", "lt_sd", "ut_shape", "ut_sd")
rownames(shape_param) <- colnames(mat)
nextreme <- 200
q <- (1 - nextreme / NROW(mat))


for (j in 1:NCOL(mat)){

  # upper tail
  tickr <- mat[, j]
  thres <- quantile(tickr, q, na.rm = T)
  out <- evd::fpot(tickr, thres, start = list(scale = 0.1, shape = 0))
  shape_param[j, 3] <- out$estimate["shape"]
  shape_param[j, 4] <- out$std.err["shape"]

  # lower tail
  tickr <- -mat[, j]
  thres <- quantile(tickr, q,  na.rm = T)
  out <- evd::fpot(tickr, thres, start = list(scale = 0.1, shape = 0))
  shape_param[j, 1] <- out$estimate["shape"]
  shape_param[j, 2] <- out$std.err["shape"]
}

sink(OUTPUT_FILE)
cat("--------------------------------------------------------------------------",
    "\n")
cat("MLE estimates of shape parameter xi (lower tail = lt, and upper tail = ut)",
    "\n")
cat("--------------------------------------------------------------------------",
    "\n")
print(shape_param)
cat("\n\n\n")
sink()

## Define crisis periods ####
event <- c("Event1", "Event2")
startDate <- c("2011-08-02", "2015-01-15")
endDate <- c("2011-12-02", "2015-05-15")

dates <- cbind(event, startDate, endDate) %>%
  as_tibble() %>%
  mutate(startDate = as.Date(startDate),
         endDate   = as.Date(endDate))


## Plot EURCHF timeseries ####
plt_eurchf <- ggplot() +
  mthemer() +
  xlab("Year") +
  ylab("EURCHF return") +
  scale_x_date(date_breaks = "1 year",
               date_labels = "%y") +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 18)) +
  geom_rect(data = dates, aes(xmin = startDate,
                              xmax = endDate,
                              ymin = -.15,
                              ymax = .1),
            fill = "red", alpha = .3) +
  geom_point(data = dat_eurchf, aes(x = date, y = eurchf),
             size = 3) +
  coord_cartesian(ylim=c(-0.15, .1))


ggsave("output/eurchf.pdf",
       plt_eurchf, width = 10, height = 6, units = c("in"))




## Causal discovery ####
# true DAG
true_dag <- rbind(c(0, 1, 1, 1),
                  rep(0, 4),
                  rep(0, 4),
                  rep(0, 4))

g <- igraph::graph_from_adjacency_matrix(true_dag)
igraph::V(g)$color <- "white"
# igraph::tkplot(g)

# data
mat <- dat_eurchf %>%
  select(-date) %>%
  as.matrix()

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
cat("EASE \n")
cat("--------------------------------------------------------------------------",
    "\n")
k <- 10
ease <- causal_discovery(mat, "ease", set_args = list(both_tails = T))
cat("SID:", causal_metrics(ll, ease)$sid, "\n")
g <- causal_tail_matrix(mat, both_tails = T, k = k)
colnames(g) <- rownames(g) <- colnames(mat); g
caus_ord <- ease(dat = mat, both_tails = T, k = k)
cat("Causal order:", colnames(mat)[caus_ord], "\n")
cat("\n\n\n")

# run ica_lingam
cat("--------------------------------------------------------------------------",
    "\n")
cat("ICA-LiNGAM \n")
cat("--------------------------------------------------------------------------",
    "\n")
ica_lingam <- causal_discovery(mat, "ica_lingam")
cat("SID:", causal_metrics(ll, ica_lingam)$sid, "\n")
cat("\n\n\n")

# run direct_lingam
cat("--------------------------------------------------------------------------",
    "\n")
cat("Pairwise LiNGAM \n")
cat("--------------------------------------------------------------------------",
    "\n")
direct_lingam <- causal_discovery(mat, "direct_lingam")
cat("SID:", causal_metrics(ll, direct_lingam)$sid, "\n")
cat("\n\n\n")

# run rank pc
cat("--------------------------------------------------------------------------",
    "\n")
cat("Rank PC \n")
cat("--------------------------------------------------------------------------",
    "\n")
pc_rank <- causal_discovery(mat, "pc_rank")
cat("SID:", causal_metrics(ll, pc_rank)$sid, "\n")
cat("\n\n\n")
sink()



## Plot rolling timeseries ####
window_size <- 1500
rolling_chf <- dat_eurchf %>%
  mutate(
    psi_eurchf_nestle = slide2_dbl( # CHF -> NESTLE
      .x = eurchf,
      .y = nestle,
      .f = causal_tail_coeff,
      k = k,
      .size = window_size),
    psi_nestle_eurchf = slide2_dbl( # NESTLE -> CHF
      .x = nestle,
      .y = eurchf,
      .f = causal_tail_coeff,
      k = k,
      .size = window_size),
    psi_eurchf_novartis = slide2_dbl( # CHF -> NOVARTIS
      .x = eurchf,
      .y = novartis,
      .f = causal_tail_coeff,
      k = k,
      .size = window_size),
    psi_novartis_eurchf = slide2_dbl( # NOVARTIS -> CHF
      .x = novartis,
      .y = eurchf,
      .f = causal_tail_coeff,
      k = k,
      .size = window_size),
    psi_eurchf_roche = slide2_dbl( # CHF -> ROCHE
      .x = eurchf,
      .y = roche,
      .f = causal_tail_coeff,
      k = k,
      .size = window_size),
    psi_roche_eurchf = slide2_dbl( # ROCHE -> CHF
      .x = roche,
      .y = eurchf,
      .f = causal_tail_coeff,
      k = k,
      .size = window_size)
  ) %>%
  select(date, starts_with("psi")) %>%
  drop_na() %>%
  gather(key = "ticker", value = "psi", -date)

plt_nestle <- plot_rolling_causal_coeff(rolling_chf, "eurchf", "nestle",
                                        TeX("$\\widehat{\\Psi}$"), span = 1e-2) +
  geom_rect(data = dates, aes(xmin = startDate, xmax = endDate, ymin = .3, ymax = 1),
            fill = "red", alpha = .3) +
  coord_cartesian(ylim=c(0.5, .975)) +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14)
  )

plt_novartis <- plot_rolling_causal_coeff(rolling_chf, "eurchf", "novartis",
                                          TeX("$\\widehat{\\Psi}$"), span = 1e-2) +
  geom_rect(data = dates, aes(xmin = startDate, xmax = endDate, ymin = .3, ymax = 1),
            fill = "red", alpha = .3) +
  coord_cartesian(ylim=c(0.5, .975)) +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14)
  )

plt_roche <- plot_rolling_causal_coeff(rolling_chf, "eurchf", "roche",
                                       TeX("$\\widehat{\\Psi}$"), span = 1e-2) +
  geom_rect(data = dates, aes(xmin = startDate, xmax = endDate, ymin = .3, ymax = 1),
            fill = "red", alpha = .3) +
  coord_cartesian(ylim=c(0.5, .975)) +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14)
  )


plt_grid <- plot_grid(plt_nestle,
                      plt_novartis,
                      plt_roche,
                      labels = c('A', 'B', 'C'),
                      label_size = 20,
                      nrow = 3)

ggsave("output/rolling_window_psi_coeff.pdf",
       plt_grid, width = 8, height = 10, units = c("in"))


## Robustness wrt to k ####
set.seed(1436) # Gutenberg's press

# helper function
compute_psi_tibble <- function(mat, v1_name, v2_name, k, root){
  ## numeric_matrix character character integer -> tibble
  ## produce a tibble with following columns
  ## k: number of observations in the tails
  ## root: fractional exponent of k
  ## v1_name: name of the variable in colnames(mat) that we condition on
  ## v2_name: name of the other variable in colnames(mat)
  ## psi: causal tail coefficient (considering both tails) of v1 -> v2
  col_names <- colnames(mat)
  idx_v1 <- which(col_names == v1_name)
  idx_v2 <- which(col_names == v2_name)
  submat <- mat[, c(idx_v1, idx_v2)]
  psi <- causal_tail_coeff(v1 = submat[, 1], v2 = submat[, 2], k = k)
  tbl <- tibble(k = k,
                root = root,
                v1_name = v1_name,
                v2_name = v2_name,
                psi = psi)
  return(tbl)
}

# data
mat <- dat_eurchf %>%
  select(-date) %>%
  as.matrix()

# varying k
len_root <- 21
n_subsample <- 1000
ticker_vec <- c("nestle", "novartis", "roche")
n_iter <- len_root * n_subsample
root_vec <- seq(.3, .7, length.out = len_root)

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

  # Prepare tibble
  tbl <- tibble(k = double(),
                v1_name = character(),
                v2_name = character(),
                psi = double(),
                psi_est = double())

  # Loop through tickers
  for (my_ticker in ticker_vec){
    psi1 <- compute_psi_tibble(mat, "eurchf", my_ticker, k, root)$psi
    psi2 <- compute_psi_tibble(mat, my_ticker, "eurchf", k, root)$psi
    tbl1 <- bind_cols(compute_psi_tibble(current_mat, "eurchf",
                                         my_ticker, k, root),
                      psi_est = psi1)
    tbl2 <- bind_cols(compute_psi_tibble(current_mat, my_ticker,
                                         "eurchf", k, root),
                      psi_est = psi2)
    tbl <- tbl %>% bind_rows(tbl1, tbl2)
  }

  # Return result
  return(tbl)
}
stopCluster(cl)

# Plot results
dat <- out_list %>%
  mutate(ticker = paste("psi_", v1_name, "_", v2_name, sep = "")) %>%
  select(root, ticker, psi, psi_est)

alpha <- .1

res <- dat %>%
  mutate(psi_diff = psi - psi_est) %>%
  group_by(root, ticker) %>%
  summarise(mean_psi = mean(psi),
            psi_est = mean(psi_est),
            se = sd(psi),
            lb = min(psi_est - quantile(psi_diff, probs = 1 - alpha/2), 1),
            ub = min(1, psi_est - quantile(psi_diff, probs = alpha/2))) %>%
  ungroup() %>%
  mutate(ticker = as.factor(ticker))


plt_nestle_rob <- plot_robustness_k(res, "eurchf", "nestle",
                  TeX("$\\widehat{\\Psi}$")) +
  ylim(c(0.4, 1.1))

plt_novartis_rob <- plot_robustness_k(res, "eurchf", "novartis",
                  TeX("$\\widehat{\\Psi}$")) +
  ylim(c(0.4, 1))

plt_roche_rob <- plot_robustness_k(res, "eurchf", "roche",
                  TeX("$\\widehat{\\Psi}$"))  +
  ylim(c(0.4, 1))

plt_grid <- plot_grid(plt_nestle_rob,
                      plt_novartis_rob,
                      plt_roche_rob,
                      labels = c('A', 'B', 'C'),
                      label_size = 20,
                      nrow = 3); plt_grid

ggsave("output/financial_k_robustness.pdf",
       plt_grid, width = 8, height = 10, units = c("in"))
