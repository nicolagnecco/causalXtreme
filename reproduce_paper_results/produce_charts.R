## Project: Causal discovery in heavy-tailed data
## Descrip: Produce charts of simulation study
## Authors: Nicola Gnecco, Nicolai Meinshausen, Jonas Peters, Sebastian Engelke


produce_charts <- function(sim0_file, sim1_file, sim2_file, sim3_file,
                           sim4_file){

  # Check whether the sim1_file and sim3_file are both demo or not
  if(dirname(sim1_file) != dirname(sim3_file)){
    stop("Simulation 1 and simulation 3 must be both demo or non-demo.")
  }

  # Set plotting theme
  theme_set(theme_bw() +
              theme(plot.background = element_blank(),
                    legend.background = element_blank()))

  ### CONSTANTS ####
  SIMULATION_K <- sim0_file
  SIMULATION <- sim1_file
  SIMULATION_TIME <- sim2_file
  SIMULATION_LINGAM <- sim3_file
  SIMULATION_RANKPC <- sim4_file
  SID_PLOTS <- "output/SID.pdf"
  SID_PLOTS_LIGHT1 <- "output/SID_alpha_25.pdf"
  SID_PLOTS_LIGHT2 <- "output/SID_alpha_35.pdf"
  TIME_PLOTS <- "output/TIME.pdf"
  SHD_PLOTS <- "output/SHD.pdf"
  LINGAM_REMOVE_BULK <- "output/table_lingam.txt"
  RANKPC_PLOTS <- "output/RANKPC_SIG.pdf"
  SID_PLOTS_POSTER <- "output/SID_POSTER.pdf"
  SETTING1 <- "output/setting1.pdf"
  SETTING2 <- "output/setting2.pdf"
  SETTING3 <- "output/setting3.pdf"
  SETTING4 <- "output/setting4.pdf"
  HTSCM <- "output/heavytailedSCM.pdf"
  HTSCMx <- "output/heavytailedSCMx.pdf"
  HTSCMy <- "output/heavytailedSCMy.pdf"
  SID_KVARYING <- "output/SID_kvarying.pdf"
  tolPalette <- c(tolBlue = "#4477AA",
                  tolRed = "#EE6677",
                  tolGreen = "#228833",
                  tolYellow = "#CCBB44",
                  tolCyan = "#66CCEE",
                  tolPurple = "#AA3377",
                  tolGrey = "#BBBBBB")


  ### FUNCTION DEFINITIONS ####
  change_method_names <- function(dat){
    ## tibble -> tibble
    ## changes names of methods

    dat %>%
      rowwise() %>%
      mutate(method = if(method == "pc_rank"){
        "Rank PC"
      } else if(method == "direct_lingam"){
        "Pairwise LiNGAM"
      } else if(method == "ica_lingam"){
        "ICA-LiNGAM"
      }else if(method == "ease"){
        "EASE"
      }else if(method == "random"){
        "Random"
      }else {
        method
      }) %>%
      mutate(method = factor(method,
                             levels = c("EASE", "Pairwise LiNGAM",
                                        "ICA-LiNGAM", "Rank PC", "Random")
      )) %>%
      ungroup()

  }

  prepare_data <- function(dat, type = c("multiple", "single")){
    ## tibble character -> tibble
    ## reshape tibble to print it

    type <- match.arg(type)

    # Filter out unneeded observations
    res <- dat %>%
      filter(method %in% c("ease", "ica_lingam", "direct_lingam",
                           "pc_rank", "random"))

    # Remove unneeded variables
    res <- res %>%
      select(-distr, -prob_connect)

    # Create dictionary with settings
    settings <- list(
      has_confounder = c(FALSE, TRUE),
      is_nonlinear = c(FALSE, TRUE),
      has_uniform_margins = c(FALSE, TRUE))
    settings <- settings %>% cross_df() %>%
      bind_cols(settings = c(1:3, 5, 4, 6:8))

    # Left join
    res <- res %>%
      left_join(settings, by = c("has_confounder", "has_uniform_margins",
                                 "is_nonlinear")) %>%
      filter(settings %in% 1:4)

    # Aggregate values
    if (type == "multiple"){
      res <- res %>%
        group_by(settings, n, p, method) %>%
        summarise(structHamming_dist = mean(shd, na.rm = T),
                  structInterv_dist = mean(sid, na.rm = T),
                  time = mean(time, na.rm = T))
    } else if (type == "single") {
      res <- res %>%
        group_by(n, p, method) %>%
        summarise(structHamming_dist = mean(shd, na.rm = T),
                  structInterv_dist = mean(sid, na.rm = T),
                  time = mean(time, na.rm = T))
    } else {
      stop("Please provide one of: single, multiple.")
    }

    # Change names of methods
    res <- change_method_names(res)

    # Add labels to variables for facets
    n_label <- paste("n =", res$n)
    res$n_label <- factor(n_label,
                          levels =
                            unique(n_label[order(res$n)])) # reoder levels

    if (type == "multiple"){
      setting_label <- paste("Setting", res$settings)
      res$setting_label <- factor(setting_label,
                                  levels =
                                    unique(setting_label
                                           [order(res$settings)]))  # reoder levels
    }
    return(res)
  }

  plot_results <- function(dat, y, type = c("multiple", "single"),
                           palette = c(tolBlue = "#4477AA",
                                       tolRed = "#EE6677",
                                       tolGreen = "#228833",
                                       tolYellow = "#CCBB44",
                                       tolCyan = "#66CCEE",
                                       tolPurple = "#AA3377",
                                       tolGrey = "#BBBBBB")){
    ## tibble character character character_vector -> ggplot
    ## produce ggplot of results according to the measure y

    type <- match.arg(type)

    # Define colorblind-friendly palette
    palette <- unname(palette)

    # Prepare y-labels
    y_lab <- if(y == "structHamming_dist"){
      "Structural Hamming Distance"
    }else if(y == "structInterv_dist"){
      "Structural Intervention Distance"
    }else if(y == "time"){
      "Time"
    }else{
      stop(paste("The column", y, "is not in the dataset."))
    }

    # Prepare axes
    x <-  "p"
    color <-  "method"
    facet_x <-  "n_label"
    facet_y <-  "setting_label"

    # Create plot
    if(type == "multiple"){
      g <- ggplot(data = dat, aes_string(x = x, y = y, color = color)) +
        geom_line(size = 1, alpha = .5) +
        geom_point(size = 3) +
        ylab(y_lab) +
        facet_grid(reformulate(facet_x, facet_y)) +
        scale_colour_manual(values = palette) +
        theme(legend.title=element_blank())

      return(g)

    } else if(type == "single"){
      # One setting only
      g <- ggplot(data = dat, aes_string(x = x, y = y, color = color)) +
        geom_line(size = 1, alpha = .5) +
        geom_point(size = 3) +
        ylab(y_lab) +
        facet_grid(cols = vars(n_label)) +
        scale_colour_manual(values = palette) +
        theme(legend.position="bottom") +
        theme(legend.title=element_blank())

      return(g)
    } else{
      stop("Please provide one of: single, multiple.")
    }
  }

  ### FIGURES FOR PAPER ####

  ## Performance of methods (alpha = 1.5) ----
  # Import data
  dat <- read_rds(SIMULATION) %>%
    filter(tail_index == 1.5)

  # Transform data
  res <- prepare_data(dat)

  # SID plots multiple settings
  g <- plot_results(res, "structInterv_dist", "multiple") +
    theme(strip.text = element_text(size = 16),
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 16),
          legend.text = element_text(size = 14)
    )

  ggsave(SID_PLOTS, g, width = 10, height = 10, units = c("in"))


  ## Performance of methods (alpha = 2.5) ----
  # Import data
  dat <- read_rds(SIMULATION) %>%
    filter(tail_index == 2.5)

  # Transform data
  res <- prepare_data(dat)

  # SID plots multiple settings
  g <- plot_results(res, "structInterv_dist", "multiple") +
    theme(strip.text = element_text(size = 16),
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 16),
          legend.text = element_text(size = 14)
    )

  ggsave(SID_PLOTS_LIGHT1, g, width = 10, height = 10, units = c("in"))


  ## Performance of methods (alpha = 3.5) ----
  # Import data
  dat <- read_rds(SIMULATION) %>%
    filter(tail_index == 3.5)

  # Transform data
  res <- prepare_data(dat)

  # SID plots multiple settings
  g <- plot_results(res, "structInterv_dist", "multiple") +
    theme(strip.text = element_text(size = 16),
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 16),
          legend.text = element_text(size = 14)
    )

  ggsave(SID_PLOTS_LIGHT2, g, width = 10, height = 10, units = c("in"))


  ## Plot computational complexity ----
  dat <- read_rds(SIMULATION_TIME) %>%
    prepare_data(type = "single") %>%
    mutate(time = log10(as.numeric(time)))

  g <- plot_results(dat, "time", type = "single") +
    ylab(TeX("log$_{10}$(Time)")) +
    theme(strip.text = element_text(size = 16),
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 16),
          legend.text = element_text(size = 14)
    )

  ggsave(TIME_PLOTS, g, width = 8, height = 3, units = c("in"))

  ## Table computational complexity ----
  dat <- read_rds(SIMULATION_TIME) %>%
    change_method_names() %>%
    group_by(method, n, p) %>%
    summarise(avg_time = round(mean(time), 2)) %>%
    spread(n, avg_time) %>%
    mutate(dim = paste("$p = ", p, "$", sep = "")) %>%
    select(method, dim, everything(), -p) %>%
    rename("\ " = method, "\ \ " = dim)


  kable(dat,
        format = "html",
        longtable = F,
        escape = F,
        booktabs = T,
        caption = "to add")  %>%
    kable_styling(latex_options = c("repeat_header", "hold_position")) %>%
    column_spec(1, bold = T) %>%
    collapse_rows(columns = c(1, 2)) %>%
    add_header_above(c(" " = 2, "Sample size" = 3))


  ## Table Lingam experiment ----
  dat_lingam <- read_rds(SIMULATION_LINGAM) %>%
    mutate(type = "tail")

  dat <- read_rds(SIMULATION) %>%
    mutate(type = "all") %>%
    filter(method %in% c("ica_lingam", "direct_lingam"),
           id %in% dat_lingam$id) %>%
    bind_rows(dat_lingam) %>%
    filter(tail_index == 3.5) %>%
    group_by(type, method, n, p) %>%
    summarise(SID = round(mean(sid, na.rm = T), 3),
              SE = round(sd(sid, na.rm = T) / sqrt(n()), 3)) %>%
    gather(key = "metric", value = "res", -type, -method, -n, -p) %>%
    unite(type_metric, c("type", "metric")) %>%
    spread(key = type_metric, value = res) %>%
    arrange(desc(method)) %>%
    mutate(Dimension = paste("$p = ", p, "$", sep = "")) %>%
    ungroup() %>%
    select(-n, -p) %>%
    select(method, Dimension, all_SID, all_SE, tail_SID, tail_SE) %>%
    change_method_names()

  kable(dat,
        format = "html",
        longtable = F,
        escape = F,
        booktabs = T,
        caption = "to add",
        col.names = c("", "",
                      "\\emph{SID}", "\\emph{SE}",
                      "\\emph{SID}", "\\emph{SE}"))  %>%
    kable_styling(latex_options = c("repeat_header", "hold_position")) %>%
    column_spec(1, bold = T) %>%
    collapse_rows(columns = c(1, 2)) %>%
    add_header_above(c(" " = 2, "\\emph{Full dataset}" = 2,
                       "\\emph{Keep tails}" = 2),
                     escape = F)

  sink(file = LINGAM_REMOVE_BULK)
  cat("--------------------------------------------------------------------------",
      "\n")
  cat("Computing time in seconds when considering all data ('all')",
      "\n",
      "or only extreme observations ('tail'), for ICA-LiNGAM and Pairwise LiNGAM",
      "\n")
  cat("--------------------------------------------------------------------------",
      "\n")
  print(dat)
  sink()
  closeAllConnections()

  ## Plot Rank PC experiment----
  dat <- read_rds(SIMULATION_RANKPC)
  dat2 <- dat %>% filter(method == "random") %>%
    select(-arg_name, -arg_value) %>%
    left_join(tibble(method = "random",
                     arg_name = rep("alpha", 4),
                     arg_value = c(5e-4, 5e-3, 5e-2, 5e-1)), by = 'method')
  dat3 <- dat %>% filter(method != "random") %>%
    select(-arg_name, -arg_value, arg_name, arg_value) %>%
    bind_rows(dat2) %>%
    rowwise() %>%
    mutate(alpha = factor(arg_value),
           method = if(method == "pc_rank"){
             "Rank PC"
           } else if(method == "random"){
             "Random"},
           method = factor(method),
           p = factor(paste("p =", p))) %>%
    group_by(method, p, alpha) %>%
    summarise(structInterv_dist = mean(sid, na.rm = T))

  g <- ggplot(data = dat3, aes(x = alpha,
                              y = structInterv_dist,
                              col = method,
                              group = method)) +
    geom_point(size = 3)  +
    stat_summary(fun=sum, geom="line", size = 1, alpha = .5) +
    ylab("Structural Intervention Distance") +
    xlab("Significance level") +
    ylim(0, 1) +
    facet_grid(rows = vars(p), scales = "free_y") +
    scale_color_manual(values = unname(tolPalette)[4:5]) +
    theme(legend.title=element_blank()); g

  ggsave(RANKPC_PLOTS, g, width = 8, height = 4, units = c("in"))


  ## Robustness of k across different settings wrt SID ----
  # Import data
  dat <- read_rds(SIMULATION_K) %>% unnest(cols = c(data)) %>%
    filter(has_confounder == F, is_nonlinear == F,
           has_uniform_margins == F)

  res <- dat %>%
    group_by(tail_index, root) %>%
    summarise(mean_sid = mean(sid),
              N = n(),
              se = sd(sid)/sqrt(N)) %>%
    ungroup() %>%
    mutate(tail_index = as.factor(tail_index))

  pd <- position_dodge(0.00)

  g <- ggplot(res) +
    geom_line(aes(x = root, y = mean_sid, color = tail_index),
              alpha = .7, size = 1, position = pd) +
    geom_point(aes(x = root, y = mean_sid, color = tail_index),
               size = 3, shape = 21, fill = "white", position = pd) +
    scale_color_manual(values = unname(tolPalette)) +
    scale_fill_manual(values = unname(tolPalette)) +
    xlab(TeX("Fractional exponent of $k_n$")) +
    ylab("Structural Intervention Distance") +
    labs(color = TeX("Tail index $\\alpha$")) +
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 16),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14)
    ) + guides(fill=FALSE); g

  ggsave(SID_KVARYING, g, width = 7.5, height = 5, units = c("in"))


  ### FIGURES FOR TALK ####

  ## Performance of methods
  # Import data
  dat <- read_rds(SIMULATION) %>% unnest()

  # Transform data
  res <- prepare_data(dat)

  # SID plots Setting 1
  setting <- 1
  g <- plot_results(res %>% filter(settings == setting),
                    "structInterv_dist", "single") +
    theme(strip.text = element_text(size = 16),
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 16),
          legend.text = element_text(size = 14)
    ) +
    ylab("SID")

  ggsave(SETTING1, g, width = 10, height = 3.5, units = c("in"))

  # SID plots Setting 2
  setting <- 2
  g <- plot_results(res %>% filter(settings == setting),
                    "structInterv_dist", "single") +
    theme(strip.text = element_text(size = 16),
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 16),
          legend.text = element_text(size = 14)
    ) +
    ylab("SID")

  ggsave(SETTING2, g, width = 10, height = 3.5, units = c("in"))

  # SID plots Setting 3
  setting <- 3
  g <- plot_results(res %>% filter(settings == setting),
                    "structInterv_dist", "single") +
    theme(strip.text = element_text(size = 16),
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 16),
          legend.text = element_text(size = 14)
    ) +
    ylab("SID")

  ggsave(SETTING3, g, width = 10, height = 3.5, units = c("in"))

  # SID plots Setting 4
  setting <- 4
  g <- plot_results(res %>% filter(settings == setting),
                    "structInterv_dist", "single") +
    theme(strip.text = element_text(size = 16),
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 16),
          legend.text = element_text(size = 14)
    ) +
    ylab("SID")

  ggsave(SETTING4, g, width = 10, height = 3.5, units = c("in"))


  ## Heavy-tailed SCM example
  set.seed(132)
  n <- 1e4
  df <- 5
  thresy <- 6
  thresx <- 5
  e1 <- rt(n, df)  # noise for variable 1
  e2 <- rt(n, df)  # noise for variable 2
  X1 <- e1
  X2 <-  X1 + e2
  colr <- tolPalette[1]

  # plot results
  df <- data.frame(X1 = X1, X2 = X2)

  g <- ggplot(df, aes(x = X1, y = X2)) +
    geom_point() +
    xlim(0, 6) +
    ylim(0, 8) +
    xlab(TeX("$X_1$")) +
    ylab(TeX("$X_2$")) +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 18)
    )

  ggsave(HTSCM, g, width = 5, height = 5, units = c("in"))

  g <- ggplot(df, aes(x = X1, y = X2)) +
    geom_vline(xintercept = thresx,  color = colr, size = 1) +
    geom_point() +
    xlim(0, 6) +
    ylim(0, 8) +
    xlab(TeX("$X_1$")) +
    ylab(TeX("$X_2$")) +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 18)
    )

  ggsave(HTSCMx, g, width = 5, height = 5, units = c("in"))

  g <- ggplot(df, aes(x = X1, y = X2)) +
    geom_hline(yintercept = thresy, color = colr, size = 1) +
    geom_point() +
    xlim(0, 6) +
    ylim(0, 8) +
    xlab(TeX("$X_1$")) +
    ylab(TeX("$X_2$")) +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 18)
    )

  ggsave(HTSCMy, g, width = 5, height = 5, units = c("in"))

  ### FIGURES FOR POSTER ####
  # Import data
  dat <- read_rds(SIMULATION) %>% unnest()

  # Transform data
  res <- prepare_data(dat) %>%
    filter(setting_label != "Setting 3") %>%
    mutate(setting_label = fct_recode(setting_label, 'Setting 3' = 'Setting 4'))


  # SID plots multiple settings
  g <- plot_results(res, "structInterv_dist", "multiple") +
    theme(strip.text = element_text(size = 16),
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 16),
          legend.text = element_text(size = 14)
    )

  ggsave(SID_PLOTS_POSTER, g, width = 15, height = 10, units = c("in"))
}
