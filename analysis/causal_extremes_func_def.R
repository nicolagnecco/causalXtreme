## Title: Function definitions for Causal Extremes project
## Authors: Sebastian Engelke, Nicolai Meinshausen, Jonas Peters, Nicola Gnecco

## Libraries
library(gtools)
library(pcalg)
library(stringr)
library(SID)
source("randomDAG.R", chdir = T)
source("randomB.R", chdir = T)
source("computeCausOrder.R", chdir = T)
source("greedy_perm_search.R", chdir = T)
source("fast_perm_search.R", chdir = T)
source("full_perm_search.R", chdir = T)
source("random_perm_search.R", chdir = T)

# !!! to check: dividing by n and look at quantile > than some threshold (gammaCoeff, unif, shifted_hockeystick)

## General functions ####
simulateData <- function(n, p, distr, df = 1.5, adj.mat = matrix(NA),
                         has.confounder = FALSE, scale = TRUE, sig_w = 0.5, nonlinear = FALSE,
                         sparsity = 5/(p-1), sparsity.confounder = 2/(3*(p-1)), transform.margins = FALSE){
  ## integer integer character numeric_matrix boolean boolean numeric boolean numeric numeric -> list
  ## produces a list with:
  ## - data: simulated data from a random causal DAG
  ## - adj.mat: weighted adjacency matrix of a random causal DAG
  ## given the parameters:
  ## - n: number of observations
  ## - p: number of variables
  ## - distr: distribution of the noise
  ## - df: degrees of freedom of distribution/tail index
  ## - adj.mat: weighted adjacency matrix of a caual DAG (optional)
  ## - has.confounder: introduce confounders?
  ## - scale: marginally scale each variable?
  ## - sig_w: proportion of signal versus noise
  ## - nonlinear: marginally transform each variable?
  ## - sparsity: probability to connect each node to one of its descendants
  ## - sparsity.confounder: probability that an unordered pair of nodes has one confounder
  ## ASSUME:
  ## 1. distr is one of:
  ##      - 'student.t': t-Student
  ##      - 'log.norm  : log-Normal
  ##      - 'norm'     : Normal
  ##      - 'unif'     : Uniform
  ## 2. adj.matr has non-negative weights
  ## 3. adj.matr has no cycles
  ## 3. sig_w is between 0 and 1

  ## helpers ####
  addRandomConfounders <- function(N, probConfouder, conf_w, alpha){
    ## numeric_matrix numeric numeric df -> numeric_matrix
    ## produces a noise matrix where each pair of variables (i, j) in N
    ## is confounded with probability probConfounder;
    ## conf_w is the signal proportion of the confounder versus the noise;
    ## alpha is the tail index of the distribution
    ## ASSUME:
    ## 1. conf_w is between 0 and 1

    # Check inputs
    if(conf_w > 1 | conf_w < 0){stop("The signal weight must be between 0 and 1.")}

    # Get variables
    n <- NROW(N)
    p <- NCOL(N)

    # Consider all possible pairs of nodes in G
    M <- which(upper.tri(diag(p)), arr.ind=TRUE)
    nPossibleConf  <- NROW(M)

    # Choose the confounders randomly
    numberConf  <- rbinom(n=1, size=nPossibleConf, prob=probConfounder)
    Confounders <- sample(x = 1:nPossibleConf, size = numberConf, replace = FALSE)
    Confounders <- M[Confounders, ]

    if(numberConf == 0){
      return(N)
    }

    # Sample the confounder weights
    S <- matrix(0, nrow = numberConf, ncol = p) # preallocate matrix with confounder-variable pairs
    confvar_pairs <- cbind(rep(1:numberConf, each=2), c(t(Confounders))) # confounder-variable pairs
    S[confvar_pairs] <- 1 # populate matrix with confounder-variable pairs
    S <- abs(randomB(t(S), lB = lB, uB = uB)) # sample random weights of confounders

    # Rescale confounder weights
    S <- scale_adjmat(S, conf_w, alpha)
    noise_w <- 1 - conf_w
    N <- scale_noise(N, S, noise_w, alpha)

    # Generate confounders
    C <- simulateNoise(n, numberConf, distr, alpha)

    # Return confounded noise
    return(N + C %*% S)
  }

  scale_adjmat <- function(adj.mat, sig_w, alpha){
    ## numeric_matrix numeric numeric -> numeric_matrix
    ## produces a weighted adjacency matrix with rescaled beta

    s <- apply(adj.mat^alpha, 2, sum)
    s[s != 0] <- 1 / (s[s != 0]) * sig_w
    S <- diag(s^(1/alpha))
    adj.mat %*% S
  }

  scale_noise <- function(N, adj.mat, noise_w, alpha){
    ## numeric_matrix numeric_matrix numeric numeric -> numeric_matrix
    ## produce a scaled version of the noise matrix N
    p <- NCOL(adj.mat)
    is_source <- apply(adj.mat, 2, sum) == 0
    d <- numeric(p)
    d[is_source] <- 1
    d[!is_source] <- noise_w^(1/alpha)
    D <- diag(d)
    (N %*% D)
  }

  simulateNoise <- function(n, p, distr, df){
    ## integer integer character numeric -> numeric_matrix
    ## produce the noise matrix from the distribution distr where df is the
    ## number of degrees of freedom.
    ## ASSUME:
    ## 1. distr is one of:
    ##      - 'student.t': t-Student
    ##      - 'log.norm'  : log-Normal
    ##      - 'norm'     : Normal
    ##      - 'unif'     : Uniform

    switch(distr,
           "student.t" = {

             N <- array(rt(n * p, df=df), dim = c(n, p))

           },
           "log.norm" = {

             N <- array(rlnorm(n * p), dim = c(n, p))

           },
           "norm" = {

             N <- array(rnorm(n * p), dim = c(n, p))
           },
           "unif" = {

             N <- array(runif(n * p), dim = c(n, p))
           },
           stop("Wrong distribution. Enter one of 'student.t', 'log.norm', 'norm', 'unif'."))

    return(N)
  }

  nonlinear_SEM <- function(adj.mat, noise){
    ## numeric.matrix numeric.matrix -> numeric.matrix
    ## produce nonlinear SEM given
    ## - adj.mat: weighted adjacency matrix
    ## - noise: (n x p) matrix with noise entries

    ## helpers ####
    shifted_hockeystick <- function(v, q=0.8){
      ## numeric.vector -> numeric.vector
      ## produces vector that keeps only the entries of v > q-quantile

      n <- length(v)
      v[which(rank(v) <= ceiling(n * q))] <- 0
      return(v)
    }
    ## function body ####
    n <- NROW(noise)
    p <- NROW(adj.mat)
    nodes <- computeCausOrder(adj.mat)

    X <- matrix(0, nrow = n, ncol = p)
    X.transf <- matrix(0, nrow = n, ncol = p)
    for(i in nodes){
      betas <- adj.mat[, i]
      eps <- noise[, i]
      X[, i] <- X.transf %*% betas + eps
      X.transf[, i] <- shifted_hockeystick(X[, i])
    }

    return(X.transf)
  }

  unif <- function(x){
    ## numeric_v -> numeric_v
    ## produce a rescaled vector with uniform margins between 0 and 1
    n <- length(x)
    rank(x)/n
  }

  ## function body ####
  # Check inputs
  if(sig_w > 1 | sig_w < 0){stop("The signal weight must be between 0 and 1.")}
  if(!missing(adj.mat)){
    if(min(adj.mat, na.rm = T) < 0){stop("The weighted adjacency matrix must have non-negative weights.")}
  }
  if(distr == "norm"){
    df <- 2 # if the distribution is Gaussian, set the tail index to 2
  }
  if(sparsity > 1/2){
    sparsity <- 1/2
  }

  # Set parameters
  probConnect <- sparsity
  lB <- 0.1
  uB <- 0.9
  probConfounder <- sparsity.confounder
  conf_w <- 0.5
  distr <- distr
  alpha <- df

  # Prepare adjacency matrix
  if (missing(adj.mat)){
    adj.mat <- randomDAG(p, probConnect, sparse = FALSE)
    adj.mat <- abs(randomB(t(adj.mat), lB = lB, uB = uB))
  }

  # Scale adjacency matrix?
  if(scale == TRUE){
    adj.mat <- scale_adjmat(adj.mat, sig_w, alpha)
  }

  # Simulate noise
  N <- simulateNoise(n, p, distr, df)

  # Add confounders?
  if(has.confounder == TRUE){
    N <- addRandomConfounders(N, probConfounder, conf_w, alpha)
  }

  # Scale noise?
  if(scale == TRUE){
    noise_w <- 1 - sig_w
    N <- scale_noise(N, adj.mat, noise_w, alpha)
  }
  # Simulate data
  if(nonlinear == FALSE){
    B <- diag(p) - t(adj.mat)
    X <- t(solve(B, t(N)))
  } else {
    X <- nonlinear_SEM(adj.mat, N)
  }


  # Marginally transform each variable?
  if(transform.margins == TRUE){
    X <- apply(X, 2, unif)
  }

  # Return list
  ll <- list()
  ll$data    <- X
  ll$adj.mat <- adj.mat
  return(ll)

}

causalDiscovery <- function(gamma, method, adj.mat = matrix(NA), dat = matrix(NA)){
  ## numeric_matrix character numeric_matrix numeric_matrix -> list
  ## produces a list given the gamma matrix, a method, one adjacency matrix
  ## (if the ground truth is known) and the matrix with data (when method = "lingam").
  ## The list is made of:
  ## - order
  ## - ancestral_dist
  ## - structInterv_dist
  ##
  ## ASSUME:
  ## 1. method is one of:
  ##      - "fast": uses delta, tries to maximize score upper.tri(delta)
  ##      - "full": uses delta, maximizes score upper.tri(delta)
  ##      - "greedy": uses delta, tries to maximize score upper.tri(delta)
  ##      - "lingam": uses data
  ##      - "maxmin": uses delta, tries to maximize score upper.tri(delta)
  ##      - "minimax": uses gamma
  ##      - "oracle": uses adjacency matrix
  ##      - "pc": uses data
  ##      - "random": uses delta
  ## 2. adj.mat has non-negative weights
  ## 3. adj.matr has no cycles

  # Check inputs
  if(!missing(adj.mat)){
    if(min(adj.mat, na.rm = T) < 0){stop("The weighted adjacency matrix must have non-negative weights.")}
  }
  if(method=="lingam" & missing(dat)){
    stop("Please provide a numeric matrix with the data.")
  }

  # build delta matrix
  delta <- gamma - t(gamma)

  # count number of variables
  p <- NROW(delta)

  # set up list
  restmp <- list(order=NA, ancestral_dist=NA, structInterv_dist = NA)

  # compute causal order
  switch(method,
         "fast" = {
           restmp$order <- fast_perm_search(delta, silent = TRUE, mode = "sum")$order
           est_adj.mat <- perm2adjmat(restmp$order)

         },
         "full" = {
           if(p < 10){
             restmp$order <- full_perm_search(delta, silent = TRUE)$order
             est_adj.mat <- perm2adjmat(restmp$order)

           }
         },
         "greedy" = {
           restmp$order <- greedy_perm_search(delta, silent = TRUE)$order
           est_adj.mat <- perm2adjmat(restmp$order)

         },
         "lingam" = {
           out <- lingam_search(dat)
           restmp$order <- out$order
           est_adj.mat <- out$adj.mat

         },
         "maxmin" = {
           restmp$order <- fast_perm_search(delta, silent = TRUE, mode = "maxmin")$order
           est_adj.mat <- perm2adjmat(restmp$order)

         },
         "minimax" = {
           restmp$order <- minimax_search(gamma)
           est_adj.mat <- perm2adjmat(restmp$order)

         },
         "oracle" = {
           restmp$order <- oracle_search(adj.mat)
           est_adj.mat <- perm2adjmat(restmp$order)

         }, "pc" = {
           out <- pc_search(dat)
           restmp$order <- out$order
           est_adj.mat <- out$adj.mat

         }, "random" = {
           restmp$order <- random_perm_search(delta)$order
           est_adj.mat <- perm2adjmat(restmp$order)

         },
         stop("Wrong method. Enter one of 'fast', 'full', 'greedy', 'lingam', 'maxmin', 'minimax', 'oracle', 'random'."))

  # compute ancestral distance (if adj.mat is available)
  if (!missing(adj.mat)){
    restmp$ancestral_dist <- ancestral_distance(adj.mat, restmp$order)
  }

  # compute structural intervention distance (if adj.mat is available)
  if (!missing(adj.mat)){
    restmp$structInterv_dist <- structInterv_distance(adj.mat, est_adj.mat)
  }
  # return list
  return(restmp)
}

causalDiscovery_wrapper <- function(alpha, p, n, distr,
                                    confounder, nonlinear, transform.margins, method_tbl){
  ## numeric integer integer character boolean boolean boolean tibble -> list
  ## produces a tibble with causal discovery performances of the methods listed
  ## in the tibble method_tbl given the parameters:
  ## - alpha: tail index
  ## - p: number of variables
  ## - n: number of observations
  ## - distr: distribution of the noise
  ## - confounder: introduce confounders?
  ## - nonlinear: introduce non-linearity into the SEM?
  ## - transform.margins: marginally transform each variable?

  # simulate data
  X <- simulateData(n = n, p = p, df = alpha, distr = distr,
                    has.confounder = confounder, nonlinear = nonlinear, transform.margins = transform.margins)

  # estimate gamma
  gamma <- gammaMatrix(X$data, k = floor(sqrt(n)))

  # perform causal discovery
  l <- method_tbl[[1]] %>%
    map(function(m){
      time1 <- Sys.time()
      restmp <- causalDiscovery(gamma, m, X$adj.mat, X$data)
      ll <- list()
      ll$order <- restmp$order
      ll$ancestral_dist <- restmp$ancestral_dist
      ll$structInterv_dist <- restmp$structInterv_dist
      ll$iscorrect <- checkCausOrder(restmp$order, X$adj.mat)
      ll$time <- Sys.time() - time1
      return(ll)
    })

  # convert into tibble
  tibble(order = map(l, "order"),
         ancestral_dist = map_dbl(l, "ancestral_dist"),
         structInterv_dist = map_dbl(l, "structInterv_dist"),
         iscorrect = map_lgl(l, "iscorrect"),
         time = map_dbl(l, "time"))

}


## Causal discovery functions ####
minimax_search <- function(Gamma){
  ## numeric_matrix -> numeric_vector
  ## produces a causal order performing minimax search on the given Gamma matrix

  # set up variables
  d <- dim(Gamma)[2]
  GammaOrig <- Gamma
  diag(Gamma) <- NA

  ## run minimax
  current.order <- add <- which.min(apply(Gamma, 2 , max, na.rm=TRUE))
  for (k in 2:d){
    Gamma[add, ] <- NA
    avail <- (1:d)[-current.order]

    add <- if(k<d) avail[which.min(apply(Gamma, 2, max, na.rm=TRUE)[avail])] else avail
    current.order <- c(current.order, add)
  }

  order <- current.order

  return(order)
}

oracle_search <- function(A){
  ## numeric_matrix <- numeric_vector
  ## produces a true causal order from the given weighted adjacency matrix A

  G <- (A != 0) * 1
  order <- computeCausOrder(G)
  return(order)

}

lingam_search <- function(dat){
  ## numeric_matrix -> list
  ## produces a list given dataset dat. The list is made of:
  ## - adj.mat: estimated adjacency matrix by LiNGAM
  ## - order: causal order obtained from adj.mat

  out <- tryCatch(
    {
      lingam.output <- lingam(X=dat)
      Bpruned <- lingam.output$Bpruned
      estim.adj.mat <- (t(Bpruned) != 0) * 1
      order <- computeCausOrder(G = estim.adj.mat)
      l <- list(adj.mat = estim.adj.mat, order = order)

      return(l)
    },
    error = function(e){
      l <- list(adj.mat = NA, order = NA)

      return(l)
    })

  return(out)
}

pc_search <- function(dat){
  ## numeric_matrix numeric_matrix -> list
  ## produces a list given a dataset dat. The list is made of:
  ## - adj.mat: a CPDAG
  ## - order: a causal order obtained by removing cycles from adj.mat

  n <- NROW(dat)
  p <- NCOL(dat)
  suffStat <- list(C = cor(dat), n = n)

  out <- tryCatch({
    pc.fit <- pc(suffStat = suffStat, indepTest = gaussCItest, p = p, alpha = 5e-2,
                 u2pd = "retry", skel.method = "stable")
    cpdag <- as(pc.fit@graph, "matrix")
    dag <- cpdag * (t(cpdag) == 0)
    order <- computeCausOrder(dag)
    l <- list(adj.mat = cpdag, order = order)

    return(l)
  }, error = function(e){
    l <- list(adj.mat = NA, order = NA)
    return(l)
  })

  return(out)
}

## Extremal coefficient (i.e., Gamma and Delta) functions ####
deltaMatrix <- function(dat, k = floor(0.5 * sqrt(n))){
  ## numeric_matrix -> numeric_matrix
  ## produces the delta matrix for the given raw data, given the number of
  ## upper order statistics n

  n <- NROW(dat) # number of observations
  gamma.matrix <- gammaMatrix(dat, k)
  delta.matrix <- gamma.matrix - t(gamma.matrix)
  # delta.matrix[is.na(delta.matrix)] <- 0
  return(delta.matrix)

}

gammaMatrix <- function(dat, k = floor(0.5 * sqrt(n))){
  ## numeric_matrix integer -> numeric_matrix
  ## produces the matrix with gamma coefficients, given the raw data and
  ## the number of upper order statistics

  ## helpers ####
  gammaMatrix.helper <- function(i, j, v1, v2, k){
    ## integer integer numeric_vector numeric_vector integer -> numeric
    ## produces gamma coefficient of v1 and v2, given their position in the matrix
    ## and the number of upper order statistics

    if(i == j){NA}
    else{gammaCoeff(v1, v2, k)}
  }

  ## function body ####
  n   <- NROW(dat) # number of observations
  l   <- NCOL(dat) # number of variables
  nms <- colnames(dat)

  # compute gamma for all combinations of indices
  g <- sapply(1:l,
              function(j){
                sapply(1:l,
                       function(i){
                         gammaMatrix.helper(i, j, dat[, i], dat[, j], k)
                       })
              })

  colnames(g) <- rownames(g) <- nms

  # return data
  return(g)
}

gammaCoeff <- function(v1, v2, k = floor(0.5 * sqrt(n))){
  ## numeric_vector numeric_vector integer -> numeric
  ## produces gamma coefficient, given two vectors and number of upper order statistics
  ## ASSUME: v1 and v2 have the same length

  n <- NROW(v1)
  u  <- 1 - k/n   # probability in the tail

  r1 <- rank(v1, ties.method = "first") # ranks of v1
  r2 <- rank(v2, ties.method = "first") # ranks of v2

  1/(k * (n + 1)) * sum(r2[r1 > n - k + 1/2]) # produces causal gamma

}

gamma_theo <- function(node1, node2, adj.mat, alpha){
  ## integer integer numeric_matrix numeric_matrix numeric -> numeric
  ## compute \Gamma(node1, node2) given an adjacency matrix and a tail index

  ## helpers ####
  compute_betas <- function(paths, adj.mat){
    ## list numeric_matrix -> numeric_vector
    ## produce a vector with the weights for each path in paths

    betas <- sapply(paths, function(x){
      b <- 1
      for(i in 1:(length(x)-1)){
        row.ind <- x[i]
        col.ind <- x[i + 1]
        b <- b * adj.mat[row.ind, col.ind]
      }
      if(identical(b, numeric(0))){b <- 0}
      return(b)
    })
    if(length(betas) == 0){betas <- 0}

    return(betas)
  }
  ## function body ####
  # Compute ancestors
  rel <- findRelatives(adj.mat)
  ancestors <- rel$ancestors

  an_node1 <- which(ancestors[node1, ] == 1)
  an_node2 <- which(ancestors[node2, ] == 1)
  nan_node2 <- which(ancestors[node2, ] == 0)

  # Check if node1 is ancestor of node2
  if(node1 == node2){
    return(NA)
  }else if(node1 %in% an_node2){
    return(1)
  }

  # Compute parts of p_12
  A <- intersect(an_node1, nan_node2)
  B <- intersect(an_node1, an_node2)

  A_vars <- sapply(X = A[-which(node1 == A)], pathFromTo, dest = node1, adj.mat = adj.mat)
  A_vars <- lapply(rapply(A_vars, enquote, how = "unlist"), eval)
  A_betas <- compute_betas(A_vars, adj.mat)

  B_vars <- sapply(X = B, pathFromTo, dest = node1, adj.mat = adj.mat)
  B_vars <- lapply(rapply(B_vars, enquote, how = "unlist"), eval)
  B_betas <- compute_betas(B_vars, adj.mat)

  # Compute p_12
  p_12 <- (1 + sum(A_betas^alpha)) / (1 + sum(A_betas^alpha) + sum(B_betas^alpha))

  # Return \Gamma(node1, node2)
  return(1 - p_12/2)
} # !!! adapt to the case where the noise is scaled by a constant \neq 1

## Graph theory functions ####
checkCausOrder <- function(caus.order, adj.mat){
  ## numeric_vector numeric_matrix -> boolean
  ## produces true if given causal order agrees with adj.mat, false otherwise

  # Check causal order
  if(all(is.na(caus.order))){
    return(NA)
  }

  # Make sure adj.mat has only 0-1 entries
  adj.mat <- (adj.mat != 0) * 1
  p <- NROW(adj.mat)

  for (i in 1:p){
    current.node <- caus.order[i]
    indegree <- sum(adj.mat[, current.node])

    if (indegree > 0) {
      return(FALSE)
    }

    adj.mat[current.node, ] <- 0
  }

  return(TRUE)
}

findRelatives <- function(adj.mat){
  ## numeric_matrix -> list
  ## produces list given an adjacency matrix. The list is made of:
  ## - ancestors: matrix with all ancestors for each node
  ## - descendants: matrix with all descendants for each node
  ## - parents: matrix with all parents for each node
  ## - children: matrix with all children for each node
  ## ASSUME:
  ##   each node is ancestor and descendant of itself

  n <- NROW(adj.mat)
  adj.mat <- (adj.mat != 0) * 1

  A0 <- adj.mat
  B0 <- adj.mat

  for (i in 1:(n-1)){
    B0 <- B0 %*% adj.mat
    A0 <- A0 + B0
  }

  A0 <- (A0 > 0) * 1

  ancestors <- t(A0) + diag(n)
  descendants <- A0  + diag(n)
  parents <- t(adj.mat)
  children <- adj.mat

  ll <- list(ancestors = ancestors,
             descendants = descendants,
             parents = parents,
             children = children)

  return(ll)
}

pathFrom <- function(node, adj.mat){
  ## integer numeric_matrix <- list
  ## produce a list with all directed paths leaving from the given node


  # find children
  children <- which(adj.mat[node, ] != 0)

  # find paths from node
  if (identical(children,integer(0))){

    dir.path <- list(c(node))
    return(dir.path)

  } else {

    dir.path <- c()
    for (ch in children){

      returned_list <- pathFrom(ch, adj.mat)
      dir.path.ch <- lapply(returned_list, function(x){c(node, x)})
      dir.path <- c(dir.path, dir.path.ch)

    }
    return(dir.path)
  }
}

pathFromTo <- function(source, dest, adj.mat){
  ## integer numeric_matrix <- list
  ## produce a list with directed paths from source to destination node

  paths <- pathFrom(source, adj.mat)

  dir.paths <- c()
  for (i in 1:length(paths)){

    vec <- paths[[i]]

    if(dest %in% vec){

      vec <- vec[1:which(vec == dest)]
      if (!Position(function(x) identical(x, vec), dir.paths, nomatch = 0) > 0){
        dir.paths <- c(dir.paths, list(vec))}
    }
  }

  return(dir.paths)

}

ancestral_distance <-function(adj.mat, order){
  ## numeric_matrix numeric_vector -> numeric
  ## produce the number of inversions to transform order in a valid permutation
  ## Let pi(i) be the permutation of a node i = 1, ..., p.
  ## An inversion is a pair of nodes (i, j) where i is ancestor of j
  ## and pi(i) > pi(j).

  p <- NROW(adj.mat)
  max_inversion <- p * (p - 1) / 2
  ancestors <- findRelatives(adj.mat)$ancestors
  ancestors_ord <- ancestors[order, order]
  return(sum(ancestors_ord[upper.tri(ancestors_ord)])/max_inversion)
}

structInterv_distance <- function(adj.mat, est_adj.mat){
  ## numeric_matrix numeric_matrix -> numeric
  ## produces the structural intervention distance between the
  ## estimated DAG (or CPDAG) and the real DAG

  if(length(est_adj.mat) == 1){
    if(is.na(est_adj.mat)){
      return(NA)
    }
  }

  p <- NROW(adj.mat)
  s <- structIntervDist(adj.mat, est_adj.mat)
  return(s$sidLowerBound/(p * (p - 1)))
}

perm2adjmat <- function(perm){
  ## numeric_vector -> numeric_matrix
  ## produces a full connected adjacency matrix based on the given permutation

  p <- length(perm)
  adj.mat <- upper.tri(x = matrix(0, nrow = p, ncol = p)) * 1
  inv_perm <- order(perm)
  adj.mat <- adj.mat[inv_perm, inv_perm]
}

## Data analysis ####
rpareto <- function(n, alpha, scale=1){
  u <- runif(n)
  # u = 1 - (scale/x)^alpha
  # (scale/x)^alpha = 1 - u
  # scale/x = (1 - u)^(1/alpha)
  # scale = x * (1 - u)^(1/alpha)
  # x = scale * (1 - u)^(-1/alpha)
  x <- scale * (1 - u)^(-1/alpha)
  return(x)
}

size_freq_plot <- function(data){
  # numeric_matrix or dataframe -> plots
  # produces size-frequency plots for the given dataset

  # Silence warnings
  oldw <- getOption("warn")
  options(warn = -1)

  # Plot data
  for(i in 1:NCOL(data)) {
    x <- data[, i]
    x <- x[x > 0]
    r <- rank(-x)
    x <- x/max(x)
    r <- r/max(r)

    plot(x, r, log = "xy", pch = 20)
    legend("topright", paste("Variable", names(data)[i]), pch = 20)

    # Ask the user for input
    readline(paste("Press Enter to advance. Plot number", i))

    # Unsilence warnings
    options(warn = oldw)

  }
}
