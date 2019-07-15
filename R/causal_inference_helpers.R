#' Random search
#'
#' Produces a random DAG with random sparsity with as many variables
#' as the number of columns in the dataset \code{dat}.
#'
#' @inheritParams causal_tail_matrix
#' @return Square binary matrix. A random DAG.
random_search <- function(dat){

  # number of variables
  p <- NCOL(dat)

  # generate random DAG
  dag <- random_dag(p = p, prob_connect = runif(1))

  # return DAG
  return(dag)
}


#' Greedy ancestral search
#'
#' Runs greedy ancestral search permutation search given the dataset
#' \code{dat}.
#'
#' @inheritParams causal_tail_matrix
#' @return  Numeric vector. The causal order estimated from the data.
greedy_ancestral_search <- function(dat, k = floor(2 * n ^ 0.4),
                                    both_tails = TRUE){
  # set up variables
  n <- NROW(dat)
  d <- NCOL(dat)

  # compute causal tail matrix
  causal_mat <- causal_tail_matrix(dat, k, both_tails)

  # run greedy ancestral search
  current_order <- add <- which.min(apply(causal_mat, 2, max, na.rm = TRUE))
  for (k in 2:d){
    causal_mat[add, ] <- NA
    avail <- (1:d)[-current_order]

    add <- if (k < d){
      avail[which.min(apply(causal_mat, 2, max, na.rm = TRUE)[avail])]
    } else {
      avail
    }
    current_order <- c(current_order, add)
  }

  order <- current_order

  # return causal order
  return(order)
}


#' Lingam search
#'
#' Runs Lingam given a dataset \code{dat}.
#'
#' @inheritParams causal_tail_matrix
#' @param contrast_fun Character. The functional form of the contrast
#' function used in the Fast-ICA step. It is one of \code{"logcosh"}
#' (the default choice) and \code{"exp"}.
#' For further details see the paper from
#' Hyvarinen, A., \url{https://ieeexplore.ieee.org/abstract/document/761722/}.
#' @return Square binary matrix (or \code{NA} in case of error).
#' The DAG estimated from the data.
lingam_search <- function(dat, contrast_fun = c("logcosh", "exp")[1]){

  out <- tryCatch({
    t.k <- estLiNGAM(dat, only.perm = T, fun = contrast_fun)$k
    lingam_output <- prune(t(dat), t.k)
    Bpruned <- lingam_output$Bpruned
    dag <- (t(Bpruned) != 0) * 1
    return(dag)
  },
  error = function(e){
    dag <- NA
    return(dag)
  })

  return(out)
}


#' PC search
#'
#' Runs PC algorithm given a dataset \code{dat}.
#'
#' The function \code{pcalg:pc()} is called with the following
#' arguments:
#' \itemize{
#' \item \code{suffStat = list(C = cor(dat), n = NROW(dat))},
#' \item \code{indepTest = gaussCItest},
#' \item \code{u2pd = "retry"} --- this ensures that the produced CPDAG
#' is extendible to a DAG,
#' \item \code{skel.method = "stable"}.
#' }
#'
#' @inheritParams causal_tail_matrix
#' @param alpha Numeric --- between 0 and 1. The significance level for the
#' individual conditional independence tests.
#' @return Square binary matrix (or \code{NA} in case of error).
#' The CPDAG estimated from the data.
pc_search <- function(dat, alpha){

  n <- NROW(dat)
  p <- NCOL(dat)

  out <- tryCatch({
    suff_stat <- list(C = cor(dat), n = n)
    pc.fit <- pcalg::pc(suffStat = suff_stat,
                        indepTest = pcalg::gaussCItest,
                        p = p, alpha = alpha, u2pd = "retry",
                        skel.method = "stable")
    cpdag <- as(pc.fit@graph, "matrix")
    return(cpdag)
  },
  error = function(e){
    cpdag <- NA
    return(cpdag)
  })

  return(out)
}


#' PC (rank) search
#'
#' Runs PC (rank) algorithm given a dataset \code{dat}.
#'
#' The PC (rank) algorithm has been proposed by Harris, N., and
#' Drton, M.,
#' \url{http://jmlr.csail.mit.edu/papers/volume14/harris13a/harris13a.pdf}
#'
#' The function \code{pcalg:pc()} is called with the following
#' arguments:
#' \itemize{
#' \item \code{suffStat = list(C = 2 * sin(cor(dat, method = "spearman")
#' * pi/6), n = nrow(dat))},
#' \item \code{indepTest = gaussCItest},
#' \item \code{u2pd = "retry"} --- this ensures that the produced CPDAG
#' is extendible to a DAG,
#' \item \code{skel.method = "stable"}.
#' }
#'
#' @inheritParams causal_tail_matrix
#' @param alpha Numeric --- between 0 and 1. The significance level for the
#' individual conditional independence tests.
#' @return Square binary matrix (or \code{NA} in case of error).
#' The CPDAG estimated from the data.
pc_rank_search <- function(dat, alpha){

  n <- NROW(dat)
  p <- NCOL(dat)

  out <- tryCatch({
    suff_stat <-  list(C = 2 * sin(cor(dat, method = "spearman") * pi/6),
                       n = nrow(dat))
    pc.fit <- pcalg::pc(suffStat = suff_stat,
                        indepTest = pcalg::gaussCItest,
                        p = p, alpha = alpha, u2pd = "retry",
                        skel.method = "stable")
    cpdag <- as(pc.fit@graph, "matrix")
    return(cpdag)
  },
  error = function(e){
    cpdag <- NA
    return(cpdag)
  })

  return(out)
}
