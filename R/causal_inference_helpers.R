#' Random search
#'
#' Produces a random DAG with random sparsity with as many variables
#' as the number of columns in the dataset \code{dat}.
#'
#' @inheritParams causal_tail_matrix
#' @return Square binary matrix. A random DAG.
#' @export
random_search <- function(dat){

  # number of variables
  p <- NCOL(dat)

  # generate random DAG
  dag <- random_dag(p = p, prob_connect = stats::runif(1))

  # return DAG
  return(dag)
}


#' Extremal Ancestral SEarch
#'
#' Runs extremal ancestral search (EASE) algorithm on the given
#' dataset \code{dat}.
#'
#' @inheritParams causal_tail_matrix
#' @return  Numeric vector. The causal order estimated from the data.
#' @export
ease <- function(dat, k = floor(n ^ 0.4),
                                    both_tails = TRUE){
  # set up variables
  n <- NROW(dat)
  d <- NCOL(dat)

  # compute causal tail matrix
  causal_mat <- causal_tail_matrix(dat, k, both_tails)

  # run Extremal Ancestral SEarch
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
#' Runs Lingam given a dataset \code{dat}. Returns the DAG estimated from the
#' data using the LiNGAM-ICA algorithm.
#'
#' This function is a wrapper around \code{estLiNGAM} from the package
#' \code{pcalg} (see \url{https://CRAN.R-project.org/package=pcalg }).
#' The function \code{estLiNGAM} is slightly modified
#' to allow for different contrast functions in the fast-ICA step of
#' LiNGAM. To modify \code{estLiNGAM}, we included in this package
#' several internal functions from\code{pcalg}.
#' All the credits go to the authors of the \code{pcalg} package:
#'
#' Markus Kalisch, Alain Hauser , Martin Maechler, Diego Colombo,
#' Doris Entner, Patrik Hoyer, Antti Hyttinen, Jonas Peters,
#' Nicoletta Andri, Emilija Perkovic, Preetam Nandy, Philipp Ruetimann,
#' Daniel Stekhoven, Manuel Schuerch, Marco Eigenmann.
#'
#' @inheritParams causal_tail_matrix
#' @param contrast_fun Character. The functional form of the contrast
#' function used in the Fast-ICA step. It is one of \code{"logcosh"}
#' (the default choice) and \code{"exp"}.
#' For further details see the paper from
#' Hyvarinen, A., \url{https://ieeexplore.ieee.org/abstract/document/761722/}.
#' @return Square binary matrix (or \code{NA} in case of error).
#' The DAG estimated from the data.
#' @export
lingam_search <- function(dat, contrast_fun = c("logcosh", "exp")){

  contrast_fun <- match.arg(contrast_fun)

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


#' Order-Lingam search
#'
#' Runs Order-Lingam given a dataset \code{dat}. Returns the causal order
#' estimated from the data using the first step of
#' the LiNGAM-ICA algorithm.
#'
#' This function is a wrapper around \code{estLiNGAM} from the package
#' \code{pcalg} (see \url{https://CRAN.R-project.org/package=pcalg }).
#' The function \code{estLiNGAM} is slightly modified
#' to allow for different contrast functions in the fast-ICA step of
#' LiNGAM. To modify \code{estLiNGAM}, we included in this package
#' several internal functions from\code{pcalg}.
#' All the credits go to the authors of the \code{pcalg} package:
#'
#' Markus Kalisch, Alain Hauser , Martin Maechler, Diego Colombo,
#' Doris Entner, Patrik Hoyer, Antti Hyttinen, Jonas Peters,
#' Nicoletta Andri, Emilija Perkovic, Preetam Nandy, Philipp Ruetimann,
#' Daniel Stekhoven, Manuel Schuerch, Marco Eigenmann.
#'
#' @inheritParams causal_tail_matrix
#' @param contrast_fun Character. The functional form of the contrast
#' function used in the Fast-ICA step. It is one of \code{"logcosh"}
#' (the default choice) and \code{"exp"}.
#' For further details see the paper from
#' Hyvarinen, A., \url{https://ieeexplore.ieee.org/abstract/document/761722/}.
#' @return Numeric vector (or \code{NA} in case of error).
#' The causal order estimated from the data.
#' @export
order_lingam_search <- function(dat, contrast_fun = c("logcosh", "exp")){

  contrast_fun <- match.arg(contrast_fun)

  out <- tryCatch({
    order <- estLiNGAM(dat, only.perm = T, fun = contrast_fun)$k
    return(order)
  },
  error = function(e){
    order <- NA
    return(order)
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
#' individual conditional independence tests. By default it is set to 0.05.
#' @return Square binary matrix (or \code{NA} in case of error).
#' The CPDAG estimated from the data.
#' @export
pc_search <- function(dat, alpha = 5e-2){

  n <- NROW(dat)
  p <- NCOL(dat)

  out <- tryCatch({
    suff_stat <- list(C = stats::cor(dat), n = n)
    pc.fit <- pcalg::pc(suffStat = suff_stat,
                        indepTest = pcalg::gaussCItest,
                        p = p, alpha = alpha, u2pd = "retry",
                        skel.method = "stable")
    cpdag <- unname(methods::as(pc.fit@graph, "matrix"))
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
#' @inheritParams pc_search
#' @return Square binary matrix (or \code{NA} in case of error).
#' The CPDAG estimated from the data.
#' @export
pc_rank_search <- function(dat, alpha = 5e-2){

  n <- NROW(dat)
  p <- NCOL(dat)

  out <- tryCatch({
    suff_stat <-  list(C = 2 * sin(stats::cor(dat, method = "spearman") *
                                     pi / 6),
                       n = n)
    pc.fit <- pcalg::pc(suffStat = suff_stat,
                        indepTest = pcalg::gaussCItest,
                        p = p, alpha = alpha, u2pd = "retry",
                        skel.method = "stable")
    cpdag <- unname(methods::as(pc.fit@graph, "matrix"))
    return(cpdag)
  },
  error = function(e){
    cpdag <- NA
    return(cpdag)
  })

  return(out)
}
