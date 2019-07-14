#' Greedy permutation search
#'
#' Runs greedy permutation search given the matrix \code{delta}.
#' Copyright (c) 2013 Jonas Peters \email{peters@@math.ku.dk}.
#' All rights reserved.
#'
#' @param delta Numeric matrix --- entries between -1 and 1. The \code{delta}
#' matrix is defined as \code{gamma - t(gamma)}, where \code{gamma}
#' is the gamma coefficient matrix. The number of rows (and columns)
#' of \code{delta} must be greater than 3.
#' @param silent Logical. Should the function communicate during the
#' execution?
#' @return List. A list made of:
#' \itemize{
#' \item \code{order} Numeric vector. The causal order estimated
#' from \code{delta},
#' \item \code{score} Numeric. The maximized score associated to \code{order}.
#' It is defined as the sum of the upper triangular entries of \code{delta},
#' after its rows and columns are sorted according to \code{order}.
#' }
greedy_perm_search <- function(delta, silent = FALSE){

  A <- delta
  d <- dim(A)[2]

  if (d < 4){
    stop("The number of rows (and columns) of delta must be greater than 3.")
  }

  var_names <- 1:d
  if (!silent){
    show("###############")
    show("current score")
    show(var_names)
    show(A)
  }

  # initialize with the first pair
  current_order <- as.vector(arrayInd(which.max(A), .dim = dim(A)))
  if (!silent){
    show("current order")
    show(current_order)
  }
  var_names <- var_names[-current_order]
  score <- max(A, na.rm = TRUE)
  current_score <- matrix(NA, d - 2, 3)
  current_score[, 1] <- A[-current_order, current_order[1]] +
    A[-current_order, current_order[2]]
  current_score[, 2] <- t(A[current_order[1], -current_order]) +
    A[-current_order, current_order[2]]
  current_score[, 3] <- t(A[current_order[1], -current_order]) +
    t(A[current_order[2], -current_order])
  if (!silent){
    show("current order")
    show(current_order)
  }

  for (k in 3:(d - 1)){
    # k-1 is the number of variables already included
    if (!silent){
      show("###############")
      show("current score")
      show(var_names)
      show(current_score)
    }
    new_max <- arrayInd(which.max(current_score), .dim = dim(current_score))
    new_var_ind <- new_max[1]
    new_var <- var_names[new_var_ind]
    new_pos <- new_max[2]

    # include the new variable into the order
    tmp <- current_order
    current_order <- rep(NA, k)
    current_order[new_pos] <- new_var
    current_order[-new_pos] <- tmp
    if (!silent){
      show("current order")
      show(current_order)
    }

    # update var_names
    var_names <- var_names[-new_var_ind]

    # update score
    score <- score + current_score[new_max]

    # update score matrix
    tmp <- current_score
    current_score <- matrix(NA, d - length(current_order), k + 1)
    current_score[, 1:new_pos] <- tmp[-new_var_ind, 1:new_pos] +
      A[-current_order, new_var]
    current_score[, (1 + new_pos):(k + 1)] <- tmp[-new_var_ind, new_pos:k] +
      matrix(A[new_var, -current_order],
             d - length(current_order),
             k - new_pos + 1)
  }

  if (!silent){
    show("###############")
    show("current score")
    show(var_names)
    show(current_score)
  }
  k <- d
  new_max <- arrayInd(which.max(current_score), .dim = dim(current_score))
  new_var_ind <- new_max[1]
  new_var <- var_names[new_var_ind]
  new_pos <- new_max[2]

  # include the new variable into the order
  tmp <- current_order
  current_order <- rep(NA, k)
  current_order[new_pos] <- new_var
  current_order[-new_pos] <- tmp
  if (!silent){
    show("current order")
    show(current_order)
  }

  # update score
  score <- score + current_score[new_max]

  return(list(order = current_order, score = score))
}


#' Fast permutation search
#'
#' Runs fast permutation search given the matrix \code{delta}.
#' Copyright (c) 2013 Jonas Peters \email{peters@@math.ku.dk}.
#' All rights reserved.
#'
#' @inheritParams greedy_perm_search
#' @param delta Numeric matrix --- entries between -1 and 1. The \code{delta}
#' matrix is defined as \code{gamma - t(gamma)}, where \code{gamma}
#' is the gamma coefficient matrix.
#' @param mode String. Is one of:
#' \itemize{
#' \item \code{"sum"} (default). At each step of the algorithm, choose
#' the node \eqn{i}, if the \emph{sum} of the entries of the \eqn{i}-th
#' row in \code{delta} is maximum across all the remaining nodes.
#' \item \code{"maxmin"}. At each step of the algorithm, choose the node
#' \eqn{i}, if the the \eqn{i}-th row of \code{delta} is the one which has
#' the maximum minimum (\emph{maxmin}) entry across all the remaining nodes.
#' }
#' @inherit greedy_perm_search return
fast_perm_search <- function(delta, mode = "sum", silent = FALSE){

  A <- delta
  d <- dim(A)[2]
  var_names <- 1:d
  if (!silent){
    show("###############")
    show("current score")
    show(var_names)
    show(A)
  }
  if (mode == "sum"){
    summ <- apply(A, 1, sum, na.rm = TRUE)
    avail <- 1:d
    current_order <- numeric(0)
    for (k in 1:d){
      add <-  avail[which.max(summ[avail])]
      current_order <- c(current_order, add)
      summ <- summ + A[add, ]
      avail <- (1:d)[ -current_order]
    }
    score <- sum(A[current_order, current_order][upper.tri(A)])
  } else if (mode == "maxmin"){
    Aorig <- A
    diag(A) <- NA
    ## initialize with the first pair
    current_order <- add <- which.max(apply(A, 1, min, na.rm = TRUE))
    for (k in 2:d){
      A[, add] <- NA
      avail <- (1:d)[-current_order]

      add <- if (k < d){
        avail[which.max(apply(A, 1, min, na.rm = TRUE)[avail])]
      } else {
        avail
      }
      current_order <- c(current_order, add)
    }

    ##  score
    Aorig <- Aorig[ current_order, current_order]
    score <- sum(Aorig[upper.tri(Aorig)])

  } else {
    stop("Wrong mode. Please enter one of the following:
              'sum', 'maxmin'.")
  }

  return(list(order = current_order, score = score))
}


#' Full permutation search
#'
#' Runs full permutation search given the matrix \code{delta}.
#' Copyright (c) 2013 Jonas Peters \email{peters@@math.ku.dk}.
#' All rights reserved.
#'
#' @inheritParams greedy_perm_search
#' @param delta Numeric matrix --- entries between -1 and 1. The \code{delta}
#' matrix is defined as \code{gamma - t(gamma)}, where \code{gamma}
#' is the gamma coefficient matrix. If the number of rows (and columns)
#' of \code{delta} is greater than 10, the function returns an error.
#' @inherit greedy_perm_search return
full_perm_search <- function(delta, silent = FALSE){

  A <- delta
  d <- dim(A)[2]

  if (d > 10){
    stop("The number of rows (and columns) of delta cannot be greater
         than 10.")
  }

  all_perms <- gtools::permutations(d, d)
  scores <- rep(NA, dim(all_perms)[1])
  for (i in 1:length(scores)){
    Atmp <- A[all_perms[i, ], all_perms[i, ]]
    scores[i] <- sum(Atmp[upper.tri(Atmp)])
  }
  ind_max <- which.max(scores)

  return(list(order = all_perms[ind_max, ], score = scores[ind_max]))
}


#' Random permutation search
#'
#' Runs random permutation search given the matrix \code{delta}.
#' Copyright (c) 2013 Jonas Peters \email{peters@@math.ku.dk}.
#' All rights reserved.
#'
#' @inheritParams compute_caus_order
#' @return Numeric vector. The causal order estimated from \code{g}.
random_perm_search <- function(g){

  A <- g
  p <- NROW(A) # number of variables
  order <-  sample(p, p, replace = FALSE)

  return(order)
}


#' Minimax permutation search
#'
#' Runs minimax permutation search given the matrix \code{gamma}.
#'
#' @param gamma Numeric matrix --- entries between 0 and 1.
#' The \code{gamma} coefficient matrix.
#' @return  Numeric vector. The causal order estimated from \code{gamma}.
minimax_search <- function(gamma){

  # set up variables
  d <- dim(gamma)[2]
  diag(gamma) <- NA

  ## run minimax
  current_order <- add <- which.min(apply(gamma, 2, max, na.rm = TRUE))
  for (k in 2:d){
    gamma[add, ] <- NA
    avail <- (1:d)[-current_order]

    add <- if (k < d){
      avail[which.min(apply(gamma, 2, max, na.rm = TRUE)[avail])]
    } else {
      avail
    }
    current_order <- c(current_order, add)
  }

  order <- current_order

  return(order)
}


#' Oracle permutation search
#'
#' Runs oracle permutation search given the true DAG \code{g}.
#'
#' @param g Numeric matrix. The weighted adjacency matrix of a DAG.
#' @return Numeric vector. The causal order estimated from \code{g}.
oracle_search <- function(g){

  g <- (g != 0) * 1
  order <- compute_caus_order(g)
  return(order)

}


#' Lingam search
#'
#' Runs Lingam given a dataset \code{mat}.
#'
#' @inheritParams causal_tail_matrix
#' @return Numeric matrix (or \code{NA} in case of error).
#' The estimated adjacency matrix of a DAG.
lingam_search <- function(mat){

  out <- tryCatch({
      lingam_output <- pcalg::lingam(X = mat)
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
#' Runs PC algorithm given a dataset \code{mat}.
#'
#' The function \code{pcalg:pc()} is called with the following
#' arguments:
#' \itemize{
#' \item \code{suffStat = list(C = cor(mat), n = NROW(mat))},
#' \item \code{indepTest = gaussCItest},
#' \item \code{alpha = 0.05},
#' \item \code{u2pd = "retry"} --- this makes sure that the produced CPDAG
#' is extendible to a DAG,
#' \item \code{skel.method = "stable"}.
#' }
#'
#' @inheritParams causal_tail_matrix
#' @return Numeric matrix (or \code{NA} in case of error).
#' The estimated adjacency matrix of a CPDAG.
pc_search <- function(mat){

  n <- NROW(mat)
  p <- NCOL(mat)

  out <- tryCatch({
    suff_stat <- list(C = cor(mat), n = n)
    pc.fit <- pcalg::pc(suffStat = suff_stat, indepTest = pcalg::gaussCItest,
                 p = p, alpha = 5e-2, u2pd = "retry", skel.method = "stable")
    cpdag <- as(pc.fit@graph, "matrix")
    return(cpdag)
    },
    error = function(e){
    cpdag <- NA
    return(cpdag)
  })

  return(out)
}
