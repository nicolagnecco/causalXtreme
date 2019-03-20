#' Greedy permutation search
#'
#' Runs greedy permutation search algorithm given the matrix \code{delta}.
#' Copyright (c) 2013 Jonas Peters \email{peters@@stat.math.ethz.ch}.
#' All rights reserved.
#'
#' @param delta Numeric matrix --- between -1 and 1. The \code{delta}
#' matrix defined as \code{gamma - t(gamma)}, where \code{gamma}
#' is the gamma coefficient matrix. The number of rows (and columns)
#' of \code{delta} must be greater than 3.
#' @param silent Logical. Should the function communicate during the
#' execution?
#' @return List. A list made of:
#' \itemize{
#' \item \code{order} Numeric vector. The causal order estimated from \code{delta},
#' \item \code{score} Numeric. The score associated to \code{order}.
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
