# Copyright (c) 2018-2018  Jonas Peters [jonas.peters@math.ku.dk]
# All rights reserved.  See the file COPYING for license terms.



greedy_perm_search <- function(A, silent = FALSE){
  d <- dim(A)[2]
  var.names <- 1:d
  if(!silent){
    show("###############")
    show("current score")
    show(var.names)
    show(A)
  }
  
  # initialize with the first pair
  current.order <- as.vector(arrayInd(which.max(A), .dim = dim(A)))
  if(!silent){
    show("current order")
    show(current.order)
  }
  var.names <- var.names[-current.order]
  score <- max(A, na.rm=TRUE)
  current.score <- matrix(NA, d-2, 3)
  current.score[,1] <- A[-current.order,current.order[1]] + A[-current.order,current.order[2]]
  current.score[,2] <- t(A[current.order[1], -current.order]) + A[-current.order,current.order[2]]
  current.score[,3] <- t(A[current.order[1], -current.order]) + t(A[current.order[2], -current.order])
  if(!silent){
    show("current order")
    show(current.order)
  }

  for(k in 3:(d-1)){ # k-1 is the number of variables already included 
    if(!silent){
      show("###############")
      show("current score")
      show(var.names)
      show(current.score)
    }
    new.max <- arrayInd(which.max(current.score), .dim = dim(current.score))
    new.var.ind <- new.max[1]
    new.var <- var.names[new.var.ind]
    new.pos <- new.max[2]
    
    # include the new variable into the order
    tmp <- current.order
    current.order <- rep(NA, k)
    current.order[new.pos] <- new.var
    current.order[-new.pos] <- tmp
    if(!silent){
      show("current order")
      show(current.order)
    }
    
    # update var.names
    var.names <- var.names[-new.var.ind]
    
    # update score
    score <- score + current.score[new.max]
    
    # update score matrix
    tmp <- current.score
    current.score <- matrix(NA, d - length(current.order), k+1)
    current.score[,1:new.pos] <- tmp[-new.var.ind,1:new.pos] + A[-current.order,new.var]
    current.score[,(1+new.pos):(k+1)] <- tmp[-new.var.ind,new.pos:k] + matrix(A[new.var,-current.order], d - length(current.order) , k-new.pos + 1) 
  }
  
  if(!silent){
    show("###############")
    show("current score")
    show(var.names)
    show(current.score)
  }
  k <- d
  new.max <- arrayInd(which.max(current.score), .dim = dim(current.score))
  new.var.ind <- new.max[1]
  new.var <- var.names[new.var.ind]
  new.pos <- new.max[2]
  
  # include the new variable into the order
  tmp <- current.order
  current.order <- rep(NA, k)
  current.order[new.pos] <- new.var
  current.order[-new.pos] <- tmp
  if(!silent){
    show("current order")
    show(current.order)
  }
  
  # update score
  score <- score + current.score[new.max]
  
  return(list(score = score, order = current.order))
}


