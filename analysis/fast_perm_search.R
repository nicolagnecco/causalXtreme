# Copyright (c) 2018-2018  Jonas Peters [jonas.peters@math.ku.dk]
# All rights reserved.  See the file COPYING for license terms.



fast_perm_search <- function(A, silent = FALSE, mode="sum"){
  d <- dim(A)[2]
  var.names <- 1:d
  if(!silent){
    show("###############")
    show("current score")
    show(var.names)
    show(A)
}
  if(mode=="sum"){
      summ <- apply(A,1,sum, na.rm=TRUE)
      avail <- 1:d
      current.order <- numeric(0)
      for (k in 1:d){
          add <-  avail[which.max(summ[avail])]
          current.order <- c(current.order, add)
          summ <- summ + A[ add,]
          avail <- (1:d)[ -current.order]
      }
      score <- sum(A[current.order,current.order][upper.tri(A)])
  }else if(mode=="maxmin"){
      Aorig <- A
      diag(A) <- NA
      ## initialize with the first pair
      current.order <- add <- which.max(apply(A, 1 , min, na.rm=TRUE))
      for (k in 2:d){
          A[, add] <- NA
          avail <- (1:d)[-current.order]
          
          add <- if(k<d) avail[which.max(apply(A, 1 , min, na.rm=TRUE)[avail])] else avail
          current.order <- c(current.order, add)
      }
  
      ##  score
      Aorig <- Aorig[ current.order,current.order]
      score <- sum(Aorig[upper.tri(Aorig)])
  }
  
  return(list(score = score, order = current.order))
}


