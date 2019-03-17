random_perm_search <- function(A){
  # Copyright (c) 2018-2018  Jonas Peters [jonas.peters@math.ku.dk]
  # All rights reserved.  See the file COPYING for license terms.
  
  p <- NROW(A) # number of variables
  order <-  sample(p, p, replace = FALSE)
  A <- A[order, order]
  score <- sum(A[upper.tri(A)])
  
  return(list(score = score, order = order))
}