library(gtools)

full_perm_search <- function(A, silent = FALSE){
  d <- dim(A)[2]
  allPerms <- permutations(d,d)
  scores <- rep(NA, dim(allPerms)[1])
  for(i in 1:length(scores)){
    Atmp <- A[allPerms[i,], allPerms[i,]]
    scores[i] <- sum(Atmp[upper.tri(Atmp)])
  }
  ind.max <- which.max(scores)
  #show(sort(scores, decreasing = TRUE))
  #show(scores[ind.max])
  return(list(score = scores[ind.max], order = allPerms[ind.max,]))
}
