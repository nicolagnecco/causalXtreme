direct_lingam <- function(X){
  # Step 1
  # prepare data
  X_orig <- X
  p <- NCOL(X)
  n <- NROW(X)

  # Center variables
  X_1 <- drop(scale(X, scale = FALSE))
  attr(X_1, "scaled:center") <- NULL

  # prepare M matrix
  M <- matrix(-1, ncol = p, nrow = p)
  diag(M) <- 0

  # prepare vector K and U - K
  K <- numeric(p)
  U_K <- 1:p

  # Step 2
  for (m in 1:(p - 1)){
    # find exogenous by using M
    exogenous <- which(colSums(t(M) == 0) == p - m + 1)

    if (identical(exogenous, integer(0))){
      endogenous <- which(colSums(t(M) == 1) > 0)
      candidates <- setdiff(U_K, endogenous)
    } else{
      candidates <- exogenous
    }


    # (a)
    # compute R, i.e., matrix containing residuals
    R <- computeR(X, candidates, U_K, M)


  }


  return(1:p) # estimated causal ordering
}

computeR <- function(X, candidates, U_K, M){

  p <- NCOL(X)
  n <- NROW(X)
  R <- array(dim = c(p, n, p))

  return(R)

}
X <- matrix(rt(100*4, 2), ncol = 4)
X_1
