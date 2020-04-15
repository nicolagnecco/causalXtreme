
# Main ####
direct_lingam <- function(X){
  # perform direct lingam on dataset X

  # transpose matrix
  X <- t(X)


  # Step 1
  # prepare data
  X_orig <- X
  p <- NROW(X)
  n <- NCOL(X)

  # center variables
  X <- centerRows(X)

  # prepare matrix M
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

    # skip exogenous finding if it is found
    if (length(candidates) == 1){
      index <- candidates
    } else {
      # find exogenous
      index <- findindex(X, R, candidates, U_K)
    }

    # (b)
    K[m] <- index
    U_K <- U_K[U_K != index]
    M[ , index] <- NA
    M[index, ] <- NA

    # (c)
    X <- R[ , , index]

  }


  # Step 3
  # update last entry of causal ordering
  K[p] <- U_K

  # return estimated causal ordering
  return(K)
}


# Helpers ####
centerRows <- function(X){
  n <- NCOL(X)
  X - matrix(rep(rowMeans(X), n), ncol = n)
}

computeR <- function(X, candidates, U_K, M){
  # compute residual matrix when regressing

  p <- NROW(X)
  n <- NCOL(X)
  R <- array(0, dim = c(p, n, p))
  Cov <- cov(t(X))

  for (j in candidates){
    for (i in setdiff(U_K, j)){
     # skip residual calculations by looking at M
      if (M[i, j] == 0) {
        R[i, , j] <- X[i, ]
      } else{
        R[i, , j] <- X[i, ] - Cov[i, j] / Cov[j, j] * X[j, ]
      }
    }
  }
  return(R)

}

findindex <- function(X, R, candidates, U_K){
  # find root node

  p <- NROW(X)
  T_vec <- rep(NaN, p)

  minT <- -1

  for (j in candidates){
    if (minT == -1 & !is.nan(minT)){

      T_vec[j] <- 0

      for (i in setdiff(U_K, j)){
        J <- mutual_info(rbind(R[i, , j], X[j, ])) # or mutual_info(rbind(X[i, ], X[j, ]))
        T_vec[j] <- T_vec[j] + J
      }

      minT <- T_vec[j]

    } else {

      T_vec[j] <- 0

      for (i in setdiff(U_K, j)){
        J <- mutual_info(rbind(R[i, , j], X[j, ])) # or mutual_info(rbind(X[i, ], X[j, ]))
        T_vec[j] <- T_vec[j] + J

        if (T_vec[j] > minT & !is.nan(minT)){

          T_vec[j] <- Inf
          break
        }
      }

      minT <- min(c(T_vec[j], minT))
    }
  }

  if (all(is.nan(T_vec))){
    index <- 1
  } else {
    index <- which.min(T_vec)
  }
  return(index)
}

mutual_info <- function(x){
  # Try kernel ica
  m <- NROW(x)
  N <- NCOL(x)

  contrast='kgv';

  if (N < 1000){
    sigma=1;
    kappa=2e-2;

  } else{
    sigma = 1/2;
    kappa = 2e-3;
  }

  kernel='gaussian';

  mc=m;
  kparam <- list()
  kparam$kappas=kappa*rep(1,mc);
  kparam$etas=kappa*1e-2*rep(1,mc);
  kparam$neigs=N*rep(1,mc);
  kparam$nchols=N*rep(1,mc);
  kparam$kernel=kernel;
  kparam$sigmas=sigma*rep(1,mc);

  J = contrast_ica(contrast,x,kparam)

  return(J)
}



