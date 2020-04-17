direct_lingam_search_r <- function(dat){
  ## numeric_matrix -> causal_order
  # perform direct lingam on dataset X

  # transpose matrix
  dat <- t(dat)


  # Step 1
  # prepare data
  dat_orig <- dat
  p <- NROW(dat)
  n <- NCOL(dat)

  # center variables
  dat <- center_rows(dat)

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
    R <- computeR(dat, candidates, U_K, M)
    # R <- computeRc(dat, candidates - 1, U_K - 1, M)

    # skip exogenous finding if it is found
    if (length(candidates) == 1){
      index <- candidates
    } else {
      # find exogenous
      index <- findindex(dat, R, candidates, U_K)
      # index <- findindexc(dat, candidates - 1, U_K - 1)
    }

    # (b)
    K[m] <- index
    U_K <- U_K[U_K != index]
    M[ , index] <- NA
    M[index, ] <- NA

    # (c)
    dat <- R[ , , index]

  }


  # Step 3
  # update last entry of causal ordering
  K[p] <- U_K

  # return estimated causal ordering
  return(K)
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
        J <- pwlingc(rbind(X[i, ], X[j, ]))  #kernel_ica(rbind(R[i, , j], X[j, ])) # or
        T_vec[j] <- T_vec[j] + J
      }

      minT <- T_vec[j]

    } else {

      T_vec[j] <- 0

      for (i in setdiff(U_K, j)){
        J <- pwlingc(rbind(X[i, ], X[j, ])) #kernel_ica(rbind(R[i, , j], X[j, ])) # or
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

kernel_ica <- function(x){
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

pwling <- function(X){
  ## numeric_matrix -> numeric
  ## computes pwling measure given a matrix 2 x n

  # compute size of matrix
  p <- NROW(X)
  n <- NCOL(X)

  # check if p == 2
  if (p != 2){
    stop(paste("Matrix X must have exactly 2 rows."))
  }

  # standardize rows
  X <- t(scale(t(X)))

  # compute covariance matrix
  C <- cov(t(X))

  # General entropy-based method, for variables of any distribution
  i <- 2
  j <- 1
  res1 <- X[j, ] - C[j, i] * X[i, ]
  res2 <- X[i, ] - C[i, j] * X[j, ]
  LR <-  mentappr(X[j, ]) - mentappr(X[i, ]) - mentappr(res1) + mentappr(res2)

  # compute score
  J = min(0, LR)^2

  # return score
  return(J)
}

mentappr <- function(x){
  ## numeric_vector -> numeric
  ## computes maximum entropy approximation

  # standardize
  xstd <- sd(x)
  x <- scale(x)

  # Constants we need
  k1 <- 36 / (8 * sqrt(3) - 9)
  gamma <- 0.37457
  k2 <- 79.047
  gaussianEntropy <- log(2 * pi) / 2 + 0.5

  # This is negentropy
  negentropy <- k2 * (mean(log(cosh(x)))-gamma)^2+k1*mean(x*exp(-x^2/2))^2

  # This is entropy
  entropy = gaussianEntropy - negentropy + log(xstd);

  return(entropy)
}

contrast_ica <- function(contrast, x, kparam){
  # character matrix (pxN) list -> numeric

  # computes constrast kernel ICA
  N <- NCOL(x)
  m <- NROW(x)
  kappas <- kparam$kappas
  etas <- kparam$etas
  sizes <- numeric(m)

  Us <- list()
  Lambdas <- list()
  Drs <- list()

  for (i in 1:m){
    res <- chol_gaussc(x[i, , drop = FALSE] / kparam$sigmas[i], 1, N * etas[i])
    G <- res$G
    # file_nm <-  paste("analysis/G", i, ".csv", sep = "") # !!!
    # G <- as.matrix(read.csv(file_nm, header = FALSE)) # !!!
    Pvec <- res$Pvec

    k <- NCOL(G)
    a <- sort(Pvec)
    Pvec <- order(Pvec)

    G <- center_cols(G[Pvec, , drop = FALSE])

    # regularization (see paper for details)
    eigen_decomp <- eigen(crossprod(G))
    A <- eigen_decomp$vectors[, k:1, drop = FALSE]
    D <- eigen_decomp$values[k:1]


    # eigen_decomp <- symm_eigen(crossprod(G))
    # A <- eigen_decomp$vectors[, k:1, drop = FALSE]
    # D <- as.numeric(eigen_decomp$values)[k:1]

    indexes <- which(D >= N * etas[i] & is.double(D));  #removes small eigenvalues
    newinds <- sort(D[indexes])
    neworder <- order(D[indexes])
    neworder <- rev(neworder)
    neig <- length(indexes)

    indexes <- indexes[neworder[seq_len(neig)]]
    if (identical(indexes, integer(0))){
      indexes = 1
    }

    D <- D[indexes]
    V = G %*% (A[, indexes] %*% diag(sqrt(1 / D), nrow = length(D)))

    Us[[i]] <- V
    Lambdas[[i]] <- D
    Dr <- D

    for (j in seq_along(D)){
      Dr[j] <- D[j] / (N * kappas[i] + D[j])
    }

    Drs[[i]] <- Dr
    sizes[i] <- length(Drs[[i]])

  }


  # calculated RKappa
  Rkappa <- diag(sum(sizes))
  starts <- cumsum(c(1, sizes))
  starts <- starts[1:m]

  for (i in 2:m){
    for (j in 1:(i-1)){
      newbottom <- diag(Drs[[i]], nrow = length(Drs[[i]])) %*% crossprod(Us[[i]], Us[[j]]) %*% diag(Drs[[j]], nrow = length(Drs[[j]]))
      ran1 <- starts[i]:(starts[i] + sizes[i] - 1)
      ran2 <- starts[j]:(starts[j] + sizes[j] - 1)
      Rkappa[ran1, ran2] <- newbottom
      Rkappa[ran2, ran1] <- t(newbottom)
    }
  }

  # output details
  D <- det(Rkappa)
  J <- -.5 * log(D)

  return(J)

}

chol_gauss <- function(x, sigma, tol){
  # incomplete cholesky decomposition of the gram matrix defined by data x

  n <- dim(x)[2]
  Pvec <- 1:n
  I <- integer(0)

  # calculates diagonal elements (all equal to 1 for gaussian kernels)
  diagG <- rep(1, n)
  i = 1
  G <- integer(0)

  while ((sum(diagG[i:n])) > tol) {
    G <- cbind(G, rep(0, n))

    # find best new element
    if (i > 1){

      diagmax <- max(diagG[i:n])
      jast <- which.max(diagG[i:n])
      jast <- jast + i - 1

      # updates permutation
      Pvec[c(i, jast)] <- Pvec[c(jast, i)]

      # updates all elements of G due to new permutation
      G[c(i, jast), 1:i] <- G[c(jast, i), 1:i, drop = FALSE]

      # do the cholesky update

    } else {
      jast <- 1
    }

    G[i, i] <- diagG[jast]
    G[i, i] <- sqrt(G[i, i, drop = FALSE])

    if (i < n){
      # calculates newAcol = A[Pvec[(i+1):n], Pvec[i]]
      newAcol = exp(-.5 / sigma^2 * sqdist(x[, Pvec[(i + 1):n], drop = FALSE],
                                           x[, Pvec[i], drop = FALSE]))

      if (i > 1){
        G[(i + 1):n, i] <- 1 / G[i, i] *
          (newAcol - G[(i + 1):n, 1:(i - 1), drop = FALSE] %*%
             t(G[i, 1:(i - 1), drop = FALSE]))
      } else {
        G[(i + 1):n, i] <- 1 / G[i, i] * newAcol
      }
    }

    # updates diagonal elements
    if (i < n){
      diagG[(i + 1):n] <- rep(1, n - i) - rowSums(G[(i + 1):n, 1:i, drop = FALSE]^2)
    }

    i <- i + 1

  }

  # return results
  res <- list()
  res$G <- G
  res$Pvec <- matrix(Pvec - 1, nrow = 1)

  return(res)
}

sqdist <- function(a, b){
  # computes squared euclidean distance matrix
  # computes a rectangular matrix of pairwise distances
  # between points in A (given in columns) and points in B
  # NOTE: very fast implementation taken from Rolang Bunschoten

  aa <- matrix(colSums(a * a), nrow = 1)
  bb <- matrix(colSums(b * b), nrow = 1)
  ab <- t(a) %*% b
  d <- abs(t(aa) + matrix(bb, nrow = dim(aa)[2]) - 2 * ab)

  return(d)
}

