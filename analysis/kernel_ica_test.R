library(kernlab)

# function contrast_ica ####
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
    res <- chol_gauss(x[i, , drop = FALSE] / kparam$sigmas[i], 1, N * etas[i])
    G <- res$G
    # file_nm <-  paste("analysis/G", i, ".csv", sep = "") # !!!
    # G <- as.matrix(read.csv(file_nm, header = FALSE)) # !!!
    Pvec <- res$Pvec

    k <- NCOL(G)
    a <- sort(Pvec)
    Pvec <- order(Pvec)

    G <- centerCols(G[Pvec, , drop = FALSE])

    # regularization (see paper for details)
    eigen_decomp <- eigen(t(G) %*% G)
    A <- eigen_decomp$vectors[, k:1, drop = FALSE]
    D <- eigen_decomp$values[k:1]

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
      newbottom <- diag(Drs[[i]], nrow = length(Drs[[i]])) %*% (t(Us[[i]]) %*% Us[[j]]) %*% diag(Drs[[j]], nrow = length(Drs[[j]]))
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

# function chol_gauss ####
centerCols <- function(X){
  n <- NROW(X)
  p <- NCOL(X)
  X - matrix(rep(colMeans(X), n), ncol = p, byrow = TRUE)
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
  res$Pvec <- Pvec

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
