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
    # !!!
    res <- chol_gauss(x[i, , drop = FALSE] / kparam$sigmas[i], 1, N * etas[i]);
    G <- res$G
    Pvec <- res$Pvec
    # !!!

    k <- NCOL(G)
    a <- sort(Pvec)
    Pvec <- order(Pvec)


    G <- drop(scale(G[Pvec,], scale = FALSE))
    attr(G, "scaled:center") <- NULL

    # regularization (see paper for details)
    eigen_decomp <- eigen(t(G) %*% G)
    D <- eigen_decomp$values[k:1]
    A <- eigen_decomp$vectors[k:1, k:1]
    indexes <- which(D >= N * etas[i]);  #removes small eigenvalues
    newinds <- sort(D[indexes])
    neworder <- order(D[indexes])
    neworder2 <- neworder[length(neworder):1]
    neig <- length(indexes)

    indexes <- indexes[neworder2(1:neig)]
    if (identical(idexes, integer(0))){
      indexes = 1
    }

    D <- D[indexes]
    V = G %*% (A[, indexes] %*% diag(sqrt(1 / D)))

    Us[[i]] <- V
    Lambdas[[i]] <- D
    Dr <- D

    for (j in seq_along(D)){
      Dr[j] <- D[j] / (N * kappas[i] + D[j])
    }

    Drs[[i]] <- Dr
    sizes[i] <- dim(Drs[[i]][1])

  }


  # calculated RKappa
  Rkappa <- diag(sum(sizes))
  starts <- cumsum(c(1, sizes))
  starts <- starts[1:m]

  for (i in 2:m){
    for (j in 1:(i-1)){
      newbottom <- diag(Drs[[i]]) %*% (t(Us[[i]]) %*% Us[[j]]) %*% diag(Drs[[j]])
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
chol_gauss <- function(x_mat, sigma, tol){
  # !!! to code
}


# Tests ####
x <- read.csv("analysis/X_mat.csv", header = F)
x <- t(x)



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

J = contrast_ica(contrast,x,kparam);

