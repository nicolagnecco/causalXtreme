# Understand causal tail coefficient with both tails
n <- 11
x1 <- runif(n)
r1 <- rank(x1)
x2 <- runif(n)
r2 <- rank(x2)
k <- 2
plot(r1 > n - k/2 | r1 <= k/2)
plot(abs(2 * r1 - n) > n - k)
plot(r1 > n - k)
sum(r1 > n - k / 2 | r1 <= k / 2)
table(abs(r2 - (n + 1)/2))
table(r2)
1 / (k * n) * sum(2 * abs(r2[r1 > n - k / 2 | r1 <= k / 2] - (n + 1)/2))

set.seed(1991)
n <- 5e2
X1 <- rt(n, df = 2.5)
X2 <- X1 + rt(n, df = 2.5)
X3 <- X1 + rt(n, df = 2.5)
dat <- cbind(X1, X2, X3)
compute_gamma_matrix(dat, both_tails = FALSE)
compute_gamma_matrix(dat, both_tails = TRUE)
pairs(dat)

# Rescale betas so that SNR = (1 - w)/w
p <- 3
g <- rbind(c(0, 1, 1), c(0, 0, 1), c(0, 0, 0))
A <- random_coeff(g, two_intervals = FALSE)
B <- t(A)
path <- solve(diag(p) - B)
has_parent <- apply(g, 2, sum) != 0
w <- .4
diag(path)[has_parent] <- sqrt(w)

path_temp <- path^2
diag(path_temp) <- 0
s <- apply(path_temp, 1, sum)
s[s != 0] <-  sqrt(1 / (s[s != 0]) * (1 - w))
s[s == 0] <- 1
mm <- matrix(rep(s, p), ncol = p)
diag(mm) <- 1
out <- path * mm
apply(out^2, 1, sum)
A_mod <- t(diag(p) - solve(out))
A
A_mod
