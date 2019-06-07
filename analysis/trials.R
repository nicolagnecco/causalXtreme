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

