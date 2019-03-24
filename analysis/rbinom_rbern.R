n <- 10
p <- 0.3
n_exp <- 1e3



rbin <- numeric(n_exp)
for (i in 1:n_exp){
  rbin[i] <- rbinom(n = 1, size = n, prob = p)
}

rbern <- matrix(0, nrow = n_exp, ncol = n)
for (i in 1:n_exp){
  rbern[i, ] <- rbinom(n = n, size = 1, prob = p)
}

hist(rbin, 5)
hist(apply(rbern, 1, sum), 5)
summary(rbin)
summary(apply(rbern, 1, sum))


difference <- numeric(100)
for (j in 1:100){
rbin <- matrix(0, nrow = n_exp, ncol = n)
for (i in 1:n_exp){
  r <- rbinom(n = 1, size = n, prob = p)
  rbin[i, sample(x = n, size = r, replace = FALSE)] <- 1
}

rbern <- matrix(0, nrow = n_exp, ncol = n)
for (i in 1:n_exp){
  rbern[i, ] <- rbinom(n = n, size = 1, prob = p)
}

difference[j] <- mean(apply(rbin, 2, sum)/ n_exp) -
  mean(apply(rbern, 2, sum)/ n_exp)
}
plot(difference, type = 'l')

