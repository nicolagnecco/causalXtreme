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

# are randomDAG and random_dag similar?
test <- matrix(0, nrow = 1e3, ncol = 2)

for (i in 1:1e3){
  # u <- sample(1e6, 1)
  # set.seed(u)
  d <- randomDAG(30, 2/29)
  # set.seed(u)
  d2 <- random_dag(30, 2/29)
  test[i, 1] <- sum(d)
  test[i, 2] <- sum(d2)
}
plot(test[, 1])
plot(test[, 2])
plot(test[, 1] - test[, 2])
sum(test[, 1] < test[, 2])
sum(test[, 1] == test[, 2])
sum(test[, 1] > test[, 2])
mean(test[, 1])
mean(test[, 2])


# inverse mirror uniform
mirror_values_2 <- function(min, max){
  sample(c(-1,1), size=1, replace = T) * runif(1, min, max)
}
n <- 1e4
mirror_values <- mirror_values2 <-  numeric(n)
uniform <- numeric(n)
for (i in 1:n){
  mirror_values[i] <- inverse_mirror_uniform(runif(1), 0.3, 0.9)
  mirror_values2[i] <- mirror_values_2(0.3, 0.9)
  uniform[i] <- runif(1)
}
mirror_values <- sapply(runif(1e4), inverse_mirror_uniform, .3, .9)

hist(mirror_values)
hist(mirror_values2)
hist(uniform)
mean(mirror_values > 0)
mean(mirror_values2 > 0)
