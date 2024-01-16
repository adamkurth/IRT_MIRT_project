

# True values of model paramters ------------------------------------------

# Person trait level, e.g., math skill
theta <- rnorm(n = 1000, mean = 0, sd = 1)

# Item parameters (a, b, c) for 4 items
a <- c(0.5, 1, 1.5, 2)
b <- c(0, -1, 1, 0)
c <- c(0.2, 0.15, 0.25, 0.2)

# Simulate data with Monte Carlo experiments ------------------------------

n.persons <- length(theta)
n.items <- length(a)
response.data <- matrix(NA, n.persons, n.items)
  
for (i in 1:n.persons){
  for (j in 1:n.items){
    p <- c[j] + (1 - c[j]) / (1 + exp(-(a[j] * theta[i] - b[j])))
    u <- runif(n = 1, min = 0, max = 1)
    if (u < p){
      response.data[i, j] <- 1
    } else {
      response.data[i, j] <- 0
    }
  }
}

colnames(response.data) <- c("I1", "I2", "I3", "I4")

# Calibrate item parameters using the mirt package, default settings ------

mirt.out <- mirt::mirt(data = response.data, model = 1, itemtype = '3PL')

# See the estimated item parameters
mirt::coef(mirt.out, simplify=T, IRTparts=T)

# Calculate estimation accuracy -------------------------------------------

a.est <- mirt::coef(mirt.out, simplify=T, IRTparts=T)$items[, 1]
RMSE.a <- sqrt(mean((a - a.est) ^ 2))
RMSE.a

b.est <- mirt::coef(mirt.out, simplify=T, IRTparts=T)$items[, 2]
RMSE.b <- sqrt(mean((b - b.est) ^ 2))
RMSE.b

c.est <- mirt::coef(mirt.out, simplify=T, IRTparts=T)$items[, 3]
RMSE.c <- sqrt(mean((c - c.est) ^ 2))
RMSE.c