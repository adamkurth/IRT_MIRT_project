# Multidimensional Item Response Theory (MIRT) simulation code

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

# "mirt.out" =  results from mirt analysis
# "data = response data" is the data matrix
# "model = 1" refers to the model should have 1 dimension
# "itemtype = '3PL'" specifies that the item in the test should be a 3 parameter logistic model

# See the estimated item parameters
mirt::coef(mirt.out, simplify=T, IRTparts=T)

# "coef()" function extracts the estimated parameters from the mirt analysis
# "simplify=T" simplifies the output if possible (into vector or matrix), otherwise a list
# " IRTparts = T " specifies the output should be a IRT parameter matrix (T = discrimination, F = original parameterization)

# OUTPUT:
# $items = df with estimated item parameters. row = items, col = parameters
  # (a1 = discrimination, d = difficulty, g = guessing, u = upper asymptote)

# $means = vector that contains the estimated mean of the latent trait (F1 in this case with a mean of 0)

# $cov = matrix that contains the estimated covariance matrix of the latent trait (F1 in this case with a variance of 1) 

# Note: only trait F1 in this case so $cov is a 1x1 matrix.
 

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

# "RMSE" = root mean square error


# Plot the estimated item parameters --------------------------------------
library(ggplot2)

plot.data <- data.frame(
  Item = rep(c("I1", "I2", "I3", "I4"), each = 2),
  Parameter = rep(c("a", "b", "c"), times = 4),
  Value = c(a, b, c, a.est, b.est, c.est),
  Type = rep(c("True", "Estimated"), each = 12)
)

ggplot(plot.data, aes(x = Item, y = Value, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Parameter, scales = "free") +
  theme_minimal() +
  labs(x = "Item", y = "Parameter Value", fill = "Type",
       title = "True vs. Estimated Item Parameters")

