# Multidimensional Item Response Theory (MIRT) simulation code
# -------------------------------------------------------------
# install.packages("mirt")
library(mirt)
# -------------------------------------------------------------

# Simulate data with Monte Carlo experiments ------------------------------
theta <- rnorm(n = 1000, mean = 0, sd = 1) # person trait level, e.g., math skill

# Item parameters (a, b, c) for 4 items
a <- c(0.5, 1, 1.5, 2)
b <- c(0, -1, 1, 0)
c <- c(0.2, 0.15, 0.25, 0.2)

n.persons <- length(theta)
n.items <- length(a)
response.data <- matrix(NA, n.persons, n.items)

for (i in 1:n.persons) {
    for (j in 1:n.items) {
        p <- c[j] + (1 - c[j]) / (1 + exp(-(a[j] * theta[i] - b[j])))
        u <- runif(n = 1, min = 0, max = 1)
        if (u < p) {
            response.data[i, j] <- 1
        } else {
            response.data[i, j] <- 0
        }
    }
}
# another option to simulate data
# simdata <- simdata(a, b, N=500, itemtype= '2PL') # simulating responses from 500 individuals

colnames(response.data) <- c("I1", "I2", "I3", "I4")

# Calibrate item parameters using the mirt package, default settings ------
# 2PL model
mirt.out <- mirt::mirt(data = response.data, model = 1, itemtype = "2PL")

# mirt.out =  results from mirt analysis
# data = response data is the data matrix
# model = 1 refers to the model should have 1 dimension
# itemtype = '2PL' specifies that the item in the test should be a 2 parameter logistic model (discrimination, difficulty)

# See the estimated item parameters
mirt::coef(mirt.out, simplify = T, IRTparts = T)

# coef() function extracts the estimated parameters from the mirt analysis
# simplify=T simplifies the output if possible (into vector or matrix), otherwise a list
# IRTparts = T specifies the output should be a IRT parameter matrix (T = discrimination, F = original parameterization)

## What are the parameters?
# discrimination parameter (a) and difficulty parameter (b)
# (a) describes how well an item (e.g. question on test) can differentiate between individuals with different levels of the underlying trait that the test is measuring
# (b) describes the level of the underlying trait that is required to have a 50% chance of answering the item correctly
# other paramaters: 
# (g) guessing parameter (measures prob. that a person with a very low level of underlying trait will respond positively/ aka. get the question right)
# (u) upper asymptote parameter (measures the maximum probability that a individual will respond positively(get question right), regardless of underlying trait)
# Note: higher the upper asymptote, the easier the item is.

# OUTPUT:
# $items = df with estimated item parameters, row = items, col = parameters
# (a1 = discrimination, d = difficulty, g = guessing, u = upper asymptote)
# $means = vector that contains the estimated mean of the latent trait (F1 in this case with a mean of 0)
# $cov = matrix that contains the estimated covariance matrix of the latent trait (F1 in this case with a variance of 1)
# Note: only trait F1 in this case so $cov is a 1x1 matrix.

# Calculate estimation accuracy -------------------------------------------

a.est <- mirt::coef(mirt.out, simplify = T, IRTparts = T)$items[, 1]
RMSE.a <- sqrt(mean((a - a.est)^2))
RMSE.a

b.est <- mirt::coef(mirt.out, simplify = T, IRTparts = T)$items[, 2]
RMSE.b <- sqrt(mean((b - b.est)^2))
RMSE.b

c.est <- mirt::coef(mirt.out, simplify = T, IRTparts = T)$items[, 3]
RMSE.c <- sqrt(mean((c - c.est)^2))
RMSE.c

# What does this mean?
# a.est extracts the estimated a parameter from the fitted irt model from (mirt.out)
# RMSE.a = calculates the root mean squared error between the true a parameter and the estimated a parameter
# (measure of differences/error between the estimated parameter values and true parameter)

# Understand the 2PL Model -------------------------------------------
# 2PL model
mirt.out <- mirt::mirt(data = response.data, model = 1, itemtype = "2PL")

# mirt.out =  results from mirt analysis
# data = response data is the data matrix
# model = 1 refers to the model should have 1 dimension
# itemtype = '2PL' specifies that the item in the test should be a 2 parameter logistic model (discrimination, difficulty)

# What does model = 1 mean?
# "dimension" refers to the underlying trait/ability that the model is trying to test
# e.g. math skill, the model is measuring this, and the model is trying to estimate the math skill of each individual

# What does itemtype = '2PL' mean?
# '2PL' specifies that the item in the test should be a 2 parameter logistic model (discrimination, difficulty)

# What does the 2PL model look like?
# 2PL model: P(Xij = 1) = cij + (1 - cij) / (1 + exp(-(aij * theta - bij))) 
# P(Xij = 1) = probability of person i answering item j correctly
# cij = guessing parameter (measures prob. that a person with a very low level of underlying trait will respond positively/ aka. get the question right)
# aij = discrimination parameter (describes how well an item (e.g. question on test) can differentiate between individuals with different levels of the underlying trait that the test is measuring)
# bij = difficulty parameter (describes the level of the underlying trait that is required to have a 50% chance of answering the item correctly)
# theta = person trait level, e.g., math skill
