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


# First, 
# Overarching Understanding ----------------------------------------------------------
# - (Objective) the goal is to estimate the parameters (a, b, c) of the test items based on the repsponses of the individuals. 
#   these parameters are then used to understand how different items funciton in measuring latent traits (ability/aptiude)
# - (Data Structure) your dataset where each row is an individual and each column is an item (question on test)
# - (Model) model fitting involves finding the best set of item parameters that explain the observed responses, this is where we optimize. 

# Second,
# Optimization in MIRT --------------------------------------------------------------- 
# Optimizing is at the heart of MIRT model, the idea is to find the paramters that maximize the likelihood of observing the data.
# - (Log-Likelihood) Statistical measure used to assess how well the model (with its given parameters) explains the data. 
#   Aim is to maximize the log-likelihood.
# - (Algorithm) like the Newton-Raphson algorithm, the EM algorithm is used for optimization. These algorithms iteratively adjust parameters to find max log-likelihood.
#   Each iteration involves calculating the log-likihood and adjusting the parameters (to maximize) based on the algorithms rules. 

# Third,
# Key Arguments in `mirt` -----------------------------------------------------------------------
# - (method) this argument type specifies the method to be used. 
#   Common methods include 'MH-RM' (Metropolis-Hastings Robbins-Monro), 'EM' (Expectation-Maximization), etc. 
# - (optimizer) this argument type specifies the numerical iterative optimizer to be used, such as Newton-Raphson, or 'BDGS', etc.
#   These algorithms have different wasy of navigating the parameter space to find the max log-likelihood.
# - (dentype) this argument refers to the density type used in the estimation, especially relaven for certain estimation methods,
#   such as 'Gaussian' (multivariate Gaussian dist), 'empiricalhist' (empirical histogram), 'Davidian-#' (semi-parametric Davidian curves), etc.
#   It affects how the log-likihood is calculated, and can influence the optimizaiton process.

# Summary
# when fitting a MIRT model, you're engaging in an optimzation process where you're looking to find the best set of item parameters to explain the data.
# This involves choosing appropriate estimation methods, an optimizer for finding the max log-likelihood, and a density type for calculating likihoods.
# Each choice affects the way the model navigates the complex parameter space to arrive at the most accurate, and reliable  estimates of item parameters.

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

# For further understanding of (1PL, 2PL, 3PL), see practice/n_PL.md file in this repo. 

# Understanding input arguments of mirt::mirt function -----------------

# further documentaiton please see practice/objective.md
?mirt::mirt # for more info on the mirt package

# some arguments
# data = 'martix' or 'df' that consists of numerically ordered data
# model = string or object from mirt.model(), that declares the model to be estimated '1', '2', etc.
#   this is the number of dimensions in the model meaning that the model is trying to estimate the ability of the individual in that dimension
# itemtype = declared as vector, type of item to be modeled, '2PL', '3PL', etc there are many more.
# ...  

# Alerable Arguments ---------------------------------------------------
# 1. method = a char object specifying estimation algorithm to be used. 
#   'EM' = expectation maximization algorithm (default)
#    - generally effective for 1-3 factors
#   'QMCEM' = quasi-Monte Carlo EM algorithm,  
#   'MCEM' = Monte Carlo EM algorithm,
#   'MHRM' = marginal maximum likelihood estimation, 
#   'SEM' = stochastic EM algorithm (first two stages of the MH-RM algorithm stage using optimizer other than sinple Newton-Raphson),
#   'BFGS' = Broyden-Fletcher-Goldfarb-Shanno algorithm
#  Note: 'QMCEM'’, ‘'MCEM'’, ‘'SEM'’, or ‘'MHRM'’ should be used when the dimensions are 3 or more.

# 2. optimizer = a char object specifying the numerical iterative optimizer to be used.
#  'BFGS' = Broyden-Fletcher-Goldfarb-Shanno algorithm (default, no upper/lower bound "box"-constraints)
#  'NR' = Newton-Raphson algorithm (more effective, but not as stable for more complex models e.g. nominal response model, or nested logit model)
#  'NR1' = Newton-Raphson algorithm but consists of only 1 update coupled with RM  Hessian (only applicable when MH-RM algorithm is used)
#  'Nelder-Mead' = Nelder-Mead algorithm (simplex algorithm, no gradient information required),

# 3. dentype = type of density form to use for latent trait parameters.
#  'Gaussian' = assumes multivariate Gaussian dist, with associated mean vector and variance-covariance matrix (default)
#  'empiricalhist' or 'EH' = estimates the latent distribution using empirical hist. Only applicable for unidimensional model estimated with the EM algorithm.
#  'empiricalhist_Woods' or 'EHW' = as above, but Woods' correction is applied to the variance estimate. 
#  'Davidian-#' = estimates semi-parametric Davidian curves, where # is the number of Davidian paramters to estimate. 
#       e.g. 'Davidian-2' estimates a 2-parameter Davidian curve.

# Report RMSE ----------------------------------------------------------

