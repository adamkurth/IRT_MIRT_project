# Multidimensional Item Response Theory (MIRT) simulation code (2PL)
# -------------------------------------------------------------
rm(list=ls())
# install.packages("mirt")
library(mirt)
library(ggplot2)
library(reshape2)
# -------------------------------------------------------------

# Simulate data with Monte Carlo experiments ------------------------------
# person trait level, e.g., math skill 
theta <- rnorm(n = 1000, mean = 0, sd = 1) 

# Item parameters (a, b, c) for 4 items
a <- c(0.5, 1, 1.5, 2)
b <- c(0, -1, 1, 0)
c <- c(0.2, 0.15, 0.25, 0.2)

n.persons <- length(theta)
n.items <- length(a)
response.data <- matrix(NA, n.persons, n.items)

for (i in 1:n.persons) {
    for (j in 1:n.items) {
        p <- c[j] + (1 - c[j]) / (1 + exp(-(a[j] * theta[i] - b[j]))) #2PL model
        u <- runif(n = 1, min = 0, max = 1) 
        if (u < p) {
            response.data[i, j] <- 1
        } else {
            response.data[i, j] <- 0
        }
    }
}

# a, b, c are the item parameters
# theta = personal trait level
# p = prob. of positive response (get question right)
# u = random number between 0 and 1 (unifiorm distribution)
# if u < p, then response.data = 1 (get question right)
# if u > p, then response.data = 0 (get question wrong)

colnames(response.data) <- c("I1", "I2", "I3", "I4")

# Calibrate item parameters using the mirt package, default settings ------
# record time to convergence and estimation accuracy ---------------------
start_time = Sys.time() #start time
mirt.out <- mirt::mirt(data = response.data, model = 1, itemtype = "2PL", storeEMhistory=TRUE) # 2PL model with history tracking
end_time = Sys.time() # end time
time_to_conv = end_time - start_time

# mirt.out <- mirt::mirt() function fits a maximum likelihood (posterier) factor analysis model to the data (dichotomus/polytomous)
# output: 
# 'M-step' = indicates that the Broyden-Fletcher-Goldfarb-Shanno (BFGS) algorithm used for optimization in M-step if the EM algorithm.
# 'EM acceleration' = indicates that ramsay's acceleration method is used to speed up the convergence of the EM algorithm.
# 'Number of rectacgular quadrature' = the number of quadrature points used in the numerical integration of the likelihood function.
# 'Latent density type' = the distribution of the latent trait (F1 in this case)

# log-likelihood = how well the model fits the data
# est parameters = number of estimated parameters 
# AIC = Akaike information criterion (measure of model fit, takes into number of parameters)
# BIC/SABIC = Bayesian information criterion/Sample-size Adjusted BIC (measure of model fit, takes into number of parameters)
# G2 (3) = 1.1, p = 0.7772 : results of the likelihood ratio test comparing the fitted model against saturated model 
# RMSEA = 0, CFI = NaN, TLI = NaN : additional measures of model fit, the Root Mean Square Error of Approximation (RMSEA),
#  the Comparative Fit Index (CFI), and the Tucker-Lewis Index (TLI)


# See the estimated item parameters
mirt::coef(mirt.out, simplify = T, IRTparts = T)

# coef() function extracts the estimated parameters from the mirt analysis
# simplify=T simplifies the output if possible (into vector or matrix), otherwise a list
# IRTparts = T specifies the output should be a IRT parameter matrix (T = discrimination, F = original parameterization)

# OUTPUT:
# $items = df with estimated item parameters, row = items, col = parameters
# (a1 = discrimination, d = difficulty, g = guessing, u = upper asymptote)
# $means = vector that contains the estimated mean of the latent trait (F1 in this case with a mean of 0)
# $cov = matrix that contains the estimated covariance matrix of the latent trait (F1 in this case with a variance of 1)
# Note: only trait F1 in this case so $cov is a 1x1 matrix.

## What are the parameters?
# discrimination parameter (a) and difficulty parameter (b)
# (a) describes how well an item (e.g. question on test) can differentiate between individuals with different levels of the underlying trait that the test is measuring
# (b) describes the level of the underlying trait that is required to have a 50% chance of answering the item correctly
# other paramaters: 
# (g) guessing parameter (measures prob. that a person with a very low level of underlying trait will respond positively/ aka. get the question right)
# (u) upper asymptote parameter (measures the maximum probability that a individual will respond positively(get question right), regardless of underlying trait)
# Note: higher the upper asymptote, the easier the item is.

# Calculate estimation accuracy and time to convergence -------------------------------------------
a.est <- mirt::coef(mirt.out, simplify = T, IRTparts = T)$items[, 1]
RMSE.a <- sqrt(mean((a - a.est)^2))

b.est <- mirt::coef(mirt.out, simplify = T, IRTparts = T)$items[, 2]
RMSE.b <- sqrt(mean((b - b.est)^2))

c.est <- mirt::coef(mirt.out, simplify = T, IRTparts = T)$items[, 3]
RMSE.c <- sqrt(mean((c - c.est)^2))

# Store RMSE values and time to convergence in dataframe
df <- data.frame("a.discrim" = a, "b.difficulty" = b, "c.guessing" = c, 
                "a.est" = a.est, "b.est" = b.est, "c.est" = c.est,
                "RMSE.a" = RMSE.a, "RMSE.b" = RMSE.b, "RMSE.c" = RMSE.c)

df$time_to_conv <- time_to_conv

# to extract history use extract.mirt() function
# params <- extract.mirt(mirt.out, 'parvec') #parameter vector


# plots:
# -log-likelihood
# -item characteristic curves (ICC) 
# -test information fucntion (TIF)
# -person-item map
# -factor loadings
# -residual analysis

# Plot log-likelihood ----------------------------------------------------------

ll.history <- extract.mirt(mirt.out, 'LLhistory')

plot(ll.history, type="b", col='blue', ylim = range(ll.history), 
    xlab = "Iteration", ylab = "Log-Likelihood", main = "Log-Likelihood Convergence")

# Plot item characteristic curves (ICC) ----------------------------------------------------------

plot(mirt.out, type = "trace", auto.key =  list(columns = 4, legend = TRUE), 
     main = "Item Characteristic Curves (ICC)")

itemplot(mirt.out, item = 1, type = "trace", auto.key =  list(columns = 4, legend = TRUE), 
         main = "Item Characteristic Curves (ICC) for Item 1")
itemplot(mirt.out, item = 2, type = "trace", auto.key =  list(columns = 4, legend = TRUE), 
         main = "Item Characteristic Curves (ICC) for Item 2")
itemplot(mirt.out, item = 3, type = "trace", auto.key =  list(columns = 4, legend = TRUE),
         main = "Item Characteristic Curves (ICC) for Item 3")
itemplot(mirt.out, item = 4, type = "trace", auto.key =  list(columns = 4, legend = TRUE), 
         main = "Item Characteristic Curves (ICC) for Item 4")

# Plot test information fucntion (TIF) ----------------------------------------------------------

plot(mirt.out, type = "info", auto.key =  list(columns = 4, legend = TRUE), main = "Test Information Function (TIF)")

# Plot person-item map ----------------------------------------------------------

plot(mirt.out, type = "map", auto.key =  list(columns = 4, legend = TRUE), main = "Person-Item Map")


# Plot factor loadings ----------------------------------------------------------

plot(mirt.out, type = "loadings", auto.key =  list(columns = 4, legend = TRUE), main = "Factor Loadings")

# Plot residual analysis ----------------------------------------------------------

plot(mirt.out, type = "residuals", auto.key =  list(columns = 4, legend = TRUE), main = "Residual Analysis")


