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
start.time = Sys.time() #start time
mirt.out <- mirt::mirt(data = response.data, model = 1, itemtype = "2PL", storeEMhistory=TRUE) # 2PL model with history tracking
end.time = Sys.time() # end time
time.to.conv = end.time - start.time

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

df$time.to.conv <- time.to.conv

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

# ICC shows the probability of a positive response to each item as a function of the latent trait. 
# x-axis = latent trait (theta) 
# y-axis = probability of positive response (get question right)
# progression of the curves across the iterations of the EM algorithm gives an idea of how the parameter estimates for each item are converging 


# Plot test information fucntion (TIF) ----------------------------------------------------------

plot(mirt.out, type = "info", auto.key =  list(columns = 4, legend = TRUE), 
    main = "Test Information Function (TIF)")

# Test Information Funciton (TIF)
# TIF provides information about the reliability of a test at different levels of a trait being measured.
# it is a functin of the item parameters and indicates where on the trait scale the test is most imformative.

# What is information?
# "information" is inversely related to the variance of the estimated ability level (theta) 
# the more information provided by a test means more percision in meaning the trait at a particular level 

# What am I looking at?
# x-axis = latent trait (theta)
# y-axis = information (reliability) of the test at a particular level of the latent trait

# Why is this useful? 
# TIF is useful for test developers to evaluate and improve tests, specifically identifying which parts of the trait continuum are well-assessed and which are not.
# test with a TIF has a high peak at certain ability leels but low information at others might be very effective for examinees around that peak, but less reliable elsewhere.
# one might aim for a TIF that provides a uniform level of high information across the range of ability levels.


# Plot person-item map ----------------------------------------------------------

# retrieving coef()$items difficutly estimates
item.difficulties <- coef(mirt.out, simplify = TRUE)$items[, 2]
item.difficulties <- data.frame(Item = paste("Item", 1:4), Difficulty = item.difficulties)

# retrieving f-scores (person ability estimates)
person.abilities <- fscores(mirt.out)
person.abilities <- data.frame(Ability = person.abilities)

item.difficulties.long <- item.difficulties[rep(seq_len(nrow(item.difficulties)), each = nrow(person.abilities)), ]
item.difficulties.long <- item.difficulties[rep(1:nrow(item.difficulties), each = length(person.abilities)), ]

# difficulties with person abilities
plot_data <- data.frame(
    Item = item.difficulties.long$Item,
    Difficulty = item.difficulties.long$Difficulty,
    Ability = person.abilities
  )

ggplot(plot_data, aes(x = F1, y = Difficulty, color = Item)) +
    geom_point() +
    labs(title = "Person-Item Map", x = "Person Ability", y = "Item Difficulty") +
    theme_minimal() +
    scale_color_brewer(palette = "Set1")

# Person-Item Map
# this map usually represents the distribution of person abilities in relation to item difficulties.
# It helps in understanding how well the item on a test are targeted to the ability of the test takers.

# Plot factor loadings ----------------------------------------------------------

# retrieves items and parameters
loadings <- coef(mirt.out, simplify = TRUE)$items[, "a1"]
loadings <- as.data.frame(loadings)
colnames(loadings) <- c("a1")

# a1 (discrimination parameter) 
# - seen as multiplier that links the latent trait (ability/theta) to the probability of a specific response pattern
# - indicates how effectively each item discriminates between two different lavels of the latent trait
# d (difficulty parameter) 

# 2PL model in IRT (where c=d=0) is a close cousin to confimatory factor analysis (CFA) with binary variables
# c = d = 0 (no guessing or upper asymptote) 

library(tidyr)
loadings.long <- loadings %>%
    tibble::rownames_to_column(var = "Item") %>%
    tidyr::pivot_longer(cols = -Item, names_to = "Parameter", values_to = "Value")

ggplot(loadings.long, aes(x = Item, y = Value, fill = Parameter)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = "Factor Loadings", x = "Item Parameters", y = "Factor Loadings") +
    theme_minimal() +
    scale_fill_brewer(palette = "Set1")

# Factor Loadings
# In multidimensional IRT, factor loadings describe the relationship(correlation) between the items and latent traits.
# High factor loadings indicate a storong association betwen an item and a particular trait.

# What Am I looking at?
# x-axis = item parameters (I1, I2, I3, I4) and (a, b, c) for each 
# y-axis = factor loadings (correlation) between the items and latent traits (F1, F2, F3, F4)

# Plot residual analysis ----------------------------------------------------------

residuals <- residuals(mirt.out)
residuals <- data.frame(residuals)
residuals <- tibble::rownames_to_column(residuals, var = "Item")
residuals.long <- reshape2::melt(residuals, id.vars = "Item")
colnames(residuals.long) <- c("Var1", "Var2", "Residual")  

ggplot(residuals.long, aes(x = Var1, y = Residual, color = abs(Residual))) +
    geom_point(alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(title = "Residual Analysis", x = "Item", y = "Residual") +
    theme_minimal() +
    scale_color_gradient(low = "blue", high = "red")

# Residual Analysis
# Examining the differences between the observed item responses and the responses predicted by the model.
# It helps identifying items that do not fit well in the model, indicating potential issues with the item/assumptions.
