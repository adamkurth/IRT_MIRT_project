# install.packages("mirt")
library(mirt)

# Simulate some data ------------------------------------------------------
set.seed(123) # set seed for reproducibility
a <- matrix(rlnorm(20), ncol=1) # discrimination parameters
b <- matrix(rnorm(20), ncol=1) # difficulty parameters
simdata <- simdata(a, b, N=500, itemtype= '2PL') # simulating responses from 500 individuals

# 2PL (2 parameter logistic model)
# 0 or 1 for 20 items with 2 parameters logistic model (discrimination and difficulty)

# Fit the model -----------------------------------------------------------
model <- mirt(simdata, 1, itemtype='3PL') # 1-factor model with 3PL items
summary(model)

# in summary we find...
# Factor loadings (F1)

# Plot the item information curves ----------------------------------------
plot(model, type='trace', item=1:5) # plot trace lines for first 5 items.
plot(model, type='info', item=1:5) # plot item information curves for first 5 items.

abilities <- fscores(model) # get factor scores
head(abilities) # show top of factor data

