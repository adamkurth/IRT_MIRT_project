# Multidimensional Item Response Theory (MIRT) simulation code (3PL)
# -------------------------------------------------------------
rm(list=ls())
# install.packages("mirt")
library(mirt)
library(ggplot2)
library(reshape2)
library(tibble)
# -------------------------------------------------------------

# Simulate data with Monte Carlo experiments ------------------------------
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
        p <- c[j] + (1 - c[j]) / (1 + exp(-(a[j] * theta[i] - b[j]))) #3PL model
        u <- runif(n = 1, min = 0, max = 1) 
        if (u < p) {
            response.data[i, j] <- 1
        } else {
            response.data[i, j] <- 0
        }
    }
}

colnames(response.data) <- c("I1", "I2", "I3", "I4")

# Calibrate item parameters using the mirt package, default settings ------
start.time = Sys.time() #start time
mirt.out <- mirt::mirt(data = response.data, model = 1, itemtype = "3PL", storeEMhistory=TRUE) # 2PL model with history tracking
end.time = Sys.time() # end time
time.to.conv = end.time - start.time

# See the estimated item parameters
mirt::coef(mirt.out, simplify = T, IRTparts = T)

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

plot(mirt.out, type = "info", auto.key =  list(columns = 4, legend = TRUE), 
    main = "Test Information Function (TIF)")

# This plot differs from the 2PL model because the 3PL model has a guessing parameter (c) that is not 0.
# This is the guessing parameter, which is not 0, and thus creates sharp peaks in the TIF.
# The TIF tells us that the test is most informative at the mean of the latent trait (theta = 0).
# The TIF tells us that the test is least informative at the extremes of the latent trait (theta = -3 and theta = 3).

# The guessing parameter 'c' in the 3PL model affects the floor of the TIF,
# This does not introduce variabillity in the sense of making the TIF more varaible at different points...
# instead, this raises the information floor, meaning that even at lower ability (where discrimination is typically low),
# the test provides some information due to the chance of guessing correctly. This is because of the guessing parameter 'c',
# ensures that the probability of a correct response is never 0, even for low-ability examinees.

plot(mirt.out, type = "info", auto.key =  list(columns = 4, legend = TRUE), 
    main = "Test Information Function (TIF)", xlim = c(-3, 3))

# most informative segment
plot(mirt.out, type = "info", auto.key =  list(columns = 4, legend = TRUE), 
    main = "Test Information Function (TIF)", xlim = c(-1, 1))
    
# least informative segment xlim=(3,1)
plot(mirt.out, type = "info", auto.key =  list(columns = 4, legend = TRUE), 
    main = "Test Information Function (TIF)", xlim = c(3, 1))


# Plot person-item map ----------------------------------------------------------

# retrieving coef()$items difficutly estimates
item.difficulties <- coef(mirt.out, simplify = TRUE)$items[, 2]
item.difficulties <- data.frame(Item = paste("Item", 1:length(item.difficulties)), Difficulty = item.difficulties)

# retrieving f-scores (person ability estimates)
person.abilities <- fscores(mirt.out)
person.abilities <- data.frame(PersonID = 1:length(person.abilities), Ability = person.abilities)

item.difficulties.long <- item.difficulties[rep(seq_len(nrow(item.difficulties)), each = nrow(person.abilities)), ]
item.difficulties.long <- item.difficulties[rep(1:nrow(item.difficulties), each = length(person.abilities)), ]

# cols= personID, Item, Difficulty, F1-score
plot_data <- expand.grid(Item = item.difficulties$Item, PersonID = person.abilities$PersonID)
plot_data <- merge(plot_data, item.difficulties, by = "Item")
plot_data <- merge(plot_data, person.abilities, by = "PersonID")

ggplot(plot_data, aes(x = F1, y = Difficulty, color = Item)) +
  geom_point() +
  labs(title = "Person-Item Map", x = "Person Ability", y = "Item Difficulty") +
  theme_minimal() +
  scale_color_brewer(palette = "Set1")

# Plot factor loadings ----------------------------------------------------------

# retrieves items and parameters

# ?? Might need adjustment

loadings <- coef(mirt.out, simplify = TRUE)$items
colnames(loadings) <- c("a1", "d", "c", "u")
loadings <- data.frame(loadings)
loadings.long <- loadings %>%
    tibble::rownames_to_column(var = "Item") %>%
    tidyr::pivot_longer(cols = -Item, names_to = "Parameter", values_to = "Value")


ggplot(loadings.long, aes(x = Item, y = Value, fill = Parameter)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = "Factor Loadings", x = "Item Parameters", y = "Factor Loadings") +
    theme_minimal() +
    scale_fill_brewer(palette = "Set1")


# parameters: a1, 

# Plot residual analysis ----------------------------------------------------------

residuals <- residuals(mirt.out)
residuals <- data.frame(residuals)
residuals <- tibble::rownames_to_column(residuals, var = "Item")
# residuals.long <- melt(residuals, variable.name = "Item", value.name = "Residual")
residuals.long <- melt(residuals, id.vars = "Item", value.name = "Residual")

residuals.long.plot <- residuals.long[!is.na(residuals.long$Residual), ]

ggplot(residuals.long.plot, aes(x = Item, y = Residual, color = abs(Residual))) +
    geom_point(alpha = 0.5) +
    geom_segment(aes(xend = Item, yend = 0), alpha = 0.5) +  # Vertical lines for each point
    geom_hline(yintercept = 0, linetype = "dashed") +  # Horizontal dashed line at y=0
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  # Vertical dashed line at x=0
    labs(title = "Residual Analysis", x = "Item", y = "Residual") +
    theme_minimal() +
    scale_color_gradient(low = "blue", high = "red")
