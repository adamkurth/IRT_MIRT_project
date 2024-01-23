# Multidimensional Item Response Theory (MIRT) simulation code (3PL)
# -------------------------------------------------------------
rm(list=ls())
# install.packages("mirt")
library(mirt)
library(ggplot2)
library(reshape2)
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
        p <- c[j] + (1 - c[j]) / (1 + exp(-(a[j] * theta[i] - b[j]))) #2PL model
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

# very sharp peak around theta = 0

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

# Plot factor loadings ----------------------------------------------------------

# retrieves items and parameters
loadings <- coef(mirt.out, simplify = TRUE)$items
loadings.long <- reshape2::melt(loadings)
colnames(loadings.long) <- c("Item", "Parameter", "Value")

ggplot(loadings.long, aes(x = Item, y = Value, fill = Parameter)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = "Factor Loadings", x = "Item Parameters", y = "Factor Loadings") +
    theme_minimal() +
    scale_fill_brewer(palette = "Set1")

# Plot residual analysis ----------------------------------------------------------

residuals <- residuals(mirt.out)
residuals <- data.frame(residuals)
residuals_long <- reshape2::melt(residuals)
colnames(residuals_long) <- c("Var1", "Residual")  

ggplot(residuals_long, aes(x = Var1, y = Residual, color = abs(Residual))) +
    geom_point(alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(title = "Residual Analysis", x = "Item", y = "Residual") +
    theme_minimal() +
    scale_color_gradient(low = "blue", high = "red")

