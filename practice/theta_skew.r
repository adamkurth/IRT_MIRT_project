# -------------------------------------------------------------
# Multidimensional Item Response Theory (MIRT) Simulation (3PL)
# -------------------------------------------------------------
rm(list=ls())
# install.packages("mirt")
library(mirt)
library(ggplot2)
library(reshape2)
library(tibble)
library(ggpubr)
library(gridExtra)
library(sn) # for skew normal distribution
# -------------------------------------------------------------
# Generate skewed theta distributions --------------------------\
n <- 500 
linspace <- seq(-5, 5, length.out = n)

# skew normal distribution
theta.left <- rsn(n, xi = 0, omega = 1, alpha = -10) # left skewed
theta.right <- rsn(n, xi = 0, omega = 1, alpha = 10) # right skewed
theta.far.left <- rsn(n, xi = 0, omega = 1, alpha = -20) # far left skewed
theta.far.right <- rsn(n, xi = 0, omega = 1, alpha = 20) # far right skewed
theta.normal <- rnorm(n, mean = 0, sd = 1) # normal

# Visualize each distribution ----------------------------------
all_distributions <- data.frame(
    theta = c(theta.left, theta.right, theta.far.left, theta.far.right, theta.normal),
    distribution = rep(c("Left Skewed", "Right Skewed", "Far Left Skewed", "Far Right Skewed", "Normal"), each = n)
)

ggplot(all_distributions, aes(x = theta, fill = distribution)) +
    geom_density(alpha = 0.5) +
    labs(x = "Theta", y = "Density") +
    scale_fill_manual(values = c("red", "blue", "green", "purple", "black")) +
    theme_minimal()


# -------------------------------------------------------------
# True Item Parameters -----------------------------------------
a <- c(0.5, 1, 1.5, 2)
b <- c(0, -1, 1, 0)
c <- c(0.2, 0.15, 0.25, 0.2)

# n.persons <- length(a)
n.items <- length(b)
n.persons <- length(theta.normal)

# function to simulate response data using 3PL model
simulate_response <- function(theta, seed = 123){
    set.seed(seed) # for all thetas
    response.data <- matrix(NA, n.persons, n.items)
    for (i in 1:n.persons) {
        for (j in 1:n.items) {
            p <- c[j] + (1 - c[j]) / (1 + exp(-(a[j] * theta[i] - b[j])))
            response.data[i, j] <- ifelse(runif(1) < p, 1, 0)
        }
    }
    return(data.frame(response.data))
} 

# store all response dataframes in a list
response.dataframes <- list(
    left = simulate_response(theta.left, seed = 123),
    right = simulate_response(theta.right, seed = 123),
    far.left = simulate_response(theta.far.left, seed = 123),
    far.right = simulate_response(theta.far.right, seed = 123),
    normal = simulate_response(theta.normal, seed = 123)
)

# renaming the items to I1, I2, I3, I4 
response.dataframes <- lapply(response.dataframes, function(df) {
  colnames(df) <- paste0("I", 1:n.items)
  return(df)
})

# What we have so far for theta dists:

# > head(response.dataframes$theta.norm)
#   I1 I2 I3 I4
# 1  1  0  0  0
# 2  0  1  0  0
# 3  0  0  0  0
# 4  0  1  1  0
# 5  1  1  1  0
# 6  0  0  0  0

# -------------------------------------------------------------
# Goal of this script: 
# see if skewness of the distribution of theta (ability) makes a difference 
# -------------------------------------------------------------

# Fit the models ------------------------------------------------
results <- list()
mirt.models <- list()

for(name in names(response.dataframes)){
    response.data <- response.dataframes[[name]]

    # Fit 3PL model
    start.time <- Sys.time()
    mirt.out <- mirt(data = response.data, model = 1, itemtype = "3PL", storeEMhistory=TRUE)
    end.time <- Sys.time()
    time.to.conv <- end.time - start.time

    # extract estimated parameters
    estimates <- coef(mirt.out, simplify = TRUE, IRTparts = TRUE)$items

    # calc RMSE
    RMSE.a <- sqrt(mean((a - estimates[, 1])^2))
    RMSE.b <- sqrt(mean((b - estimates[, 2])^2))
    RMSE.c <- sqrt(mean((c - estimates[, 3])^2))

    # store results and model
    df <- data.frame(a.discrim = a, b.difficulty = b, c.guessing = c, 
                     a.est = estimates[, 1], b.est = estimates[, 2], c.est = estimates[, 3],
                     RMSE.a = RMSE.a, RMSE.b = RMSE.b, RMSE.c = RMSE.c,
                     time.to.conv = as.numeric(time.to.conv))
    results[[name]] <- df
    mirt.models[[name]] <- mirt.out
}

# -------------------------------------------------------------
true_a <- c(0.5, 1, 1.5, 2)
true_b <- c(0, -1, 1, 0)
true_c <- c(0.2, 0.15, 0.25, 0.2)

plots.a <- list()
plots.b <- list()
plots.c <- list()

for(name in names(results)){
    estimates <- results[[name]][, c("a.est", "b.est", "c.est")]
    df <- data.frame(Item = paste("Item", seq_len(nrow(estimates))),
                     TrueDiscrimination = true_a,
                     TrueDifficulty = true_b,
                     TrueGuessing = true_c,
                     Discrimination = estimates[, "a.est"],
                     Difficulty = estimates[, "b.est"],
                     Guessing = estimates[, "c.est"])
    # for plotting 
    df_long <- reshape2::melt(df, id.vars = "Item", 
                            measure.vars = c("TrueDiscrimination", "Discrimination", 
                                            "TrueDifficulty", "Difficulty", 
                                            "TrueGuessing", "Guessing"))
                                                
    # Plot for Discrimination
    p_a <- ggplot(df_long[df_long$variable %in% c("TrueDiscrimination", "Discrimination"), ], 
                  aes(x = Item, y = value, fill = variable)) +
           geom_bar(stat = "identity", position = position_dodge()) +
           scale_fill_manual(values = c("TrueDiscrimination" = "darkblue", "Discrimination" = "skyblue")) +
           labs(title = paste("Discrimination (a) -", name), x = "Item", y = "Value")

    # Plot for Difficulty
    p_b <- ggplot(df_long[df_long$variable %in% c("TrueDifficulty", "Difficulty"), ], 
                  aes(x = Item, y = value, fill = variable)) +
           geom_bar(stat = "identity", position = position_dodge()) +
           scale_fill_manual(values = c("TrueDifficulty" = "darkred", "Difficulty" = "lightcoral")) +
           labs(title = paste("Difficulty (b) -", name), x = "Item", y = "Value")

    # Plot for Guessing
    p_c <- ggplot(df_long[df_long$variable %in% c("TrueGuessing", "Guessing"), ], 
                  aes(x = Item, y = value, fill = variable)) +
           geom_bar(stat = "identity", position = position_dodge()) +
           scale_fill_manual(values = c("TrueGuessing" = "darkgreen", "Guessing" = "lightgreen")) +
           labs(title = paste("Guessing (c) -", name), x = "Item", y = "Value")

    plots.a[[name]] <- p_a
    plots.b[[name]] <- p_b
    plots.c[[name]] <- p_c
}

all_plots_a <- ggarrange(plotlist = plots.a, common.legend = TRUE, legend = "bottom")
all_plots_b <- ggarrange(plotlist = plots.b, common.legend = TRUE, legend = "bottom")
all_plots_c <- ggarrange(plotlist = plots.c, common.legend = TRUE, legend = "bottom")

print(all_plots_a)
print(all_plots_b)
print(all_plots_c)
# -------------------------------------------------------------

# Plot time to convergence ----------------------------------------------------------
conv.times <- data.frame(
    Theta = names(results),
    TimeToConvergence = sapply(results, function(df) df$time.to.conv[1])
    # NumIterations = sapply(mirt.models, function(model) length(extract.mirt(model, 'EMhistory')$deviance))
)

# Plot time to convergence ----------------------------------------------------------

p <- ggplot(conv.times, aes(x = Theta, y = TimeToConvergence)) +
    geom_col(fill = "blue") +
    labs(title = "Time to Convergence for Each Theta Distribution",
        x = "Theta Distribution", y = "Time to Convergence (seconds)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p)

# Plot log-likelihood ----------------------------------------------------------

plots.ll <- list()
for(name in names(mirt.models)){
    ll.history <- extract.mirt(mirt.models[[name]], 'LLhistory')
    df <- data.frame(Iteration = 1:length(ll.history), LogLikelihood = ll.history)
    
    p <- ggplot(df, aes(x = Iteration, y = LogLikelihood)) +
        geom_line(color = 'blue') +
        geom_point(color = 'red') +
        labs(title = paste("Log-Likelihood Convergence for", name),
             x = "Iteration", y = "Log-Likelihood") +
        theme_minimal()
    
    plots.ll[[name]] <- p
}

combined.plots <- ggarrange(plotlist = lapply(plots.ll, function(p) p), 
                            ncol = 2, nrow = ceiling(length(plots.ll) / 2))
print(combined.plots)

print(plots.ll$left)
print(plots.ll$right)
print(plots.ll$far.left)
print(plots.ll$far.right)
print(plots.ll$normal)

# Plot item characteristic curves (ICC) ----------------------------------------------------------

plots.icc <- list()
gen.plot <- function(name, model){
    p <- plot(model, type = "trace", auto.key =  list(columns = 4, legend = TRUE), 
            main = paste("Item Characteristic Curves (ICC) for", name))
    return(p)
}

for(name in names(mirt.models)){
    mirt.out <- mirt.models[[name]]

    p <- gen.plot(name, mirt.out)

    plots.icc[[name]] <- p
}

grid.arrange(grobs = plots.icc, ncol = 2)

print(plots.icc$left)
print(plots.icc$right)
print(plots.icc$far.left)
print(plots.icc$far.right)
print(plots.icc$normal)


# Plot test information fucntion (TIF) ----------------------------------------------------------

plots.tif <- list()

for(name in names(mirt.models)){
    mirt.out <- mirt.models[[name]]

    p <- plot(mirt.out, type = "info", auto.key =  list(columns = 4, legend = TRUE), 
            main = "Test Information Function (TIF)")

    plots.tif[[name]] <- p
}

print(plots.tif$left)
print(plots.tif$right)
print(plots.tif$far.left)
print(plots.tif$far.right)
print(plots.tif$normal)

# Plot person-item map ----------------------------------------------------------

plots.person.item <- list()

for (name in names(mirt.models)) {
    mirt.out <- mirt.models[[name]]

    # retrieving coef()$items difficulty estimates
    item.difficulties <- coef(mirt.out, simplify = TRUE)$items[, 2]
    item.difficulties <- data.frame(Item = paste("Item", seq_along(item.difficulties)), Difficulty = item.difficulties)

    # retrieving f-scores (person ability estimates)
    person.abilities <- fscores(mirt.out)
    person.abilities <- data.frame(PersonID = seq_along(person.abilities), Ability = person.abilities)

    item.difficulties.long <- item.difficulties[rep(seq_len(nrow(item.difficulties)), each = nrow(person.abilities)), ]

    # cols = personID, Item, Difficulty, F1-score
    plot_data <- expand.grid(Item = item.difficulties$Item, PersonID = person.abilities$PersonID)
    plot_data <- merge(plot_data, item.difficulties, by = "Item")

    p <- ggplot(plot_data, aes(x = PersonID, y = Difficulty, color = Difficulty)) +
        geom_point() +
        geom_line(aes(group = Item)) +
        labs(title = paste("Person-Item Map for", name),
                 x = "Person ID", y = "Item Difficulty") +
        theme_minimal()

    plots.person.item[[name]] <- p
}

print(plots.person.item$left)
print(plots.person.item$right)
print(plots.person.item$far.left)
print(plots.person.item$far.right)
print(plots.person.item$normal)

# Plot factor loadings ----------------------------------------------------------

factor.loadings.df <- data.frame(Theta = character(),
                                 Item = character(),
                                 Discrimination = numeric(),
                                 stringsAsFactors = FALSE)
                            
for(theta in names(results)){
    a.est <- results[[theta]]$a.est

    temp.df <- data.frame(Theta = rep(theta, length(a.est)),
                        Item = paste("Item", seq_along(a.est)), 
                        Discrimination = a.est)

    factor.loadings.df <- rbind(factor.loadings.df, temp.df)
}

factor.loadings.df$Theta <- factor(factor.loadings.df$Theta, levels = names(mirt.models))
factor.loadings.df$Item <- factor(factor.loadings.df$Item)

p <- ggplot(factor.loadings.df, aes(x = Item, y = Discrimination, fill = Theta)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    labs(title = "Discrimination Parameter Estimates Across Theta Distributions",
         x = "Item", y = "Discrimination Parameter (a)") +
    scale_fill_brewer(palette = "Set1") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p)


# Plot residual analysis ----------------------------------------------------------

all.residuals <- data.frame(Theta = character(),
                            Item = character(),
                            Residual = numeric(),
                            stringsAsFactors = FALSE)

for(theta in names(mirt.models)){
    mirt.model <- mirt.models[[theta]]
    # extract residuals
    residuals.df <- data.frame(residuals(mirt.model))
    residuals.df <- tibble::rownames_to_column(residuals.df, var = "Item")

    # convert to long format
    residuals.long <- melt(residuals.df, id.vars = "Item", value.name = "Residual")
    
    # add theta dist name
    residuals.long$Theta <- theta

    # add to all residuals df
    all.residuals <- rbind(all.residuals, residuals.long)
}
all.residuals <- all.residuals[!is.na(all.residuals$Residual), ]

p <- ggplot(all.residuals, aes(x = Item, y = Residual, color = abs(Residual))) +
    geom_point(alpha = 0.5) +
    geom_segment(aes(xend = Item, yend = 0), alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed") +  
    labs(title = "Residual Analysis Across Theta Distributions", x = "Item", y = "Residual") +
    theme_minimal() +
    scale_color_gradient(low = "blue", high = "red") +
    facet_wrap(~ Theta) +   
    theme(legend.position = "bottom", legend.justification = "right")

print(p)

#  Plot RMSE ----------------------------------------------------------

rmse.plots <- list()
library(tidyverse)
for(theta in names(results)){
    df <- results[[theta]]

    rmse.df <- df %>%
        select(contains("RMSE")) %>%
        pivot_longer(everything(), names_to = "Parameter", values_to = "RMSE") %>%
        mutate(Item = row_number())


    p <- ggplot(rmse.df, aes(x = as.factor(Item), y = RMSE, fill = Parameter)) +
        geom_bar(stat = "identity", position = position_dodge()) +
        labs(title = paste("RMSE for", theta),
             x = "Item", y = "RMSE") +
        scale_fill_brewer(palette = "Set1") +
        theme_minimal()

    rmse.plots[[theta]] <- p
}

combined.plots <- ggarrange(plotlist = lapply(rmse.plots, function(p) p), 
                                ncol = 2, nrow = ceiling(length(rmse.plots) / 2))
print(combined.plots)



# Misc. plots ----------------------------------------------------------
# might be useful later?
rmse.df$Item <- as.factor(rmse.df$Item)
rmse.df$Parameter <- gsub("RMSE.", "", rmse.df$Parameter)  # Clean up parameter names

p_bar <- ggplot(rmse.df, aes(x = Item, y = RMSE, fill = Parameter)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    labs(title = "RMSE by Parameter and Item", x = "Item", y = "RMSE") +
    theme_minimal()

p_line <- ggplot(rmse.df, aes(x = Item, y = RMSE, color = Parameter, group = Parameter)) +
    geom_line() +
    labs(title = "RMSE Trends by Parameter", x = "Item", y = "RMSE") +
    theme_minimal()

p_box <- ggplot(rmse.df, aes(x = Parameter, y = RMSE)) +
    geom_boxplot() +
    labs(title = "Distribution of RMSEs by Parameter", x = "Parameter", y = "RMSE") +
    theme_minimal()

print(p_bar)
print(p_line)
print(p_box)
