# Multidimensional Item Response Theory (MIRT) simulation code (3PL)
# -------------------------------------------------------------
rm(list=ls())
# install.packages("mirt")
library(mirt)
library(ggplot2)
library(reshape2)
library(tibble)
library(ggpubr)
# -------------------------------------------------------------
# Simulate data with Monte Carlo experiments ------------------
theta.exp <- rexp(n = 500, rate = 1) 
theta.gamma <- rgamma(n = 500, shape = 1, rate = 1)
theta.unif <- runif(n = 500, min = -1, max = 1)
theta.norm <- rnorm(n = 500, mean = 0, sd = 1)
theta.beta <- rbeta(n = 500, shape1 = 1, shape2 = 1)
theta.chisq <- rchisq(n = 500, df = 1)
theta.t <- rt(n = 500, df = 1)
theta.f <- rf(n = 500, df1 = 1, df2 = 1)
# -------------------------------------------------------------
a <- c(0.5, 1, 1.5, 2)
b <- c(0, -1, 1, 0)
c <- c(0.2, 0.15, 0.25, 0.2)

n.persons <- length(theta.exp) # holds for all theta
n.items <- length(a)

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
    theta.exp = simulate_response(theta.exp, seed = 123),
    theta.gamma = simulate_response(theta.gamma, seed = 123),
    theta.norm = simulate_response(theta.norm, seed = 123),
    theta.unif = simulate_response(theta.unif, seed = 123),
    theta.beta = simulate_response(theta.beta, seed = 123),
    theta.chisq = simulate_response(theta.chisq, seed = 123),
    theta.t = simulate_response(theta.t, seed = 123),
    theta.f = simulate_response(theta.f, seed = 123)
)

# renaming the items to I1, I2, I3, I4
response.dataframes <- lapply(response.dataframes, function(df) {
  colnames(df) <- paste0("I", 1:n.items)
  return(df)
})
# -------------------------------------------------------------
# 1. see if distribution of theta (ability) makes a difference 

# 2. skew the theta value to simulate a gifted/screening test scenario
# -------------------------------------------------------------

# 1. Test Distribution ------------------------------------------
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

# Plot estimates ----------------------------------------------------------
plot.est.param <- list()
a_true <- c(0.5, 1, 1.5, 2)
b_true <- c(0, -1, 1, 0)
c_true <- c(0.2, 0.15, 0.25, 0.2)

plots.a <- list()
plots.b <- list()
plots.c <- list()

# Loop through each fitted model
for(name in names(mirt.models)) {
    mirt.out <- mirt.models[[name]]
    estimates <- coef(mirt.out, simplify = TRUE, IRTparts = TRUE)$items
    df <- data.frame(Item = paste("Item", seq_len(nrow(estimates))),
                    Discrimination = estimates[, "a1"],
                    Difficulty = estimates[, "d"],
                    Guessing = estimates[, "g"])
    # Extract estimated parameters
    p_a <- ggplot(df, aes(x = Item, y = Discrimination)) +
        geom_bar(stat = "identity", fill = "blue") +
        geom_hline(yintercept = a_true, linetype = "dashed", color = "black") +
        labs(title = paste("Discrimination (a) -", name),
             x = "Item", y = "Discrimination")

    p_b <- ggplot(df, aes(x = Item, y = Difficulty)) +
        geom_bar(stat = "identity", fill = "red") +
        geom_hline(yintercept = b_true, linetype = "dashed", color = "black") +
        labs(title = paste("Difficulty (b) -", name),
             x = "Item", y = "Difficulty")

    p_c <- ggplot(df, aes(x = Item, y = Guessing)) +
        geom_bar(stat = "identity", fill = "green") +
        geom_hline(yintercept = c_true, linetype = "dashed", color = "black") +
        labs(title = paste("Guessing (c) -", name),
             x = "Item", y = "Guessing")

    plots.a[[name]] <- p_a
    plots.b[[name]] <- p_b
    plots.c[[name]] <- p_c
}
plot_a_combined <- ggarrange(plotlist = plots.a, ncol = 2, nrow = ceiling(length(plots.a) / 2))
plot_b_combined <- ggarrange(plotlist = plots.b, ncol = 2, nrow = ceiling(length(plots.b) / 2))
plot_c_combined <- ggarrange(plotlist = plots.c, ncol = 2, nrow = ceiling(length(plots.c) / 2))

print(plot_a_combined)
print(plot_b_combined)
print(plot_c_combined)

# Don't know if this is a helpful visualization

## plot estimates already in the dataframe



# Plot time to convergence ----------------------------------------------------------
conv.times <- sapply(results, function(df) df$time.to.conv)
conv.times <- data.frame(Theta = names(conv.times), TimeToConvergence = conv.times)

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

ggpubr::ggarrange(plots.ll$theta.exp, plots.ll$theta.gamma, plots.ll$theta.norm, plots.ll$theta.unif,
          plots.ll$theta.beta, plots.ll$theta.chisq, plots.ll$theta.t, plots.ll$theta.f,
          ncol = 2, nrow = 4)

print(plots.ll$theta.exp)
print(plots.ll$theta.gamma)
print(plots.ll$theta.norm)
print(plots.ll$theta.unif)
print(plots.ll$theta.beta)
print(plots.ll$theta.chisq)
print(plots.ll$theta.t)

# Plot item characteristic curves (ICC) ----------------------------------------------------------

plots.icc <- list()

for(name in names(mirt.models)){
    mirt.out <- mirt.models[[name]]

    p <- plot(mirt.out, type = "trace", auto.key =  list(columns = 4, legend = TRUE), 
            main = paste("Item Characteristic Curves (ICC) for", name))
    plots.icc[[name]] <- p
}

print(plots.icc$theta.exp)
print(plots.icc$theta.gamma)
print(plots.icc$theta.norm)
print(plots.icc$theta.unif)
print(plots.icc$theta.beta)
print(plots.icc$theta.chisq)
print(plots.icc$theta.t)

# Plot test information fucntion (TIF) ----------------------------------------------------------

plots.tif <- list()

for(name in names(mirt.models)){
    mirt.out <- mirt.models[[name]]

    p <- plot(mirt.out, type = "info", auto.key =  list(columns = 4, legend = TRUE), 
            main = "Test Information Function (TIF)")

    plots.tif[[name]] <- p
}

print(plots.tif$theta.exp)
print(plots.tif$theta.gamma)
print(plots.tif$theta.norm)
print(plots.tif$theta.unif)
print(plots.tif$theta.beta)
print(plots.tif$theta.chisq)
print(plots.tif$theta.t)

# Plot person-item map ----------------------------------------------------------

plots.person.item <- list()

for(name in names(mirt.models)){
    mirt.out <- mirt.models[[name]]

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

    p <- ggplot(plot_data, aes(x = PersonID, y = Difficulty, color = Difficulty)) +
        geom_point() +
        geom_line(aes(group = Item)) +
        labs(title = paste("Person-Item Map for", name),
             x = "Person ID", y = "Item Difficulty") +
        theme_minimal()

    plots.person.item[[name]] <- p
}

print(plots.person.item$theta.exp)
print(plots.person.item$theta.gamma)
print(plots.person.item$theta.norm)
print(plots.person.item$theta.unif)
print(plots.person.item$theta.beta)
print(plots.person.item$theta.chisq)
print(plots.person.item$theta.t)


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
factor.loadings.df$Item <- factor(factor_loadings_df$Item)

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
