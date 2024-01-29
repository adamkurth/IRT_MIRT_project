load <- function(){
    library(mirt)
    library(ggplot2)
    library(reshape2)
    library(tibble)
    library(ggpubr)
    library(gridExtra)
    library(sn)
    library(dplyr)
    library(moments)
} # end load

generate_skewed_distribitions <- function(n, seed, visualize = FALSE) {
    set.seed(seed) # Set the seed for reproducibility
    
    theta_left <- rsn(n, xi = 0, omega = 1, alpha = -10) # left skewed
    theta_right <- rsn(n, xi = 0, omega = 1, alpha = 10) # right skewed
    theta_far_left <- rsn(n, xi = 0, omega = 1, alpha = -20) # far left skewed
    theta_far_right <- rsn(n, xi = 0, omega = 1, alpha = 20) # far right skewed
    theta_normal <- rnorm(n, mean = 0, sd = 1) # normal

    all_distributions <- data.frame(
        theta = c(theta_left, theta_right, theta_far_left, theta_far_right, theta_normal),
        distribution = rep(c("Left Skewed", "Right Skewed", "Far Left Skewed", "Far Right Skewed", "Normal"), each = n)
    )
    # format the dataframe for 'theta' and 'distribution'.
    format_df <- function(df){
        # ensure that 'theta', and 'distribution' are the correct data types
        if (!all(c("theta", "distribution") %in% names(df))) {
            stop("The dataframe must have columns named 'theta' and 'distribution'")
        } else {
            df$theta <- as.numeric(df$theta)
            df$distribution <- as.factor(df$distribution)
        }
        return(df)
    }
    
    all_distributions <- format_df(all_distributions)

    if (visualize) {
        p <- plot_skewed_distributions(all_distributions)
        print(p) # Explicitly print the plot
    }
    return(all_distributions)
} # end generate_skewed_distribitions

plot_skewed_distributions <- function(all_distributions) {
    theta <- all_distributions$theta
    distribution <- all_distributions$distribution

    # number of unique distributions
    num_distributions <- length(unique(distribution))

    # color palette
    colors <- RColorBrewer::brewer.pal(min(num_distributions, 9), "Set1")

    p <- ggplot(all_distributions, aes(x = theta, fill = distribution)) +
            geom_density(alpha = 0.5) +
            labs(x = "Theta", y = "Density") +
            scale_fill_manual(values = colors) +
            theme_minimal()
    return(p)
} # end plot_skewed_distributions

simulate_response_data <- function(all_distributions, a, b, c, seed = 123) {
    # helper function to simulate response data
    
    # checks for skewed distribution implementation
    if (!("theta" %in% names(all_distributions)) && !("distribution" %in% names(all_distributions))) {
        stop("The dataframe must contain 'theta' and 'distribution' columns")
    }
    
    all_distributions$theta <- as.numeric(all_distributions$theta)
    all_distributions$distribution <- as.factor(all_distributions$distribution)
    
    simulate_response <- function(theta) {
        n.persons <- length(theta)
        n.items <- length(a)
        set.seed(seed)
        
        response.data <- matrix(NA, n.persons, n.items)
        for (i in 1:n.persons) {
            for (j in 1:n.items) {
                # Calculate probability p for each person-item combination
                p <- c[j] + (1 - c[j]) / (1 + exp(-(a[j] * theta[i] - b[j])))
                # Ensure p is a single numeric value
                if (length(p) != 1) {
                    stop("Probability calculation did not return a single value.")
                }
                response.data[i, j] <- ifelse(runif(1) < p, 1, 0)
            }
        }
        response.df <- as.data.frame(response.data)
        colnames(response.df) <- paste0("I", 1:n.items)
        return(response.df)
    }
    response.dataframes <- lapply(split(all_distributions$theta, all_distributions$distribution), simulate_response)
    
    return(response.dataframes)
}

calculate_descriptive_stats <- function(response.dataframes, all.distributions){
    # names of dist strings 'left' 'normal', etc.
    unique.dists <- unique(all.distributions$distribution)

    # use summarise function for all items 
    calculate_stats <- function(df){
        df %>%
            summarise(across(starts_with("I"),
                list(mean = mean, sd = sd, skewness = moments::skewness, 
                    kurtosis = moments::kurtosis, min = min, max = max)))
    }
    
    stats.tables <- lapply(response.dataframes, calculate_stats)
    names(stats.tables) <- names(response.dataframes)

    comb.stats <- do.call("cbind", lapply(stats.tables, t))
    comb.stats <- as.data.frame(comb.stats)
    colnames(comb.stats) <- unique.dists

    for (i in 1:length(response.dataframes)) {
        print(paste("Dimensions for response.dataframes", i, ":", nrow(response.dataframes[[i]]), ncol(response.dataframes[[i]])))
    }
    print(paste("Dimensions for comb.stats:", nrow(comb.stats), ncol(comb.stats)))
    # colnames(comb.stats) <- new.colnames
    return(comb.stats) 

} # end calculate_descriptive_stats


plot_distributions <- function(n){
    # code
}



main <- function(n){
    # developing the script with theta_skew.r
    load()
    n <- 500
    all_distributions <- generate_skewed_distribitions(n, seed=123, visualize = FALSE)
    # all_distributions <- generate_skewed_distribitions(n, seed=123, visualize = TRUE)
    
    # checks
    # head(all_distributions)
    # print(all_distributions)
    # colnames(all_distributions)

    # true parameters
    a <- c(0.5, 1, 1.5, 2)
    b <- c(0, -1, 1, 0)
    c <- c(0.2, 0.15, 0.25, 0.2)
    
    response.dataframes <- simulate_response_data(all_distributions, a, b, c)
    
    # checks
    # head(response.dataframes)
    # print(response.dataframes)
    # names(response.dataframes)

    # type is what type of dataframe we are analyzing

    calculate_descriptive_stats(response.dataframes, all_distributions)

}

main(500)



