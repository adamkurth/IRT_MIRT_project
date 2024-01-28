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
}

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
}

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
}

simulate_response_data <- function(n){
    # code
}

calculate_descriptive_stats <- function(n){
    # code
}

plot_distributions <- function(n){
    # code
}


main <- function(n){
    # developing the script with theta_skew.r
    load()
    n <- 500
    all_distributions <- generate_skewed_distribitions(n, seed=123, visualize = TRUE)
    # all_distributions <- generate_skewed_distribitions(n, seed=123, visualize = TRUE)

    head(all_distributions)
    # print(all_distributions) # check
    
}

main(500)
