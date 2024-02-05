load <- function(){
    library(mirt)
    library(ggplot2)
    library(reshape2)
    library(tibble)
    library(ggpubr)
    library(gridExtra)
    library(sn)
    library(tidyr)å
    library(dplyr)
    library(moments)
} # end load

generate_skewed_distribitions <- function(n, seed, visualize = FALSE) {
    require(sn, quietly = TRUE)
    set.seed(seed) # Set the seed for reproducibility
    
    theta_left <- rsn(n, xi = 0, omega = 1, alpha = -10) # left skewed
    theta_right <- rsn(n, xi = 0, omega = 1, alpha = 10) # right skewed
    theta_normal <- rnorm(n, mean = 0, sd = 1) # standard normal

    all_distributions <- data.frame(
        theta = c(theta_left, theta_right, theta_normal),
        distribution = rep(c("left.skew", "right.skew", "stnd.norm"), each = n)
    )

    # Format the dataframe for 'theta' and 'distribution'
    all_distributions$theta <- as.numeric(all_distributions$theta)
    all_distributions$distribution <- as.factor(all_distributions$distribution)
    # check formatting
    print(paste("Dimensions for all_distributions dataframe:", nrow(all_distributions), ncol(all_distributions)))
    print(head(all_distributions))
    print(tail(all_distributions))  

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

simulate_response_data <- function(all_distributions, params, seed = 123) {
    require(mirt, quietly = TRUE)
    # set seed for reproducibility
    set.seed(seed)

    # helper function to simulate response data
    a <- params$a
    b <- params$b

    # Structure: rows = items, columns = parameters (a, b)
    param.matrix <- as.matrix(params)

    # ensure all dist. have necessary columns: checks for skewed distribution implementation
    if (!("theta" %in% names(all_distributions)) || !("distribution" %in% names(all_distributions))) {
        stop("The dataframe must contain 'theta' and 'distribution' columns")
    }
    
    all_distributions$theta <- as.numeric(all_distributions$theta)
    all_distributions$distribution <- as.factor(all_distributions$distribution)
    
    simulate_response <- function(theta, params) {
        n.persons <- length(theta)
        n.items <- nrow(params)
        set.seed(seed)

        response.data <- matrix(NA, n.persons, n.items)
        for (i in 1:n.persons) {
            for (j in 1:n.items) {
                # Extract item params for the jth item
                a_j <- params[j, "a"]
                b_j <- params[j, "b"]

                # Calculate probability p for each person-item combination using the 2PL model
                p <- 1 / (1 + exp(-(a_j * (theta[i] - b_j))))
                response.data[i, j] <- ifelse(runif(1) < p, 1, 0)
            }
        }
        response.df <- as.data.frame(response.data)
        colnames(response.df) <- paste0("I", 1:n.items)
        return(response.df)
    }

    # apply simulate_response to each group of theta values
    response.dataframes <- lapply(split(all_distributions$theta, all_distributions$distribution), simulate_response, params = param.matrix)
    
    return(response.dataframes)
}

fit_mirt_models <- function(response.dataframes, all.distributions, true_params){
    results <- list()
    mirt.models <- list()

    for(name in names(response.dataframes)){
        #extract response data
        response.data <- response.dataframes[[name]]

        # fit the 2PL model
        start.time <- Sys.time()
        mirt.model <- mirt(data = response.data, model = 1, itemtype = '2PL', storeEMhistory = TRUE)
        end.time <- Sys.time()
        time.to.conv <- end.time - start.time

        # extract estimated parameters
        param.estimates <- coef(mirt.model, simplify = TRUE)$items
        colnames(param.estimates) <- c("a1", "b1", "c", "d")
        print(param.estimates)

        # calculate RMSE
        rmse.a <- sqrt(mean((param.estimates[, "a1"] - true_params$a)^2))
        rmse.b <- sqrt(mean((param.estimates[, "b1"] - true_params$b)^2))
        
        # store results
        df <- data.frame(
            distribution = name,
            a.discrim = true_params$a,
            b.diff = true_params$b,
            a.est = param.estimates[, "a1"], 
            b.est = param.estimates[, "b1"],
            rmse.a = rmse.a,
            rmse.b = rmse.b,
            time.to.conv = as.numeric(time.to.conv),
            stringsAsFactors = FALSE
        )

        results[[name]] <- df
        mirt.models[[name]] <- mirt.model 
    }   

    return(list(Results = results, MirtModels = mirt.models))
} # end mirt.models


calculate_descriptive_stats_as_tibble <- function(response.dataframes, all.distributions){
    # names of dist strings 'left' 'normal', etc.
    # extracting names of all.distributions dataframe
    df.list <- lapply(names(response.dataframes), function(name){
        df <- response.dataframes[[name]]
        response.vector <- as.vector(t(df))  # Corrected to use as.vector
        temp.df <- data.frame(
            distribution = rep(name, times = length(response.vector) / ncol(df)), 
            item = rep(colnames(df), each = nrow(df)), 
            response = response.vector
        ) 
        return(temp.df)
    })

    # combine all temp dataframes into one
    df.long <- do.call(rbind, df.list)

    # ensure df.long treated as dataframe
    df.long <- as.data.frame(df.long)

    # set distribution column as a factor
    df.long$distribution <- factor(df.long$distribution, levels = unique(all.distributions$distribution))
    
    stats <- df.long %>%
        group_by(item, distribution) %>%
        summarise(
            mean = mean(response, na.rm = TRUE),
            sd = sd(response, na.rm = TRUE),
            skewness = moments::skewness(response, na.rm = TRUE),
            kurtosis = moments::kurtosis(response, na.rm = TRUE),
            min = min(response, na.rm = TRUE),
            max = max(response, na.rm = TRUE), 
            .groups = 'drop'
        ) %>%
        pivot_wider(names_from = distribution, values_from = c(mean, sd, skewness, kurtosis, min, max))
    
    #transpose 
    stats <- t(stats)
    
    print(paste("Dimensions for stats dataframe:", nrow(stats), ncol(stats)))

  return(stats)
}# end calculate_descriptive_stats

calculate_descriptive_stats_base <- function(response.dataframes, all.distributions){
    # Initialize an empty list to store transformed data frames
    df.list <- lapply(names(response.dataframes), function(name){
        df <- response.dataframes[[name]]
        response.vector <- as.vector(t(df))
        # Create a temporary data frame for each distribution
        temp.df <- data.frame(
            distribution = rep(name, times = length(response.vector) / ncol(df)), 
            item = rep(colnames(df), each = nrow(df)), 
            response = response.vector
        ) 
        return(temp.df)
    })

    # Combine all temporary data frames into one
    df.long <- do.call(rbind, df.list)

    # Convert to a data frame and set the distribution column as a factor
    df.long <- as.data.frame(df.long)
    df.long$distribution <- factor(df.long$distribution, levels = unique(all.distributions$distribution))

    # Calculate descriptive statistics
    stats <- df.long %>%
        group_by(item, distribution) %>%
        summarise(
            mean = mean(response, na.rm = TRUE),
            sd = sd(response, na.rm = TRUE),
            skewness = moments::skewness(response, na.rm = TRUE),
            kurtosis = moments::kurtosis(response, na.rm = TRUE),
            min = min(response, na.rm = TRUE),
            max = max(response, na.rm = TRUE),
            .groups = 'drop'
        ) %>%
        pivot_wider(names_from = distribution, values_from = c(mean, sd, skewness, kurtosis, min, max))

    # Convert the result to a standard data frame to avoid tibble output
    stats <- as.data.frame(stats)
    # Optional: Print dimensions of the stats data frame
    print(paste("Dimensions for stats dataframe:", nrow(stats), ncol(stats)))
    
    return(stats)
}



main <- function(n){
    # developing the script with theta_skew.r
    load()
    n <- 500
    all_distributions <- generate_skewed_distribitions(n, seed=123, visualize = FALSE)
    # all_distributions <- generate_skewed_distribitions(n, seed=123, visualize = TRUE)
    
    # checks
    # print("Checking all_distributions-------------------")
    # # print("Head:")
    # print(head(all_distributions))
    # print("Print:")
    # print(all_distributions)
    # print(paste("colnames:", paste(colnames(all_distributions), collapse = ", ")))
    # print(paste("Dimensions:", paste(dim(all_distributions), collapse = ", ")))

    # true parameters: 2PL model
    true_params <- data.frame(
        a = c(0.5, 1, 1.5, 2), 
        b = c(0, -1, 1, 0)
    ) 

    
    # check
    print("Parameter dataframe:")
    print(true_params)
    
    response.dataframes <- simulate_response_data(all_distributions, true_params, seed = 123)
    
    # checks
    print("Checking response.dataframes-------------------")
    # print("Head:")
    # print(head(response.dataframes))
    print("Print Names:")
    print(names(response.dataframes))


    true_params <- list(a = true_params$a, b = true_params$b) # to list
    model.results <- fit_mirt_models(response.dataframes, all_distributions, true_params)

    # checks
    print("Checking model.results-------------------")
    print("Print Names:")
    print(names(model.results))
    print("Print Results:")
    print(model.results$Results)
    print("Print MirtModels:")
    print(model.results$MirtModels)

#should we incorperate the all.distributions? in these models results? That is the primary interest here. 



    # type is what type of dataframe we are analyzing

    # might need to fix

    # descriptive.stats <- calculate_descriptive_stats_as_tibble(response.dataframes, all_distributions) 
    # descriptive.stats <- t(descriptive.stats)
    # descriptive.stats <- calculate_descriptive_stats_base(response.dataframes, all_distributions) 

    # print("Checking descriptive.stats-------------------")
    # print("Head:")
    # print(head(descriptive.stats))
    # print("Print:")
    # print(descriptive.stats)
    # View(descriptive.stats)
    # print(paste("colnames:", paste(colnames(descriptive.stats), collapse = ", ")))
    # print(paste("Dimensions:", paste(dim(descriptive.stats), collapse = ", ")))


}

main(500)



