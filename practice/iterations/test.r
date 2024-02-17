rm(list=ls())

load <- function(){
    packages <- c("mirt", "ggplot2", "reshape2", "tibble", "ggpubr", "gridExtra", "sn", "tidyr", "dplyr", "moments", "parallel", "progressr")

    sapply(packages, require, character.only = TRUE)
} # end load

# DONE & working
generate.skewed.distribitions <- function(n, seed=123) {
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

    return(all_distributions)
} # end generate_skewed_distribitions


# DONE & working
simulate.response.data <- function(all_distributions, cross.param, seed = 123) {
    require(mirt, quietly = TRUE)
    set.seed(seed)

    # ensure all dist. have necessary columns: checks for skewed distribution implementation
    if (!("theta" %in% names(all_distributions)) || !("distribution" %in% names(all_distributions))) {
        stop("The dataframe must contain 'theta' and 'distribution' columns")
    }
    
    all_distributions$theta <- as.numeric(all_distributions$theta)
    all_distributions$distribution <- as.factor(all_distributions$distribution)
    
    simulate_response <- function(theta, cross.param) {
        n.persons <- length(theta)
        n.items <- nrow(cross.param)
        a <- cross.param$a
        b <- cross.param$b

        theta.matrix <- matrix(rep(theta, each = n.items), nrow = n.persons, byrow = TRUE)
        a.matrix <- matrix(rep(a, n.persons), nrow = n.persons, byrow = FALSE)
        b.matrix <- matrix(rep(b, n.persons), nrow = n.persons, byrow = FALSE)

        p.matrix <- 1 / (1 + exp(-(a.matrix * (theta.matrix - b.matrix))))

        response.matrix <- ifelse(runif(n.persons * n.items) < p.matrix, 1, 0)

        response.df <- as.data.frame(response.matrix)
        colnames(response.df) <- paste0("I", seq_len(n.items))
        return(response.df)
    }
    # apply simulate_response to each group of theta values
    response.dataframes <- list()
    distributions <- unique(all.distributions$distribution)
    for (distribution in distributions) {
        theta.values <- all_distributions$theta[all_distributions$distribution == distribution]
        response.dataframes[[distribution]] <- simulate_response(theta.values, cross.param)
    }
    return(response.dataframes)
} # end simulate.response.data

fit.mirt <- function(response.data, cross.param, methods, dentypes, n.replications = 100) {
    require(mirt, quietly = TRUE)

    if(!is.data.frame(response.data)) {
        stop("response.data must be a dataframe.")
    }

    # Assuming response.data is correctly formatted and passed to the function
    n.items <- ncol(response.data)  # Number of items
    # n.persons <- nrow(response.data)  # Number of persons

    # Initialize an empty list to store results from all methods and dentypes
    results <- list()
    

    for (method in methods) {
        for (dentype in dentypes) {
            message(paste("Running method:", method, "with dentype:", dentype))
            
            est.params.a <- matrix(NA, nrow = n.replications, ncol = n.items)
            est.params.b <- matrix(NA, nrow = n.replications, ncol = n.items)
            total.time <- numeric(n.replications)
            
            successful_run <- TRUE
        
            for (i in 1:n.replications) {
                cat("\n")
                message(paste("Replication", i, ": Starting"))
                start <- Sys.time()
                
                # Try to fit the model with the current method and dentype
                fit <- tryCatch({
                mirt(data = response.data, model = 1, itemtype = '2PL', storeEMhistory = TRUE, method = method, dentype = dentype)
                }, error = function(e) {
                message(paste("Error in replication", i, ":", e$message))
                    successful_run <- FALSE
                    return(NULL)
                })
                
                if (is.null(fit)) break  # Exit the replication loop if the model fit was unsuccessful
                
                end <- Sys.time()
                total.time[i] <- as.numeric(end - start)
                message(paste("Replication", i, ": Completed", "- Time to convergence:", total.time[i]))
                
                coefs <- coef(fit, simplify = TRUE, IRTpars = TRUE)$items
                est.params.a[i,] <- coefs[, "a"]
                est.params.b[i,] <- coefs[, "b"]
            }
        
            if (successful_run) {
                    avg.time <- mean(total.time, na.rm = TRUE)
                    avg.a <- rowMeans(est.params.a, na.rm = TRUE)
                    avg.b <- rowMeans(est.params.b, na.rm = TRUE)
                    
                    metrics <- data.frame(item = 1:n.items, 
                                        true.a = cross.param$a, true.b = cross.param$b, 
                                        avg.a = avg.a, avg.b = avg.b, avg.time = avg.time)
                    
                    metrics$bias.a <- metrics$avg.a - metrics$true.a
                    metrics$bias.b <- metrics$avg.b - metrics$true.b
                    metrics$rmse.a <- sqrt(mean((metrics$avg.a - metrics$true.a)^2))
                    metrics$rmse.b <- sqrt(mean((metrics$avg.b - metrics$true.b)^2))
                    
                    # Store metrics in the results list with a unique name
                    results[[paste("method", method, "dentype", dentype)]] <- metrics
            } else {
                message(paste("Method:", method, "with dentype:", dentype, "failed. Skipping."))
            }
        }
    }
  
  return(results)
}

fit.mirt.parallel <- function(all.distributions, cross.param, methods, dentypes) {
    require(parallel, quietly = TRUE)
    require(mirt, quietly = TRUE)
    
    # Generate or load response dataframes for each distribution type
    response.dataframes <- simulate.response.data(all.distributions, cross.param)  # Placeholder; adjust as needed

    log <- "calc_averages_parallel.log"
    write(paste("Starting parallel processing at", Sys.time()), file = log)

    numCores <- detectCores() - 1
    cl <- makeCluster(numCores)
    on.exit(stopCluster(cl), add = TRUE)
    
    dist.types <- unique(all.distributions[,2])  # Assuming this correctly identifies unique distribution types

    # Prepare the environment for the cluster
    clusterExport(cl, c("response.dataframes", "fit.mirt", "cross.param", "methods", "dentypes", "dist.types"))
    clusterEvalQ(cl, {
        library(mirt)
        # Ensure simulate.response.data and any other necessary functions are defined or loaded here
    })

    # Function to process each distribution type by its index
    processCombination <- function(i) {
        distType <- dist.types[i]  # Directly use 'i' to access the current distType
        response.data <- response.dataframes[[i]]  # Access response data by index directly
        tryCatch({
            metrics.df <- fit.mirt(response.data, cross.param, methods, dentypes, 100)
            return(metrics.df)
        }, error = function(e) {
            message(paste("Error with distribution type", distType, ":", e$message))
            return(NULL)
        })
    }    
    # Process each distribution type in parallel and store results
    metrics.list <- parLapply(cl, seq_along(dist.types), processCombination)
    
    # Combine results into a single dataframe, handling potential NULLs
    metrics.df <- do.call(rbind, lapply(metrics.list, function(x) if (is.null(x)) data.frame() else x))

    return(metrics.df)
}


load()
n <- 300
all.distributions <- generate.skewed.distribitions(n, seed=123)

# cros.param has 20 rows and 3 columns (of true param values for a, b, and d)
cross.param <- expand.grid(d = c(-2.5, -1.25, 0, 1.25, 2.5), a = c(0.5, 1, 1.5, 2.5))
cross.param$b <- with(cross.param, -d/a)

# simulate response data
# for each dist, has 300 rows of 20 item responses
response.dataframes <- simulate.response.data(all.distributions, cross.param, seed = 123)

methods <- c("BL")
dentypes <- c("Gaussian")
dist.types <- c("stnd.norm")


# Corrected function call to specifically use 'stnd.norm' dataframe
# metrics <- fit.mirt(response.dataframes$stnd.norm, cross.param, methods, dentypes, 100) ## WORKING! 

metrics <- fit.mirt.parallel(all.distributions, cross.param, methods, dentypes) 
