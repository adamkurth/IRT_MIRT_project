load <- function(){
    rm(list=ls())
    library(mirt)
    library(ggplot2)
    library(reshape2)
    library(tibble)
    library(ggpubr)
    library(gridExtra)
    library(sn)
    library(tidyr)
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

fit_mirt_models <- function(response.dataframes, true_params, estimation, dentype){
    results <- list()
    mirt.models <- list()

    for(name in names(response.dataframes)){
        #extract response data
        response.data <- response.dataframes[[name]]

        # fit the 2PL model
        start.time <- Sys.time()
        mirt.model <- mirt(data = response.data, model = 1, itemtype = '2PL', storeEMhistory = TRUE, method=estimation, dentype=dentype)
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

compare.mirt <- function(response.dataframes, true.params, methods, dentypes, dist.types){
    N <- length(true.params$a)
    true.params.matrix <- matrix(unlist(true.params), nrow=N, byrow=TRUE)

    results <- list()
    
    for(method in methods){
        for(distribution in dist.types){
            for(dentype in dentypes){

                response.data <- response.dataframes[[distribution]]

                if (is.null(response.data)) {
                    print(paste("No response data for distribution", distribution, "in method", method))
                    next
                }

                tryCatch({
                    # fit the 2PL model
                    start.time <- Sys.time()
                    mirt.model <- mirt(data = response.data, model = 1, itemtype = '2PL', 
                                        storeEMhistory = TRUE, method=method, dentype=dentype)
                    end.time <- Sys.time()

                    # extract estimated parameters & rename
                    # converged <- mirt.model@converged
                    param.estimates <- coef(mirt.model, simplify = TRUE)$items
                    relevant.params <- param.estimates[, c("a1", "d")] # drop the "g/c" and "u/d" parameters
                    colnames(relevant.params) <- c("a1", "b1")

                    # assign numeric values
                    conv.time <- end.time - start.time
                    # rmse
                    rmse.a <- sqrt(mean((relevant.params[, "a1"] - true.params$a)^2))
                    rmse.b <- sqrt(mean((relevant.params[, "b1"] - true.params$b)^2))
                    # bias 
                    bias.a <- mean(relevant.params[, "a1"] - true.params.matrix[, 1])
                    bias.b <- mean(relevant.params[, "b1"] - true.params.matrix[, 2])

                    # store results
                    results[[paste(method, dentype, distribution, sep=".")]] <- list(
                        convergence.time = conv.time,
                        term.conv = NA,  # FIX
                        rmse.a = rmse.a,
                        rmse.b = rmse.b,
                        bias.a = bias.a,
                        bias.b = bias.b,
                        method = method,
                        distribution = distribution,
                        dentype = dentype
                    )

                }, error = function(e) {
                    message(paste("Error with method", method, "and dentype", dentype, 
                                  "in distribution", distribution, ":", e$message))
                    })
                }
            }
        }

    reorganize <- function(results) {
        organized.results <- do.call(rbind, lapply(names(results), function(key) {
            result <- results[[key]]
            cbind(method = result$method,
                  dentype = result$dentype,
                  distribution = result$distribution,
                  convergence.time = as.numeric(result$convergence.time),
                  rmse.a = result$rmse.a,
                  rmse.b = result$rmse.b,
                  bias.a = result$bias.a,
                  bias.b = result$bias.b)
        }))
        
        return(organized.results)
    }
    
    organized.results <- reorganize(results)
    return(organized.results)
} # end compare.methods


# MAIN SCRIPT
load()
n <- 500
all.distributions <- generate_skewed_distribitions(n, seed=123, visualize = FALSE)

# true parameters: 2PL model
true.params <- data.frame(
    a = c(0.5, 1, 1.5, 2), 
    b = c(0, -1, 1, 0)
) 

response.dataframes <- simulate_response_data(all.distributions, true.params, seed = 123)

true.params <- list(a = true.params$a, b = true.params$b) # to list

# ESTIMATION METHODS
# 1. Block & Lieberman approach (BL)
# 2. Expectation- Maximization (EM) algorithm
# 3. Metropolis-Hastings Robbins-Monro (MHRM)
# 4. Monte Carlo EM (MCEM)
# 5. Stochastic EM
# 6. Quasi-Monte Carlo EM

dist.types <- c("left.skew", "right.skew", "stnd.norm")
method.types <- c("BL", "EM", "MHRM", "MCEM", "SEM", "QMCEM")
dentype.types <- c("Gaussian", "EH", "EHW","Davidian-2", "Davidian-4")


# model.results <- fit_mirt_models(response.dataframes, true_params, esimation=est.methods, dentype=default)
# results <- model.results$Results

# response.dataframes, true.params, methods, dentypes, dist.types

mirt.results <- compare.mirt(response.dataframes, true.params, 
                methods = method.types, dentypes = dentype.types, dist.types = dist.types)
View(mirt.results)






#should we incorperate the all.distributions? in these models results? That is the primary interest here. 

# Q: Should we incorperate the all.distributions when running the fit_mirt_models function?
# A: No, the all.distributions dataframe is not used in the fit_mirt_models function. The function only uses the response dataframes and the true parameters. The all.distributions dataframe is only used to generate the response dataframes.

# DIRECTIONS: 
# EXAMINE THE SKEWED DISTRIBUTIONS:
# 1. Skew right
# 2. Skew left
# 3. Standard normal

# EXAMINE ESTIMATION MEHTODS:
# 1. Block & Lieberman approach
# 2. Expectation- Maximization (EM) algorithm
# 3. Metropolis-Hastings Robbins-Monro (MHRM)
# 4. Monte Carlo EM (MCEM)
# 5. Stochastic EM
# 6. Quasi-Monte Carlo EM

# OUTCOME METRICS:
# 1. RMSE (for each item, averaged over 100 replications)
# 2. Bias (for each item, averaged over 100 replications)
# 3. Coveragece time



