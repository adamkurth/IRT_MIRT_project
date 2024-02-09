load <- function(){
    packages <- c("mirt", "ggplot2", "reshape2", "tibble", "ggpubr", "gridExtra", "sn", "tidyr", "dplyr", "moments", "parallel", "progressr")
    sapply(packages, require, character.only = TRUE)
} # end load

generate.skewed.distribitions <- function(n, seed, visualize = FALSE) {
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

# revise 
compare.mirt <- function(response.dataframes, true.params, methods, dentypes, dist.types){
    # N <- length(true.params$a)
    # true.params.matrix <- matrix(unlist(true.params), nrow=N, byrow=TRUE)

    results <- list()
    
    for(method in methods){
        for(distribution in dist.types){
            for(dentype in dentypes){
                response.data <- response.dataframes[[distribution]]

                if (is.null(response.data)) {
                    print(paste("No response data for distribution", distribution, "in method", method))
                    next
                } # end if null


                tryCatch({
                    # fit the 2PL model
                    start.time <- Sys.time()
                    mirt.model <- mirt(data = response.data, model = 1, itemtype = '2PL', 
                                        storeEMhistory = TRUE, method=method, dentype=dentype)
                    end.time <- Sys.time()

                    # extract estimated parameters & rename
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
                }) # end tryCatch
            } # end dentype
        }   # end distribution
    } # end method

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
        })) # end do.call
        
        return(organized.results)
    } # end reorganize

    return(reorganize(results))
} # end compare.mirt


# called by calc.averages
fit.mirt.model <- function(response.data, cross.param, method, dentype) {
    message(paste("Starting: " , method, dentype, sep=" "))
    start.time <- Sys.time()
    mirt.model <- tryCatch({
        mirt(data = response.data, model = 1, itemtype = '2PL', storeEMhistory = TRUE, method=method, dentype=dentype)
    }, error = function(e) {
        message(paste("Error with method", method, "and dentype", dentype, ":", e$message))
        return(NULL) # return NULL if error
    })
    end.time <- Sys.time()
    time.to.conv <- end.time - start.time

    if(is.null(mirt.model)){
        return(list(estimated.params = NULL, conv.time = NULL))
    }

    # extract estimated parameters & rename
    param.estimates <- coef(mirt.model, simplify = TRUE)$items
    relevant.param.estimates <- param.estimates[, c("a1", "d")] # drop the "g/c" and "u/d" parameters
    colnames(relevant.param.estimates) <- c("a1", "b1")
    
    estimated.params <- list(a = relevant.param.estimates[, "a1"], b = relevant.param.estimates[, "b1"])
    
    message(paste("Completed:", method, dentype, sep=" "))
    out <- list(estimated.params = estimated.params, conv.time = time.to.conv)
    return(out)
}

# called by calc.averages
calc.metrics <- function(estimated.params, cross.param){
    # estimated.params$a <- as.numeric(estimated.params$a)
    # estimated.params$b <- as.numeric(estimated.params$b)
    rmse.a <- sqrt(mean((estimated.params$a - cross.param$a)^2))
    rmse.b <- sqrt(mean((estimated.params$b - cross.param$b)^2))
    bias.a <- mean(estimated.params$a - cross.param$a)
    bias.b <- mean(estimated.params$b - cross.param$b)
    time <- estimated.params$conv.time
    out <- list(rmse.a = rmse.a, rmse.b = rmse.b, bias.a = bias.a, bias.b = bias.b, time = time)
    return(out)
}

# master function
calc.averages <- function(all.distributions, cross.param, methods, dentypes, dist.types){
    aggregated.results <- list()

    for(method in methods){
        for(dentype in dentypes){
            for(distribution in dist.types){
                metrics.list <- list(rmse.a = numeric(0),
                                     rmse.b = numeric(0),
                                     bias.a = numeric(0),
                                     bias.b = numeric(0),
                                     conv.time = numeric(0)) 
                
                for (rep in 1:100){
                    print(paste("Method:", method, "Dentype:", dentype, "Distribution:", distribution, "Rep:", rep))

                    set.seed(rep)
                    sim.data <- simulate.response.data(all.distributions, cross.param, seed = rep)
                    response.data <- sim.data[[distribution]]

                    fit.results <- fit.mirt.model(response.data, cross.param, method, dentype)
                    metrics <- calc.metrics(fit.results$estimated.params, cross.param)

                    metrics.list$rmse.a <- c(metrics.list$rmse.a, metrics$rmse.a)
                    metrics.list$rmse.b <- c(metrics.list$rmse.b, metrics$rmse.b)
                    metrics.list$bias.a <- c(metrics.list$bias.a, metrics$bias.a)
                    metrics.list$bias.b <- c(metrics.list$bias.b, metrics$bias.b)
                    metrics.list$conv.time <- c(metrics.list$conv.time, fit.results$time.to.conv)
                }

                aggregated.results[[paste(method, dentype, distribution, sep=".")]] <- list(
                    method = method,
                    dentype = dentype,
                    distribution = distribution,
                    rmse.a = mean(metrics.list$rmse.a),
                    rmse.b = mean(metrics.list$rmse.b),
                    bias.a = mean(metrics.list$bias.a),
                    bias.b = mean(metrics.list$bias.b),
                    conv.time = mean(metrics.list$time)
                )
                
            }
        }
    }

    return(aggregated.results)
} # end fit.mirt

calc.averages.parallel <- function(all.distributions, cross.param, methods, dentypes, dist.types) {
    require(parallel, quietly = TRUE)
    require(mirt, quietly = TRUE)

    # Initialize logging
    logFile <- "calc_avg_parallel.log"
    write(paste("Computation complete.", Sys.time()), file = logFile, append = TRUE)

    # Setup parallel cluster
    numCores <- detectCores() - 1
    cl <- makeCluster(numCores)
    on.exit({
        stopCluster(cl)
        writeLines(paste("Computation complete.", Sys.time()), logFile)
    }, add = TRUE)

    # Prepare environment in each cluster node
    clusterExport(cl, list("simulate.response.data", "fit.mirt.model", "calc.metrics", "all.distributions", "cross.param", "methods", "dentypes", "dist.types"))
    clusterEvalQ(cl, library(mirt))

    # Generate combinations of parameters to process
    combinations <- expand.grid(method = methods, dentype = dentypes, dist = dist.types, stringsAsFactors = FALSE)
    totalCombinations <- nrow(combinations)

    # Process combinations in batches
    processCombination <- function(index) {
        comb <- combinations[index, ]
        tryCatch({
            simData <- simulate.response.data(all.distributions, cross.param, seed = index)
            responseData <- simData[[comb$dist]]
            fitResults <- fit.mirt.model(responseData, cross.param, comb$method, comb$dentype)
            metrics <- calc.metrics(fitResults$estimated.params, cross.param)
            result <- c(comb, metrics)
            return(data.frame(t(result)))  # Ensure the return is a data frame
        }, error = function(e) {
            message <- sprintf("Error in %s %s %s: %s", comb$method, comb$dentype, comb$dist, e$message)
            return(data.frame(error = message))  # Return a data frame with error message
        })
    }

    # Execute the processing in parallel
    results <- parLapply(cl, seq_len(totalCombinations), processCombination)

    # Filter out results with errors for logging
    errors <- Filter(function(x) "error" %in% names(x), results)
    if (length(errors) > 0) {
        errorMessages <- sapply(errors, function(e) e$error)
        writeLines(paste(errorMessages, logFile, append = TRUE))
    }

    # Combine valid results into a single data frame
    validResults <- do.call(rbind, Filter(NROW, results))
    if (nrow(validResults) == 0) return(NULL)  # Return NULL if no valid results
    return(validResults)
} # end 



# MAIN SCRIPT
load()
n <- 500
all.distributions <- generate.skewed.distribitions(n, seed=123, visualize = FALSE)

# true parameters: 2PL model
true.params <- data.frame(
    a = c(0.5, 1, 1.5, 2), 
    b = c(0, -1, 1, 0)
) 

# 4x4 grid of true parameters (dimensions: 16 row, 2 columns = a, b)
cross.param <- expand.grid(a = true.params$a, b = true.params$b)

response.dataframes <- simulate.response.data(all.distributions, cross.param, seed = 123)

methods <- c("BL", "EM", "MHRM", "MCEM", "SEM", "QMCEM")
dentypes <- c("Gaussian", "EH", "EHW","Davidian-2", "Davidian-4")
dist.types <- c("left.skew", "right.skew", "stnd.norm")

# averages.data <- calc.averages(all.distributions, cross.param, methods = methods, dentypes = dentypes, dist.types = dist.types)
# View(averages.data)

averages.data <- calc.averages.parallel(all.distributions, cross.param, methods = methods, dentypes = dentypes, dist.types = dist.types)
View(averages.data)
