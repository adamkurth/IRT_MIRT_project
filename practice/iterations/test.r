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


calc.metrics <- function(estimated.params, cross.param){
  # Assuming `estimated.params` contains averaged estimates across replications for each item
  
  # True parameters
  true.a <- cross.param$a
  true.b <- cross.param$b
  
  # Averaged estimated parameters
  avg.est.a <- mean(estimated.params$a)
  avg.est.b <- mean(estimated.params$b)
  
  # Calculating bias for averaged estimates
  bias.a <- avg.est.a - true.a
  bias.b <- avg.est.b - true.b
  
  # Calculating RMSE for averaged estimates
  rmse.a <- sqrt(mean((estimated.params$a - true.a)^2))
  rmse.b <- sqrt(mean((estimated.params$b - true.b)^2))
  
  # Assuming `time` is an average over replications if relevant
  # time <- mean(estimated.params$conv.time)
  
  # Compile into a single-row dataframe
  metrics.df <- data.frame(
    true.a = true.a, 
    true.b = true.b, 
    avg.est.a = avg.est.a, 
    avg.est.b = avg.est.b, 
    bias.a = bias.a, 
    bias.b = bias.b, 
    rmse.a = rmse.a, 
    rmse.b = rmse.b
    # ,time = time  # Uncomment if time is relevant and available
  )
  
  return(metrics.df)
}


fit.mirt.models.all <- function(response.data, cross.param, method, dentype, dist.types) {
    require(mirt, quietly = TRUE)
    n.replications <- 100
    all.est.params <- vector("list", n.replications)
    all.times <- numeric(n.replications)
    
    

    for (i in 1:n.replications){
        cat("\n")
        
        message(paste("Replication", i, "-", method, dentype, ": Starting"))
        
        current.data <- response.data[[i]]
        
        # Diagnostic check for current.data
        if(is.null(current.data) || length(current.data) == 0) {
            message(paste("Replication", i, "-", method, dentype, ": current.data is empty or NULL. Skipping."))
            next  # Skip this iteration
        }

        current.data <- response.data[[i]]
        start.time <- Sys.time()
        # Fit model
        mirt.model <- tryCatch({
            mirt(data = current.data, model = 1, itemtype = '2PL', storeEMhistory = TRUE, method = method, dentype = dentype)
        }, error = function(e) {
            message(paste("Replication", i, "-", method, dentype, ": Error -", e$message))
            return(NULL)  # Return NULL if error
        })
        end.time <- Sys.time()
        time.to.conv <- as.numeric(end.time - start.time)
        all.times[i] <- time.to.conv

        if (!is.null(mirt.model)) {
            para.estimates <- coef(mirt.model, simplify = TRUE, IRTpars = TRUE)$items
            relevant.param.estimates <- para.estimates[, c("a", "b")]
            all.est.params[[i]] <- list(a = relevant.param.estimates[, "a"], b = relevant.param.estimates[, "b"])
        } else {
            all.est.params[[i]] <- list(a = NA, b = NA)
        }
        message(paste("Replication", i, "-", method, dentype, ": Completed"))
    }

    avg.times <- mean(all.times, na.rm = TRUE)
    avg.a <- mean(unlist(lapply(all.est.params, function(x) x$a), na.rm = TRUE))
    avg.b <- mean(unlist(lapply(all.est.params, function(x) x$b), na.rm = TRUE))
    
    aggregated.results <- list(estimated.params = list(a = avg.a, b = avg.b), avg.conv.time = avg.times)
    return(aggregated.results)
}

calc.averages.parallel <- function(all.distributions, cross.param, methods, dentypes, dist.types){
    require(parallel, quietly = TRUE)
    require(mirt, quietly = TRUE)

    log <- "calc_averages_parallel.log"
    write(paste("Starting parallel processing at", Sys.time()), log)

    numCores <- detectCores() - 1
    cl <- makeCluster(numCores)
    on.exit({
        stopCluster(cl)
        writeLines(paste("Computation complete.", Sys.time()), log)
    }, add = TRUE)

    # environment prep for cluster
    clusterExport(cl, list("simulate.response.data", "fit.mirt.model", "calc.metrics", "all.distributions", "cross.param", "methods", "dentypes", "dist.types"))
    clusterEvalQ(cl, library(mirt))

    # generate combinations of parameters to process
    combinations <- expand.grid(method=methods, dentype=dentypes, dist.type=dist.types)
    n.combinations <- nrow(combinations)

    # function to process each combination
    processCombination <- function(i){
        comb <- combinations[i, ]

        tryCatch({
            simData <- simulate.response.data(all.distributions, cross.param)
            responseData <- simData[[comb$dist.type]]
            fit <- fit.mirt.model(responseData, cross.param, comb$method, comb$dentype)
            metrics.df <- calc.metrics(fit$estimated.params, cross.param)
            metrics.df$method <- comb$method
            metrics.df$dist.type <- comb$dist.type
            metrics.df$dentype <- comb$dentype
            return(metrics.df)
        }, error = function(e) {
            message(paste("Error with combination", i, ":", e$message))
            return(NULL)
        })
    }
    
    # process each combination in parallel and store results 
    metrics.list <- parLapply(cl, seq_len(n.combinations), processCombination) # outputs a list of dataframes 
    
    # combine list of dataframes into one dataframe
    metrics.df <- do.call(rbind, Filter(NROW, metrics.list)) # remove NULL elements and combine into 1 dataframe
    
    return(metrics.df)
} # end calc.averages.parallel 


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

# averages.data <- calc.averages.parallel(all.distributions, cross.param, methods = methods, dentypes = dentypes, dist.types = dist.types)
# View(averages.data)


combinations <- expand.grid(method=methods, dentype=dentypes, dist.type=dist.types)
fit <- fit.mirt.models.all(response.dataframes$stnd.norm, cross.param, "BL", "Gaussian")

# fit.mirt.model <- function(response.data, cross.param, method, dentype, n.replications) {

test_combination <- function() {
    
    i <- 1 # Assuming 1 is a valid index for a test
    comb <- combinations[i, ]
    dist <- as.character(comb$dist.type)
    method <- as.character(comb$method)
    dentype <- as.character(comb$dentype)

    simData <- simulate.response.data(all.distributions, cross.param)
    responseData <- simData[[dist]]
    fit <- fit.mirt.models.all(responseData, cross.param, method, dentype)
    metrics.df <- calc.metrics(fit$estimated.params, cross.param)
    print(metrics.df)
    metrics.df$method <- method 
    metrics.df$dist.type <- dist
    metrics.df$dentype <- dentype
    return(metrics.df)
}

test.df <- test_combination()
View(test.df)














fit.mirt <- function(response.dataframes, cross.param) {
    require(mirt, quietly = TRUE)
    response.data <- response.dataframes$stnd.norm  # Only for standard normal distribution

    n.replications <- 100
    n.items <- ncol(response.data)  # Number of items
    total.time <- numeric(n.replications)  # Store time for each replication
    
    # Initialize matrices to store parameter estimates across replications
    est.params.a <- matrix(NA, nrow = n.replications, ncol = n.items)
    est.params.b <- matrix(NA, nrow = n.replications, ncol = n.items)

    for (i in 1:n.replications) {
        cat("\n")
        message(paste("Replication", i, ": Starting"))
        start <- Sys.time()
        
        fit <- mirt(data = response.data, model = 1, itemtype = '2PL', storeEMhistory = TRUE)
        
        end <- Sys.time()
        total.time[i] <- as.numeric(end - start)  # Store time for this replication
        message(paste("Replication", i, ": Completed", "- Time to convergence:", total.time[i]))
        
        coefs <- coef(fit, simplify = TRUE, IRTpars = TRUE)$items
        est.params.a[i,] <- coefs[, "a"]
        est.params.b[i,] <- coefs[, "b"]
    }
    
    # Average parameter estimates and total time across replications
    avg.a <- rowMeans(est.params.a, na.rm = TRUE)
    avg.b <- rowMeans(est.params.b, na.rm = TRUE)
    avg.time <- mean(total.time, na.rm = TRUE)
    
    # Calculate metrics based on averaged parameter estimates
    metrics <- data.frame(item = 1:n.items, 
                          true.a = cross.param$a, true.b = cross.param$b, 
                          avg.a = avg.a, avg.b = avg.b, avg.time = avg.time)
    
    metrics$bias.a <- metrics$avg.a - metrics$true.a
    metrics$bias.b <- metrics$avg.b - metrics$true.b
    metrics$rmse.a <- sqrt(mean((metrics$avg.a - metrics$true.a)^2))
    metrics$rmse.b <- sqrt(mean((metrics$avg.b - metrics$true.b)^2))

    return(metrics)
}

metrics <- fit.mirt(response.dataframes, cross.param)

