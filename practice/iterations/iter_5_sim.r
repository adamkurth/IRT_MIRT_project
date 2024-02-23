rm(list=ls())

load <- function(){
    packages <- c("mirt", "ggplot2", "reshape2", "tibble", "ggpubr", "gridExtra", "sn", "tidyr", "dplyr", "moments", "parallel", "progressr")

    sapply(packages, require, character.only = TRUE)
} # end load

# Method 1
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

# Method 2
### revised functions 
quick.gen.dist <- function(n, dist.type = 'stnd.norm', seed=123){
    require(sn, quietly = TRUE)
    set.seed(seed) # Set the seed for reproducibility   
    data <- data.frame(theta = numeric(n), distribution = character(n))
    
    if (dist.type == 'stnd.norm') {
        data$theta <- rnorm(n, mean = 0, sd = 1)  # Standard normal
        data$distribution <- 'stnd.norm'
    } else if (dist.type == 'left.skew') {
        data$theta <- rsn(n, xi = 0, omega = 1, alpha = -10)  # Left skewed
        data$distribution <- 'left.skew'
    } else if (dist.type == 'right.skew') {
        data$theta <- rsn(n, xi = 0, omega = 1, alpha = 10)  # Right skewed
        data$distribution <- 'right.skew'
    } else {
        stop("dist.type must be 'stnd.norm', 'left.skew', or 'right.skew'")
    }
    
    return(data)
}  # end quick.gen.dist

quick.sim.response <- function(theta.values, cross.param, seed=123) {
    set.seed(seed)
    n.persons <- length(theta.values) # 300 persons
    n.items <- nrow(cross.param) 
    a <- cross.param$a
    b <- cross.param$b

    # create matrices for theta, a, and b
    theta.matrix <- matrix(rep(theta.values, each = n.items), nrow = n.persons, byrow = TRUE)
    a.matrix <- matrix(rep(a, n.persons), nrow = n.persons, byrow = FALSE)
    b.matrix <- matrix(rep(b, n.persons), nrow = n.persons, byrow = FALSE)
    
    # calculate the probability matrix
    p.matrix <- 1 / (1 + exp(-(a.matrix * (theta.matrix - b.matrix))))
    response.matrix <- ifelse(runif(n.persons * n.items) < p.matrix, 1, 0)

    response.df <- as.data.frame(response.matrix)
    colnames(response.df) <- paste0("I", seq_len(n.items))
    return(response.df)
} # end quick.sim.response

# for reading a specific response data file 
read.data <- function(dist.type='stnd.norm', rep=1, read.dir='response_data'){
    dist.type.pattern <- gsub("\\.", "_", dist.type)
    file.pattern <- sprintf("response_%s_rep_%02d.txt", dist.type.pattern, rep)
    
    file.path <- file.path(read.dir, file.pattern)
    if (!file.exists(file.path)) {
        stop(paste("File not found:", file.path))
    }
    response.data <- read.table(file.path, header = TRUE, sep = "\t")
    print(paste("Read response data from file:", file.path))
    return(response.data)
} # end read.data

export.data <- function(cross.param, output.dir='response_data', n=300, replications=10, seed=123){
    dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)
    dist.types <- c('stnd.norm', 'left.skew', 'right.skew')

    # for each distribution, simulate response data and write to .txt file
    for (dist.type in dist.types){
        for (rep in 1:replications){
            # generate theta values
            theta.data <- quick.gen.dist(n, dist.type = dist.type, seed = seed + rep) # seed is incremented by rep
            # simulate response data
            response.data <- quick.sim.response(theta.data$theta, cross.param, seed=123+rep) # seed is incremented by rep
            # generate filename with replication number 
            filename <- sprintf("%s/response_%s_rep_%02d.txt", output.dir, gsub("\\.", "_", dist.type), rep)
            write.table(response.data, filename, sep = "\t", row.names = FALSE, col.names = TRUE)
            message("Wrote response data to file: ", filename)
        }
    }

} # end export.data

#### Method 2
fit.mirt <- function(dist.type, rep, cross.param, methods, dentypes, n.replications=10) {
    require(mirt, quietly = TRUE)
    read.dir <- 'response_data'

    # initialize results list
    results <- list()
    n.items <- length(cross.param$a)

    for (method in methods) {
        for (dentype in dentypes) {
            message(paste("Running method:", method, "with dentype:", dentype))
            
            est.params.a <- data.frame(matrix(NA, nrow = n.replications, ncol = n.items))
            est.params.b <- data.frame(matrix(NA, nrow = n.replications, ncol = n.items))

            total.time <- numeric(n.replications)
            successful.run <- TRUE
        
            for (rep in 1:n.replications) {
                cat("\n")
                message(sprintf("Processing: Method = %s, Replication = %d", method, rep))
                
                # read response data for current replication
                response.data <- read.data(dist.type = dist.type, rep = rep, read.dir = read.dir)

                message(paste("Replication", rep, ": Starting"))
                start.time <- Sys.time()
                
                tryCatch({
                    # model fit
                    fit <- mirt(model = 1, data = response.data, itemtype="2PL", method=method, dentype=dentype, verbose = FALSE)
                    end.time <- Sys.time()
        
                    total.time[rep] <- as.numeric(end.time - start.time)
                    message(paste("Replication", rep, ": Completed", "- Time to convergence:", total.time[rep]))

                    # coefficient 
                    coefs <- coef(fit, simplify = TRUE)$items # working
                    est.params.a[rep, ] <- coefs[, "a1"] # for a  
                    est.params.b[rep, ] <- coefs[, "d"] # for b

                }, error = function(e) {
                    message(sprintf("Error during replication %d: %s", rep, e$message))
                    successful.run <- FALSE
                })                
                if (is.null(fit)) {
                    message(paste("Replication", rep, ": Failed. Skipping."))
                    next # skip if failed to fit
                }
                
            } # end for rep loop
        
            if (successful.run) {
                metrics <- data.frame(item = 1:n.items, true.a = cross.param$a, true.b = cross.param$b)
                metrics$avg.time <- mean(total.time, na.rm = TRUE)

                # bias calculation for each item
                metrics$bias.a <- colMeans(est.params.a, na.rm = TRUE) - metrics$true.a
                metrics$bias.b <- colMeans(est.params.b, na.rm = TRUE) - metrics$true.b

                # RMSE calculation for each item
                metrics$rmse.a <- sqrt(colMeans((est.params.a - rep(metrics$true.a, each = n.replications))^2, na.rm = TRUE))
                metrics$rmse.b <- sqrt(colMeans((est.params.b - rep(metrics$true.b, each = n.replications))^2, na.rm = TRUE))

                # Store metrics in the results list with a unique name
                results[[paste("method", method, "dentype", dentype)]] <- metrics
                
                # calculating the average bias and RMSE for each item
                # summarizing over replications for each item not accross all items!
            } else {
                message(paste("Method:", method, "with dentype:", dentype, "failed. Skipping."))
            }
        }
    }
  
  return(results)
} # end fit.mirt

fit.mirt.parallel <- function(all.distributions, cross.param, methods, dentypes) {
    require(parallel, quietly = TRUE)
    require(mirt, quietly = TRUE)
    
    # Generate or load response dataframes for each distribution type
    response.dataframes <- simulate.response.data(all.distributions, cross.param)  # Placeholder; adjust as needed

    log <- "parallel_fit_mirt.log"
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
} # end fit.mirt.parallel

read.data.all <- function(){
    # read in response data from .txt files in the "response_data" directory
    response.dataframes <- list()
    filenames <- list.files(path = "response_data", pattern = "response_.*\\.txt", full.names = TRUE)
    
    for (f in filenames){
        response.data <- read.table(f, header = TRUE, sep = "\t")

        name <- basename(f) # remove path, leaving only the filename
        name <- gsub("^response_", "", name)  # Remove 'response_' prefix
        name <- gsub(".txt", "", name)  # Remove '.txt' suffix
        name <- gsub("_", ".", name)  # Replace '_' with '.'
        # modified name as the key name
        response.dataframes[[name]] <- response.data
    }
    return(response.dataframes)
} # end read.data.all



load()
n <- 300
######
# Method 1: initial implementation
all.distributions <- generate.skewed.distribitions(n, seed=123)

# cros.param has 20 rows and 3 columns (of true param values for a, b, and d)
cross.param <- expand.grid(d = c(-2.5, -1.25, 0, 1.25, 2.5), a = c(0.5, 1, 1.5, 2.5))
cross.param$b <- with(cross.param, -d/a)

# simulate response data: for each dist, has 300 rows of 20 item responses
response.dataframes <- simulate.response.data(all.distributions, cross.param, seed = 123)
# output here: response.dataframes$left.skew, response.dataframes$right.skew, response.dataframes$stnd.norm

######
# Method 2: revised implementation 
dist.data <- quick.gen.dist(300, dist.type='stnd.norm', seed=123)
response.data <- quick.sim.response(dist.data$theta, cross.param, seed=123) # test call 
export.data(cross.param, output.dir='response_data', n=300, replications=10, seed=123)

methods <- c("BL")
dentypes <- c("Gaussian")
dist.types <- c("stnd.norm")


# Corrected function call to specifically use 'stnd.norm' dataframe
metrics <- fit.mirt(dist.type='stnd.norm', rep=1, cross.param, methods, dentypes, 10) # for quick testing

# metrics.revised <- fit.mirt(den.type='stnd.norm', rep, cross.param, methods, dentypes, 100) 
# metrics.revised # view

# fit.mirt <- function(dist.type, rep, cross.param, methods, dentypes, n.replications = 10) {
