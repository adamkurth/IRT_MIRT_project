rm(list=ls())

load <- function(){
    packages <- c("mirt", "ggplot2", "reshape2", "tibble", "ggpubr", "gridExtra", "sn", "tidyr", "dplyr", "moments", "parallel", "dplyr", "stringr", "tidyverse")

    sapply(packages, require, character.only = TRUE)
} # end load

# ----------------------------------------------------
# Simulation Data ------------------------------------
# ----------------------------------------------------

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
    
    simulate.response <- function(theta, cross.param) {
        n.persons <- length(theta)
        n.items <- nrow(cross.param)
        a <- cross.param$a
        b <- cross.param$b

        theta.matrix <- matrix(rep(theta, each = n.items), nrow = n.persons, byrow = TRUE) 
        a.matrix <- matrix(rep(a, n.persons), nrow = n.persons, byrow = TRUE) # BYROW = TRUE
        b.matrix <- matrix(rep(b, n.persons), nrow = n.persons, byrow = TRUE) # BYROW = TRUE  

        p.matrix <- 1 / (1 + exp(-(a.matrix * (theta.matrix - b.matrix))))

        response.matrix <- ifelse(runif(n.persons * n.items) < p.matrix, 1, 0)
        response.df <- as.data.frame(response.matrix)
        colnames(response.df) <- paste0("I", seq_len(n.items))
        
        # added: return probabilities for plotting
        # return probabilities for plotting 
        prob.df <- as.data.frame(p.matrix)
        colnames(prob.df) <- paste0("P", seq_len(n.items))

    return(list(response = response.df, probabilities = prob.df))
    }
    # apply simulate_response to each group of theta values
    response.dataframes <- list()
    prob.dataframes <- list()
    distributions <- unique(all.distributions$distribution)
    for (distribution in distributions) {
        theta.values <- all_distributions$theta[all_distributions$distribution == distribution]
        simulated.data <- simulate.response(theta.values, cross.param)
        response.dataframes[[distribution]] <- simulated.data$response # added
        prob.dataframes[[distribution]] <- simulated.data$probabilities # added
    }
    return(list(response = response.dataframes, probabilities = prob.dataframes))
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
    a.matrix <- matrix(rep(a, n.persons), nrow = n.persons, byrow = TRUE) # BYROW = TRUE adjusted 
    b.matrix <- matrix(rep(b, n.persons), nrow = n.persons, byrow = TRUE) # BYROW = TRUE adjusted 
    
    # calculate the probability matrix
    p.matrix <- 1 / (1 + exp(-(a.matrix * (theta.matrix - b.matrix))))
    response.matrix <- ifelse(runif(n.persons * n.items) < p.matrix, 1, 0)

    response.df <- as.data.frame(response.matrix)
    prob.df <- as.data.frame(p.matrix)
    colnames(prob.df) <- paste0("P", seq_len(n.items))
    colnames(response.df) <- paste0("I", seq_len(n.items))
    return(list(response = response.df, probabilities = prob.df))
} # end quick.sim.response

# ----------------------------------------------------
# Read/Write -----------------------------------------
# ----------------------------------------------------

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
            sim <- quick.sim.response(theta.data$theta, cross.param, seed=123+rep) # seed is incremented by rep
            response.data <- sim$response # only want pure response data
            # generate filename with replication number 
            filename <- sprintf("%s/response_%s_rep_%02d.txt", output.dir, gsub("\\.", "_", dist.type), rep)
            write.table(response.data, filename, sep = "\t", row.names = FALSE, col.names = TRUE)
            message("Wrote response data to file: ", filename)
        }
    }

} # end export.data

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

# ----------------------------------------------------
# fit.mirt functions ---------------------------------
# ----------------------------------------------------

fit.mirt <- function(dist.type,  cross.param, methods, dentypes, n.replications=10) {
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
            extreme.case <- FALSE
            failed.due.to.data.0 <- 0
            failed.due.to.data.1 <- 0
            failed.due.to.error <- 0
            
            for (rep in 1:n.replications) {
                cat("\n")
                message(sprintf("Processing: Method = %s, Replication = %d", method, rep))
                
                # read response data for current replication
                response.data <- read.data(dist.type = dist.type, rep = rep, read.dir = read.dir)
                should.continue <- TRUE

                # Check for extreme scenarios: all 0's or all 1's in any item
                if(any(colSums(response.data) == 0) || any(colSums(response.data) == nrow(response.data))) {
                    if(any(colSums(response.data) == 0)) {
                        message(paste("Replication", rep, ": Skipping due to all 0 responses"))
                        failed.due.to.data.0 <- failed.due.to.data.0 + 1
                    } else {
                        message(paste("Replication", rep, ": Skipping due to all 1 responses"))
                        failed.due.to.data.1 <- failed.due.to.data.1 + 1
                    }
                    next  # Skip this replication
                    # if uncommented, this error occurs and continues the script
                    # "Error during replication 10: The following items have only one response category and cannot be estimated: I4" 
                    should.continue <- FALSE
                }
                # start replication
                message(paste("Replication", rep, ": Starting"))
                
                tryCatch({
                    # model fit
                    start.time <- Sys.time()
                    fit <- mirt(model = 1, data = response.data, itemtype="2PL", method=method, dentype=dentype)
                    end.time <- Sys.time()
        
                    total.time[rep] <- as.numeric(end.time - start.time)
                    message(paste("Replication", rep, ": Completed", "- Time to convergence:", total.time[rep]))

                    # coefficient 
                    coefs <- coef(fit, simplify = TRUE, IRTpars = TRUE )$items

                    est.params.a[rep, ] <- coefs[, "a"]
                    est.params.b[rep, ] <- coefs[, "b"]

                }, error = function(e) {
                    message(sprintf("Error during replication %d: %s", rep, e$message))
                    failed.due.to.error <- failed.due.to.error + 1
                    should.continue <- FALSE  # Update control variable
                    
                })
                if (!should.continue) {
                    next  # Skip to the next replication if an error was caught
                }                

                if (is.null(fit)) {
                    message(paste("Replication", rep, ": Failed. Skipping."))
                    # next # skip if failed to fit
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

                metrics$est.a <- colMeans(est.params.a, na.rm = TRUE)
                metrics$est.b <- colMeans(est.params.b, na.rm = TRUE)

                # Store metrics in the results list with a unique name
                results[[paste0("method_", method, "_dentype_", dentype)]] <- metrics
                message(paste("\nMethod:", method, "with dentype:", dentype, "-- Successful iteration"))
            } else {
                message(paste("Method:", method, "with dentype:", dentype, "had issues. Check warnings."))
            }
            # Print summary of failed replications
            if(extreme.case == TRUE) {
                message(paste('\t', failed.due.to.data.0, "replications failed due to all 0's."))
                message(paste('\t', failed.due.to.data.1, "replications failed due to all 1's."))
                message(paste('\t', failed.due.to.error, "replications failed due to error.\n"))
            }
            
        } # end for dentype loop
    } # end for method loop
  
  return(results)
} # end fit.mirt

fit.mirt.parallel <- function(all.distributions, cross.param, methods, dentypes) {
    require(parallel, quietly = TRUE)
    require(mirt, quietly = TRUE)
    
    # Generate or load response dataframes for each distribution type
    sim <- simulate.response.data(all.distributions, cross.param)  # Placeholder; adjust as needed
    response.dataframes <- sim$response  # added 
    probabilities.dataframes <- sim$probabilities  # added

    # FIX THIS: Need to ensure response dataframes are implemented correctly

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

# ----------------------------------------------------
# Plotting Functions----------------------------------
# ----------------------------------------------------

# Line Charts (Seperately) ---------------------------
# Plotting RMSE for a/b parameters SEPERATELY
plot.rmse.a.diff <- function(metrics_list) {
    # Combine all metrics into a single dataframe with 'method' and 'dentype' columns
    combined_data <- bind_rows(lapply(names(metrics_list), function(name) {
        df <- metrics_list[[name]]
        method_dentype <- strsplit(name, "_", fixed = TRUE)[[1]]
        method <- method_dentype[2]
        dentype <- method_dentype[4]
        df %>% mutate(method = method, dentype = dentype, item = as.factor(item))
    }), .id = "id") %>% 
        select(-id) %>% 
        mutate(method = factor(method), dentype = factor(dentype))

    # Ensure 'rmse.a' is available for plotting
    if (!"rmse.a" %in% names(combined_data)) {
        stop("The dataframe does not contain 'rmse.a' column.")
    }

    # Optionally, sort by 'item' if needed
    combined_data <- combined_data %>% arrange(item)

    # Dynamically calculate the number of columns for a symmetric grid
    n_plots <- length(unique(combined_data$item))
    n_cols <- ceiling(sqrt(n_plots))

    # Plot RMSE differences by Method and Dentype across a flexible grid
    p <- ggplot(combined_data, aes(x = item, y = rmse.a, group = interaction(method, dentype), color = method)) +
        geom_point() +
        geom_line() +
        facet_wrap(~method + dentype, scales = "free_y", ncol = n_cols) +  # Using facet_wrap with dynamic ncol
        labs(title = "RMSE differences for parameter 'a' across methods ",
             x = "Item",
             y = "RMSE for 'a'") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "bottom") +
        scale_color_manual(values = c("BL" = "blue", "EM" = "red", "MCEM" = "green")) # add more colors as needed

    return(p)
}
plot.rmse.b.diff <- function(metrics_list) {
    # Combine all metrics into a single dataframe with 'method' and 'dentype' columns
    combined_data <- bind_rows(lapply(names(metrics_list), function(name) {
        df <- metrics_list[[name]]
        method_dentype <- strsplit(name, "_", fixed = TRUE)[[1]]
        method <- method_dentype[2]
        dentype <- method_dentype[4]
        df %>% mutate(method = method, dentype = dentype, item = as.factor(item))
    }), .id = "id") %>% 
        select(-id) %>% 
        mutate(method = factor(method), dentype = factor(dentype))

    # Ensure 'rmse.b' is available for plotting
    if (!"rmse.b" %in% names(combined_data)) {
        stop("The dataframe does not contain 'rmse.b' column.")
    }

    # Optionally, sort by 'item' if needed
    combined_data <- combined_data %>% arrange(item)

    # Dynamically calculate the number of columns for a symmetric grid
    n_plots <- length(unique(combined_data$item))
    n_cols <- ceiling(sqrt(n_plots))

    # Plot RMSE differences by Method and Dentype across a flexible grid
    p <- ggplot(combined_data, aes(x = item, y = rmse.b, group = interaction(method, dentype), color = method)) +
        geom_point() +
        geom_line() +
        facet_wrap(~method + dentype, scales = "free_y", ncol = n_cols) +  # Using facet_wrap with dynamic ncol
        labs(title = "RMSE differences for parameter 'b' across methods ",
                 x = "Item",
                 y = "RMSE for 'b'") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
                    legend.position = "bottom") +
        scale_color_manual(values = c("BL" = "blue", "EM" = "red", "MCEM" = "green")) # add more colors as needed

    return(p)
}
# Plotting BIAS for a/b parameters SEPERATELY
plot.bias.a.diff <- function(metrics_list) {
    combined_data <- bind_rows(lapply(names(metrics_list), function(name) {
        df <- metrics_list[[name]]
        method_dentype <- strsplit(name, "_", simplify = TRUE)[[1]]
        method <- method_dentype[2]
        dentype <- method_dentype[4]
        df %>% mutate(method = method, dentype = dentype, item = as.factor(item))
    }), .id = "id") %>%
        select(-id) %>%
        mutate(method = factor(method), dentype = factor(dentype))

    if (!"bias.a" %in% names(combined_data)) {
        stop("The dataframe does not contain 'bias.a' column.")
    }

    combined_data <- combined_data %>% arrange(item)

    n_plots <- length(unique(combined_data$item))
    n_cols <- ceiling(sqrt(n_plots))

    p <- ggplot(combined_data, aes(x = item, y = bias.a, group = interaction(method, dentype), color = method)) +
        geom_point() +
        geom_line() +
        facet_wrap(~method + dentype, scales = "free_y", ncol = n_cols) +
        labs(title = "Bias differences for parameter 'a' across methods ",
                 x = "Item",
                 y = "Bias for 'a'") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
                    legend.position = "bottom") +
        scale_color_manual(values = c("BL" = "blue", "EM" = "red", "MCEM" = "green")) # add more colors as needed

    return(p)
}
plot.bias.b.diff <- function(metrics_list) {
    combined_data <- bind_rows(lapply(names(metrics_list), function(name) {
        df <- metrics_list[[name]]
        method_dentype <- strsplit(name, "_", simplify = TRUE)[[1]]
        method <- method_dentype[2]
        dentype <- method_dentype[4]
        df %>% mutate(method = method, dentype = dentype, item = as.factor(item))
    }), .id = "id") %>%
        select(-id) %>%
        mutate(method = factor(method), dentype = factor(dentype))

    if (!"bias.b" %in% names(combined_data)) {
        stop("The dataframe does not contain 'bias.b' column.")
    }

    combined_data <- combined_data %>% arrange(item)

    n_plots <- length(unique(combined_data$item))
    n_cols <- ceiling(sqrt(n_plots))

    p <- ggplot(combined_data, aes(x = item, y = bias.b, group = interaction(method, dentype), color = method)) +
        geom_point() +
        geom_line() +
        facet_wrap(~method + dentype, scales = "free_y", ncol = n_cols) +
        labs(title = "Bias differences for parameter 'b' across methods ",
                 x = "Item",
                 y = "Bias for 'b'") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
                    legend.position = "bottom") +
        scale_color_manual(values = c("BL" = "blue", "EM" = "red", "MCEM" = "green")) # add more colors as needed

    return(p)
}
# Line Charts (Together) ------------------------------
# Designed to work with entire list of dataframes (metrics_list)
extract.info.from.list <- function(metrics.list) {
    # Use lapply to iterate over all items in the list
    info_list <- lapply(names(metrics.list), function(name) {
        # Split the name to extract method and dentype
        name_parts <- strsplit(name, "_")[[1]]
        
        # Use robust indexing to handle variable naming conventions
        method_index <- which(name_parts == "method")
        dentype_index <- which(name_parts == "dentype")
        
        # Extract method and dentype assuming they follow 'method' and 'dentype' labels
        # This assumes that 'method' and 'dentype' are followed by their respective values
        method <- if (length(method_index) > 0 && (method_index + 1) <= length(name_parts)) {
            name_parts[method_index + 1]
        } else {
            NA  # Return NA if the pattern doesn't match
        }
        
        dentype <- if (length(dentype_index) > 0 && (dentype_index + 1) <= length(name_parts)) {
            name_parts[dentype_index + 1]
        } else {
            NA  # Return NA if the pattern doesn't match
        }
        
        # Return a list containing the method, dentype, and the dataframe
        list(
            method = method,
            dentype = dentype,
            data = metrics.list[[name]]
        )
    })
    
    # Return the list of extracted information
    return(info_list)
}

# Plot RMSE for a/b parameters TOGETHER
plot.rmse.diff <- function(metrics_list, parameter = "a") {
    if (!parameter %in% c("a", "b")) {
        stop("Parameter must be either 'a' or 'b'.")
    }
    # Use the extract.info.from.list function to structure the input data
    structured_list <- extract.info.from.list(metrics_list)
  
    # Combine all structured data into a single dataframe
    combined_data <- bind_rows(lapply(structured_list, function(info) {
        df <- info$data
        rmse_col <- paste("rmse", parameter, sep = ".")
        if (!rmse_col %in% names(df)) {
            stop(paste("The dataframe does not contain '", rmse_col, "' column.", sep = ""))
        }
        df %>%
            select(item, !!rmse_col) %>%
            rename(rmse = !!rmse_col) %>%
            mutate(method = info$method, dentype = info$dentype, item = as.factor(item))
    }), .id = "id") %>%
    select(-id) %>%
    mutate(method = factor(method, levels = c("BL", "EM", "MCEM")))  # Reorder methods so "BL" is plotted first 

    
    # Dynamically calculate the number of columns for a symmetric grid if needed
    n_plots <- length(unique(combined_data$item))
    n_cols <- ceiling(sqrt(n_plots))
    
    # Plotting RMSE differences across methods for a single parameter
    p <- ggplot(combined_data, aes(x = item, y = rmse, group = interaction(method, dentype), color = method)) +
        geom_line(aes(linewidth = method)) +
        geom_point(aes(shape = method)) +
        scale_linewidth_manual(values = c("EM" = 1, "BL" = 2.5, "MCEM" = 1)) +  # Adjust line widths
        facet_wrap(~ dentype, scales = "free_y", ncol = n_cols) +
        labs(title = paste("RMSE differences for parameter '", parameter, "' across methods", sep = ""),
             x = "Item", y = paste("RMSE for '", parameter, "'", sep = "")) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "bottom") +
        scale_color_manual(values = c("EM" = "blue", "BL" = "red", "MCEM" = "green"))

    return(p)
}
# Bias for a/b parameters -------------------------------
plot.bias.diff <- function(metrics_list, parameter = "a") {
    if (!parameter %in% c("a", "b")) {
        stop("Parameter must be either 'a' or 'b'.")
    }
    # Use the extract.info.from.list function to structure the input data
    structured_list <- extract.info.from.list(metrics_list)
  
    # Combine all structured data into a single dataframe
    combined_data <- bind_rows(lapply(structured_list, function(info) {
        df <- info$data
        # Select only the relevant bias column based on param
        bias_col <- paste("bias", parameter, sep = ".")
        if (!bias_col %in% names(df)) {
            stop(paste("The dataframe does not contain '", bias_col, "' column.", sep = ""))
        }
        # Keep only necessary columns and add method and dentype
        df %>%
            select(item, !!bias_col) %>%
            rename(bias = !!bias_col) %>%
            mutate(method = info$method, dentype = info$dentype, item = as.factor(item))
    }), .id = "id") %>%
    select(-id) %>%
    mutate(method = factor(method), dentype = factor(dentype))
    
    # Dynamically calculate the number of columns for a symmetric grid if needed
    n_plots <- length(unique(combined_data$item))
    n_cols <- ceiling(sqrt(n_plots))
    
    # Plotting Bias differences across methods for a single parameter
    p <- ggplot(combined_data, aes(x = item, y = bias, group = interaction(method, dentype), color = method)) +
        geom_line(aes(size = method)) +  # Control line size based on method
        scale_size_manual(values = c("EM" = 1, "BL" = 2, "MCEM" = 1)) +  # Customize line width
        geom_point(aes(shape = method)) +
        facet_wrap(~ dentype, scales = "free_y", ncol = n_cols) +
        labs(title = paste("Bias differences for parameter '", parameter, "' across methods", sep = ""),
             x = "Item", y = paste("Bias for '", parameter, "'", sep = "")) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "bottom") +
        scale_color_manual(values = c("EM" = "blue", "BL" = "red", "MCEM" = "green")) +
        guides(size = "none")  # Hide the size guide since it's being used to differentiate lines
    
    return(p)
}
# testing 
# # rmse
# p1 <- plot.rmse.diff(metrics.stnd, parameter = "a")
# p2 <- plot.rmse.diff(metrics.stnd, parameter = "b")
# p.rmse <- ggarrange(p1, p2, ncol = 1, nrows = 2)
# # bias
# p1 <- plot.bias.diff(metrics.stnd, parameter = "a")
# p2 <- plot.bias.diff(metrics.stnd, parameter = "a")
# p.bias <- ggarrange(p1, p2, ncol = 1, nrows = 2)

# Bar Charts --------------------------------------------
# RMSE for a/b parameters -------------------------------
plot.rmse.bar <- function(metrics_1, metrics_2) {
    metrics_1$item <- as.character(metrics_1$item)
    metrics_2$item <- as.character(metrics_2$item)

    cross_param <- metrics_1 %>%
        select(item, true_a = true.a, true_b = true.b) %>%
        distinct(item, .keep_all = TRUE)

    combined_metrics <- full_join(metrics_1, metrics_2, by = "item", suffix = c("_1", "_2"))

    combined_long <- combined_metrics %>%
        pivot_longer(cols = c("rmse.a_1", "rmse.a_2", "rmse.b_1", "rmse.b_2"),
                                 names_to = "Parameter_Method", values_to = "Value") %>%
        mutate(Method = ifelse(str_detect(Parameter_Method, "_1"), "1", "2"),
                     Parameter = ifelse(str_detect(Parameter_Method, "a"), "a", "b")) %>%
        left_join(cross_param, by = "item") %>%
        mutate(AB_Label = factor(paste("a =", true_a, "b =", true_b), levels = unique(paste("a =", true_a, "b =", true_b))))

    # Calculate the maximum RMSE value and add a buffer for labels
    max_value <- max(combined_long$Value, na.rm = TRUE)
    label_space <- max_value * 0.1 # Adjust this factor to increase/decrease space for labels
    
    p1 <- ggplot(combined_long, aes(x = Parameter, y = Value, fill = Method)) +
        geom_col(position = position_dodge(width = 0.8)) +
        facet_wrap(~ AB_Label, scales = "fixed", ncol = length(unique(combined_long$true_a))) +
        geom_text(aes(label = Parameter), position = position_dodge(width = 0.8), vjust = -0.25, size = 4.5, angle = 0) +
        scale_fill_manual(values = c("blue", "red")) +
        labs(title = "RMSE for Each Item", x = "Parameter", y = "RMSE") +
        theme_minimal() +
        theme(legend.position = "bottom",
              strip.background = element_blank(),
              strip.text.x = element_text(size = 8)) +
        ylim(0, max_value + label_space) # Dynamically set y-axis limits

    print(p1)
    return(p1)
}
# test:
# plot.rmse.bar(metrics.stnd[[1]], metrics.stnd[[1]])
# plot.bias.bar(metrics.stnd[[1]], metrics.stnd[[1]])

# Bias for a/b parameters -------------------------------
plot.bias.bar <- function(metrics_1, metrics_2) {
    metrics_1$item <- as.character(metrics_1$item)
    metrics_2$item <- as.character(metrics_2$item)

    cross_param <- metrics_1 %>%
        select(item, true_a = true.a, true_b = true.b) %>%
        distinct(item, .keep_all = TRUE)

    combined_metrics <- full_join(metrics_1, metrics_2, by = "item", suffix = c("_1", "_2"))
    combined_long <- combined_metrics %>%
        pivot_longer(cols = c("bias.a_1", "bias.a_2", "bias.b_1", "bias.b_2"),
                                 names_to = "Parameter_Method", values_to = "Value") %>%
        mutate(Method = ifelse(str_detect(Parameter_Method, "_1"), "1", "2"),
                     Parameter = case_when(
                         str_detect(Parameter_Method, "a_") ~ "a",
                         str_detect(Parameter_Method, "b_") ~ "b",
                         TRUE ~ NA_character_)) %>%
        left_join(cross_param, by = "item") %>%
        mutate(AB_Label = factor(paste("a =", true_a, "b =", true_b), levels = unique(paste("a =", true_a, "b =", true_b))),
                     Text_Y_Pos = ifelse(Value < 0, 0, Value))  # Determine Y position based on Value
    
    # Calculate the maximum Bias value and add a buffer for labels
    max_value <- max(abs(combined_long$Value), na.rm = TRUE) # Using abs() to account for negative values as well
    label_space <- max_value * 0.1 # Adjust this factor to increase/decrease space for labels

    p2 <- ggplot(combined_long, aes(x = Parameter, y = Value, fill = Method)) +
        geom_col(position = position_dodge(width = 0.8)) +
        facet_wrap(~ AB_Label, scales = "fixed", ncol = length(unique(combined_long$true_a))) +
        geom_text(aes(label = Parameter, y = Text_Y_Pos),
                  position = position_dodge(width = 0.8), vjust = -0.75, size = 4.5, angle = 0) +
        scale_fill_manual(values = c("blue", "red")) +
        labs(title = "Bias for Each Item", x = "Parameter", y = "Bias") +
        theme_minimal() +
        theme(legend.position = "bottom",
              strip.background = element_blank(),
              strip.text.x = element_text(size = 8)) +
        expand_limits(y = c(-max_value * 0.1, max_value + label_space)) # Adjust y-axis limits
    
    print(p2)
    return(p2)
}
# test:
# p.rmse <- plot.rmse.bar(metrics.stnd[[1]], metrics.stnd[[1]])
# p.bias <- plot.bias.bar(metrics.stnd[[1]], metrics.stnd[[1]])

# ----------------------------------------------------
# old functions --------------------------------------
# ----------------------------------------------------

# plots rmse and bias differences between different methods
plot.rmse.bias.differences <- function(metrics_list) {    
 # Extract information for each dataframe in the list
  metrics_info <- extract.info.from.list(metrics_list)
  
  # Combine all dataframes into one, ensuring 'Method' and 'Dentype' are included
  combined_metrics <- do.call(rbind, lapply(metrics_info, function(info) {
    cbind(info$data, Method = info$method, Dentype = info$dentype)
  }))
  
  # Make sure all dataframe elements have a uniform 'item' column
  combined_metrics$item <- factor(combined_metrics$item)
  
  # Pivot the data to a long format suitable for plotting
  combined_long <- combined_metrics %>%
    pivot_longer(cols = starts_with("rmse") | starts_with("bias"),
                 names_to = c(".value", "Parameter"), 
                 names_pattern = "(.*)\\.(.*)") 
  
  # Separate the parameter and the metric type (RMSE/Bias)
  combined_long <- combined_long %>%
    separate(Parameter, into = c("Metric_Type", "Parameter"), sep = "_")
  
  # Filter out the Parameter you want to plot (e.g., "a" or "b")
  combined_long <- combined_long %>%
    filter(Parameter == "a") # Change to "b" if you want to plot "b"
  
  # Create a label for each bar with the method, metric type, and parameter
  combined_long$Label <- paste(combined_long$Method, combined_long$Metric_Type, combined_long$Parameter, sep = ": ")
  
  # Plot the RMSE and Bias for each parameter-method-dentype combination
  p <- ggplot(combined_long, aes(x = Metric_Type, y = Value, fill = Method)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.7)) +
    geom_text(aes(label = Label, y = ifelse(Value > 0, Value + 0.02, Value - 0.02)), 
              position = position_dodge(width = 0.7),
              size = 3, angle = 90, hjust = 1, vjust = 1) +
    facet_wrap(~ item, scales = "free_y", ncol = 4) +
    scale_fill_brewer(palette = "Set1") +
    labs(title = "RMSE and Bias for Each Item by Method and Dentype",
         x = "", y = "Value") +
    theme_minimal() +
    theme(legend.position = "bottom",
          strip.background = element_blank(),
          strip.text.x = element_text(size = 8, angle = 45, hjust = 1),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          axis.ticks.y = element_line(color = "black", size = 1.2)) +
    scale_y_continuous(limits = c(NA, NA), oob = scales::rescale_none) # Use oob to handle out of bounds

  return(p)
}
# plots single theta distribution (dist.data)
plot.dist <- function(dist.data) {
    # Ensure 'theta' and 'distribution' columns are present in the dataframe
    if (!("theta" %in% names(dist.data)) | !("distribution" %in% names(dist.data))) {
        stop("dist.data must contain 'theta' and 'distribution' columns")
    }

    # Generate the plot
    p <- ggplot(dist.data, aes(x = theta, color = distribution)) +
            geom_density() +
            labs(x = "Theta", y = "Density") +
            theme_minimal() +
            theme(legend.title = element_blank()) # Hide the legend title if desired

    return(p)
} # end plot.dist
# show all three distributions (dist.list) in one plot
plot.dist.combined <- function(dist.list) {
    # Check if the list is not empty
    if (length(dist.list) == 0) {
        stop("The input list is empty.")
    }
  
    # Combine all distribution data into one dataframe
    combined.data <- do.call(rbind, lapply(seq_along(dist.list), function(i) {
        # Ensure 'theta' column is present
        if (!"theta" %in% names(dist.list[[i]])) {
            stop("Each dataframe in dist.list must contain a 'theta' column.")
        }
        # Add or overwrite the 'distribution' column with the correct name
        dist.data <- dist.list[[i]]
        dist.data$distribution <- names(dist.list)[i]
        return(dist.data)
    }))
    
    # Ensure 'distribution' column is treated as a factor
    combined.data$distribution <- factor(combined.data$distribution, levels = names(dist.list))
    
    # Plot
    p <- ggplot(combined.data, aes(x = theta, color = distribution)) +
        geom_density() +
        labs(x = "Theta", y = "Density", color = "Distribution") +
        scale_color_brewer(palette = "Set1") +
        theme_minimal() +
        theme(legend.title = element_blank())
    
    print(p)
} # end plot.dist.combined
# ----------------------------------------------------

load()
n <- 300
# Method 1 ----------------------------------------------------
all.distributions <- generate.skewed.distribitions(n, seed=123)

# cros.param has 20 rows and 3 columns (of true param values for a, b, and d)
cross.param <- expand.grid(a = c(0.5, 1, 1.5, 2.5),b = c(-2.5, -1.25, 0, 1.25, 2.5))
# cross.param$b <- with(cross.param, -d/a)

# simulate response data: for each dist, has 300 rows of 20 item responses
# response.dataframes <- simulate.response.data(all.distributions, cross.param, seed = 123) # obsolete
sim <- simulate.response.data(all.distributions, cross.param)
response.dataframes <- sim$response
probabilities.dataframes <- sim$probabilities
# output here: response.dataframes$left.skew, response.dataframes$right.skew, response.dataframes$stnd.norm

# Method 2 ----------------------------------------------------
dist.data <- quick.gen.dist(300, dist.type='stnd.norm', seed=123)
dist.stnd <- dist.data
dist.right <-  quick.gen.dist(300, dist.type='right.skew', seed=123)
dist.left <- quick.gen.dist(300, dist.type='left.skew', seed=123)
# so far, dist.stnd, dist.right, dist.left are dataframes with 300 rows of theta values and distribution type
sim <- quick.sim.response(dist.data$theta, cross.param, seed=123)
response.data <- sim$response
probabilities.data <- sim$probabilities
export.data(cross.param, output.dir='response_data', n=300, replications=10, seed=123) # calls quick.sim.response and quick.gen.dist

# Methods/Dentypes/Distributions -------------------------------
methods <- c("BL", "EM")
dentypes <- c("Gaussian")
dist.types <- c("stnd.norm", "left.skew", "right.skew")

# MIRT Models ---------------------------------------------------
# Methods:
#   BL: Bayesian likelihood: common MLE
#   EM: Expectation Maximization (EM): common MLE
#   MCEM: Monte Carlo EM: variant of the EM algorithm that uses Monte Carlo simulation for the E-step.
#   QMCEM: Quasi Monte Carlo EM: quasi-random sequences to improve efficiency of the Monte Carlo integration in the E-step
#   MHRM: Metropolis Hastings Random Walk Metropolis: stochastic, combines the Robbins-Monro(RH) stochastic approximation with Metropolis-Hastings (MH) sampling.
#   SEM: Stochastic EM: stochastic version EM algorithm
methods <- c("BL", "EM", "MCEM", "QMCEM", "MHRM", "SEM")

# Dentypes:
# Gaussian ("Gaussian"): normal distribution: assumes a multivariate Gaussian w/ assoc. mean vector and var-covariance matrix
# empiricalhist ("EH"): empirical histogram: estimates latent distribution using the empirical histogram of the data
# empiricalhist_Woods ("EHW"): est. latent distribution using the empirical histogram of the data w/ Woods correction (extrapolation-interpolation)
#   - improve stability in the presence of extreme response styles (i.e., all highest or lowest in each item)
# Davidian-# ("Davidian-2"): Davidian distribution: 2-parameter distribution
#  - estimates semi-parametric Davidian curves described by Woods
#  - the '#' placeholder, assumes a #-parameter distribution for the latent trait 
dentypes <- c("Gaussian", "EH", "EHW", "Davidian-2", "Davidian-3")

# Simulations ---------------------------------------------------
metrics.stnd <- fit.mirt(dist.type=dist.types[1], cross.param, methods, dentypes, 10) 
metrics.left <- fit.mirt(dist.type=dist.types[2], cross.param, methods, dentypes, 10) 
metrics.right <- fit.mirt(dist.type=dist.types[3], cross.param, methods, dentypes, 10)

# Error:
    # throwing error for response_data/response_left_skew_rep_01.txt where "I20" contains all 0's
    # Solution: 
    # Bayesian Priors: https://cran.r-project.org/web/packages/mirt/vignettes/priors.html
    # Handle Error: https://stackoverflow.com/questions/64669882/how-to-handle-error-in-r-when-reading-data-from-a-file

# all metrics
metrics <- list(stnd.norm = metrics.stnd, left.skew = metrics.left, right.skew = metrics.right)
metrics 

# RMSE ----------------------------------------------------
p1 <- plot.rmse.diff(metrics.stnd, parameter="a")
p2 <- plot.rmse.diff(metrics.stnd, parameter="b")
plot.rmse <- ggpubr::ggarrange(p1, p2, nrow = 2, ncol = 1)
print(plot.rmse)
# ggsave("~/Desktop/plot_rmse.png", plot.rmse)

# BIAS ----------------------------------------------------
p1 <- plot.bias.diff(metrics.stnd, parameter="a")
p2 <- plot.bias.diff(metrics.stnd, parameter="b")
plot.bias <- ggpubr::ggarrange(p1, p2, nrow = 2, ncol = 1)
print(plot.bias)
# ggsave("~/Desktop/plot_bias.png", plot.bias)

# HISTOGRAMS ----------------------------------------------------
# visualize the rmse and bias differences between different methods 
p.rmse <- plot.rmse.bar(metrics.stnd[[1]], metrics.stnd[[1]])
p.bias <- plot.bias.bar(metrics.stnd[[1]], metrics.stnd[[1]])
# ggsave("~/Desktop/plot_rmse_bar.png", p.rmse)
# ggsave("~/Desktop/plot_bias_bar.png", p.bias)

# ----------------------------------------------------