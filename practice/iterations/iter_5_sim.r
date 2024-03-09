rm(list=ls())

load <- function(){
    packages <- c("mirt", "ggplot2", "reshape2", "tibble", "ggpubr", "gridExtra", "sn", "tidyr", "dplyr", "moments", "parallel", "dplyr", "stringr")

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

#### Method 2
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
                    coefs <- coef(fit, simplify = TRUE, IRTpars = TRUE )$items # adjusted 

                    est.params.a[rep, ] <- coefs[, "a"]  # Corrected from "a1" to "a"
                    est.params.b[rep, ] <- coefs[, "b"]  # Corrected from "d" to "b"

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

                metrics$est.a <- colMeans(est.params.a, na.rm = TRUE)
                metrics$est.b <- colMeans(est.params.b, na.rm = TRUE)

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

# visuals: 439, 442, 445 ?

# plot.characteristic.curves <- function(dist.data, probabilities.data, metrics){
#     # Ensure the presence of necessary columns
#     if (!"theta" %in% names(dist.data)) {
#         stop("dist.data must contain a 'theta' column.")
#     }
    
#     # Ensure the number of rows match
#     if (nrow(dist.data) != nrow(probabilities.data)) {
#         stop("The number of rows in dist.data and probabilities must be the same.")
#     }

#     # print the distribution in dist.data to be aware
#     message(paste("Please check to ensure that dist.type is correct:", names(all.distributions)))


#     # combine the theta and probabilites into 1 dataframe
#     combined <- cbind(dist.data, probabilities.data) %>%
#         pivot_longer(cols = starts_with("P"), names_to = "item", values_to = "probability") %>%
#         mutate(item = as.numeric(gsub("P", "", item)))

#     # subset to only include the items corresponding to "b" parameter 
#     message(paste("Fixing parameter a = 1.5"))
#     items.to.plot <- cross.param[cross.param$a == 1.5, "b"]
#     combined <- combined %>% filter(item %in% items.to.plot)
    
#     calc_prob_2pl <- function(theta, a, b) {
#         1 / (1 + exp(-a * (theta - b)))
#     }
    
#     curve_data <- metrics %>%
#         rowwise() %>%
#         do({
#             data.frame(
#                 theta = dist.data$theta,
#                 probability_true = calc_prob_2pl(dist.data$theta, .$true.a, .$true.b),
#                 probability_est = calc_prob_2pl(dist.data$theta, .$est.a, .$est.b),
#                 item = as.factor(.$item)
#             )
#         }) %>%
#         ungroup()

#     p <- ggplot() +
#         geom_line(data = combined, aes(x = theta, y = probability, group = item, color = item), alpha = 0.3) +
#         geom_line(data = curve_data, aes(x = theta, y = probability_true, color = item), linewidth = 1.2) +
#         geom_line(data = curve_data, aes(x = theta, y = probability_est, color = item), linewidth = 0.6, linetype = "dashed") +
#         scale_color_manual(values = rainbow(length(unique(curve_data$item)))) +
#         labs(title = "Item Characteristic Curves: True vs. Estimated", x = "Theta", y = "Probability") +
#         theme_minimal() +
#         theme(legend.position = "none")

#     print(p)
    
# } # end characteristic.curves

# plot.icc <- function(cross.param, metrics, dist.data, probabilities.data) {

# }


plot.method.comparison.rmse <- function(metrics.BL, metrics.EM){
    # ensure arguments are of the form metrics[[1]] for proper handling of dataframes
    metrics.BL <- metrics.BL %>% mutate(Method = "BL")
    metrics.EM <- metrics.EM %>% mutate(Method = "EM")
    combined.metrics <- rbind(metrics.BL, metrics.EM)

    # convert items into factors for proper ordering
    combined.metrics$item <- as.factor(combined.metrics$item)

    # into long format
    long.metrics <- combined.metrics %>% 
        gather(key="Parameter", value="RMSE", rmse.a, rmse.b) %>%
        mutate(Parameter = gsub("rmse.a.", "", Parameter)) %>%
        mutate(Difference = ifelse(Method == "BL", RMSE, -RMSE))

    p <- ggplot(long.metrics, aes(x = Parameter, y = Difference, fill = Method)) + 
        geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
        facet_wrap(~ item, scales = "free_y") +
        labs(title = "RMSE Comparison for Each Item: BL vs. EM", 
             x = "Parameter", y = "RMSE Difference") +
        theme_minimal() +
        theme(legend.position = "bottom",
              axis.text.x = element_text(angle = 45, hjust = 1),
              strip.background = element_blank(),
              strip.text.x = element_text(size = 8))
        
    print(p)
    
}

plot.rmse.differences.with.true.values <- function(metrics.BL, metrics.EM) {
    combined.metrics <- full_join(metrics.BL, metrics.EM, by = "item", suffix = c("_BL", "_EM"))
  
    # Convert to long format for easier plotting
    combined.long <- combined.metrics %>%
        pivot_longer(cols = c("rmse.a_BL", "rmse.a_EM", "rmse.b_BL", "rmse.b_EM"), 
                     names_to = "Parameter_Method", values_to = "RMSE") %>%
        mutate(Method = ifelse(str_detect(Parameter_Method, "_BL"), "BL", "EM"),
               Parameter = ifelse(str_detect(Parameter_Method, "a_"), "a", "b"),
               RMSE_Diff = RMSE * ifelse(Method == "BL", 1, -1)) %>%
        select(-Parameter_Method)

    # Create a single row per item and parameter in true.values for annotations
    true.values <- combined.metrics %>%
        select(item, starts_with("true")) %>%
        pivot_longer(cols = starts_with("true"), names_to = "True_Param", values_to = "True_Value") %>%
        mutate(Parameter = ifelse(str_detect(True_Param, "a"), "a", "b")) %>%
        distinct(item, Parameter, .keep_all = TRUE)  # Ensure unique item-parameter pairs

    # Join the true values for plotting
    combined.long <- left_join(combined.long, true.values, by = c("item", "Parameter"))

    # Plot the RMSE differences with true values annotated
    p <- ggplot(combined.long, aes(x = Parameter, y = RMSE_Diff, fill = Method)) +
        geom_col(position = position_dodge(width = 0.8)) +
        facet_wrap(~ item, scales = "free_y", labeller = label_bquote(Item~.x)) +  # Add 'Item' label to the top of each plot
        geom_text(aes(label = sprintf("True %s: %.2f", Parameter, True_Value)),
                  position = position_dodge(width = 0.8), check_overlap = TRUE,
                  vjust = -0.5, size = 2.5) +
        scale_fill_manual(values = c("blue", "red")) +
        labs(title = "RMSE Comparison for Each Item: BL vs. EM", x = "", y = "RMSE Difference") +
        theme_minimal() +
        theme(legend.position = "bottom",
              strip.background = element_blank(),
              strip.text.x = element_text(size = 8)) +
        ylim(-0.5, 0.5)  # Set y-axis limits to -0.5 to 0.5


    print(p)
}

load()
n <- 300
# Method 1: initial implementation
all.distributions <- generate.skewed.distribitions(n, seed=123)

# cros.param has 20 rows and 3 columns (of true param values for a, b, and d)
cross.param <- expand.grid(d = c(-2.5, -1.25, 0, 1.25, 2.5), a = c(0.5, 1, 1.5, 2.5))
cross.param$b <- with(cross.param, -d/a)

# simulate response data: for each dist, has 300 rows of 20 item responses
# response.dataframes <- simulate.response.data(all.distributions, cross.param, seed = 123)
sim <- simulate.response.data(all.distributions, cross.param)
response.dataframes <- sim$response
probabilities.dataframes <- sim$probabilities

# output here: response.dataframes$left.skew, response.dataframes$right.skew, response.dataframes$stnd.norm

######
# Method 2: revised implementation 
dist.data <- quick.gen.dist(300, dist.type='stnd.norm', seed=123)

sim <- quick.sim.response(dist.data$theta, cross.param, seed=123)
response.data <- sim$response
probabilities.data <- sim$probabilities
export.data(cross.param, output.dir='response_data', n=300, replications=10, seed=123) # calls quick.sim.response and quick.gen.dist

methods <- c("BL", "EM")
dentypes <- c("Gaussian")
dist.types <- c("stnd.norm", "left.skew", "right.skew")


# Corrected function call to specifically use 'stnd.norm' dataframe
# no rep parameter
metrics.stnd <- fit.mirt(dist.type=dist.types[1], cross.param, methods, dentypes, 10) 
metrics.left <- fit.mirt(dist.type=dist.types[2], cross.param, methods, dentypes, 10) 
metrics.right <- fit.mirt(dist.type=dist.types[3], cross.param, methods, dentypes, 10)

metrics <- list(stnd.norm = metrics.stnd, left.skew = metrics.left, right.skew = metrics.right)
metrics # view

# testing
# plot.icc(cross.param, metrics.stnd, dist.data, probabilities.data)
plot.method.comparison.rmse(metrics.stnd[[1]], metrics.stnd[[2]])

plot.rmse.differences.with.true.values(metrics.stnd[[1]], metrics.stnd[[2]])
