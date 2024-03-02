a <- c(1, 2, 3, 4)
b <- c(5, 6, 7, 8)

params <- data.frame(a, b)

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
        # p <- plot_skewed_distributions(all_distributions)
        # print(p) # Explicitly print the plot
    }
    return(all_distributions)
} # end generate_skewed_distribitions


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


fit_mirt_models <- function(response.dataframes, true_params) {
    results <- list()
    mirt.models <- list()

    for (name in names(response.dataframes)) {
        # extract response data
        response.data <- response.dataframes[[name]]

        # fit the 2PL model
        start.time <- Sys.time()
        mirt.model <- mirt(data = response.data, model = 1, itemtype = "2PL", storeEMhistory = TRUE)
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


all.distributions <- generate_skewed_distribitions(500, 123, FALSE)
response.dataframes <- simulate_response_data(all.distributions, params, 123)
names(response.dataframes)

fit <- fit_mirt_models(response.dataframes, params)

results <- fit$Results
models <- fit$MirtModels
