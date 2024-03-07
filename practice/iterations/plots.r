# Plots all item characteristic curves for the 2PL model
# MESSY
plot.characteristic.curves <- function(cross.param, dist.data, probabilities.data){

    # check if theta in dist.data 
    if (!"theta" %in% names(dist.data)) {
        stop("dist.data must contain a 'theta' column.")
    }
    # print the distribution in dist.data to be aware
    message(paste("Please check to ensure that dist.type is correct:", names(all.distributions)))

    # ensure nrows of dist.data and probabilities.data are equal
    if (nrow(dist.data) != nrow(probabilities.data)) {
        stop("The number of rows in dist.data and probabilities must be the same.")
    }

    # combine the theta and probabilites into 1 dataframe
    combined <- cbind(dist.data, probabilities.data) %>% 
        pivot_longer(
        cols = starts_with("P"), 
        names_to = "item", 
        values_to = "probability"
        ) %>% 
        mutate(item = as.numeric(gsub("P", "", item)))  # Convert factor levels to numeric if needed

    p <- ggplot(combined, aes(x = theta, y = probability, group = item, color = as.factor(item))) +
        geom_line() +
        labs(title = "Item Characteristic Curves", x = "Theta", y = "Probability", color = "Item") +
        theme_minimal() +
        theme(legend.position = "none")  # Remove the legend if too many lines
        
    print(p)
    return(p)
} # end characteristic.curves


