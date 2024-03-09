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


# 3/9/23 RMSE/Bias plot compares method for each item

# working
plot.rmse.differences <- function(metrics.BL, metrics.EM) {
    # join metrics 
    combined.metrics <- full_join(metrics.BL, metrics.EM, by = "item", suffix = c("_BL", "_EM"))
  
    # convert to longer format
    combined.long <- combined.metrics %>%
        pivot_longer(cols = c("rmse.a_BL", "rmse.a_EM", "rmse.b_BL", "rmse.b_EM"), 
                     names_to = "Parameter_Method", values_to = "RMSE") %>% 
        mutate(Method = ifelse(str_detect(Parameter_Method, "_BL"), "BL", "EM"), # extract method
               Parameter = ifelse(str_detect(Parameter_Method, "a_"), "a", "b"), # extract parameter
               RMSE_Diff = RMSE * ifelse(Method == "BL", 1, -1)) %>% # calc difference
        select(-Parameter_Method) # remove intermediate column

    # separate the true values for annotations
    true.values <- combined.metrics %>% 
        select(item, starts_with("true")) %>% 
        pivot_longer(cols = starts_with("true"), names_to = "True_Param", values_to = "True_Value") %>% 
        mutate(Parameter = ifelse(str_detect(True_Param, "a"), "a", "b")) %>% # extract parameter
        distinct(item, Parameter, .keep_all = TRUE)  # ensure unique item-parameter pairs

    # join the true values for plotting
    combined.long <- left_join(combined.long, true.values, by = c("item", "Parameter")) 

    # plot the RMSE differences with true values annotated
    p <- ggplot(combined.long, aes(x = Parameter, y = RMSE_Diff, fill = Method)) +
        geom_col(position = position_dodge(width = 0.8)) +
        facet_wrap(~ item, scales = "fixed", labeller = label_both) +  # Set fixed scales and label each plot with "Item"
        geom_text(aes(label = sprintf("True %s: %.2f", Parameter, True_Value)),
                  position = position_dodge(width = 0.8), check_overlap = TRUE,
                  vjust = -0.5, size = 2.5) +
        scale_fill_manual(values = c("blue", "red")) +
        labs(title = "RMSE Comparison for Each Item: BL vs. EM", x = "Parameter", y = "RMSE Difference") +
        theme_minimal() +
        theme(legend.position = "bottom",
              strip.background = element_blank(),
              strip.text.x = element_text(size = 8)) + 
        ylim(-0.5, 0.5) # Set fixed y-axis limits

    print(p)
}

plot.bias.differences <- function(metrics.BL, metrics.EM){
    # combine metrics
    combined.metrics <- full_join(metrics.BL, metrics.EM, by = "item", suffix = c("_BL", "_EM"))
  
    # convert to long format for easier plotting of bias
    combined.long <- combined.metrics %>%
        pivot_longer(cols = c("bias.a_BL", "bias.a_EM", "bias.b_BL", "bias.b_EM"), 
                     names_to = "Parameter_Method", values_to = "Bias") %>%
        mutate(Method = ifelse(str_detect(Parameter_Method, "_BL"), "BL", "EM"),
               Parameter = ifelse(str_detect(Parameter_Method, "a_"), "a", "b"), # extract parameter
               Bias_Diff = Bias * ifelse(Method == "BL", 1, -1)) %>% # calculate difference
        select(-Parameter_Method) # remove intermediate column

    # prep true values for annotations
    true.values <- combined.metrics %>%
        select(item, starts_with("true")) %>%
        pivot_longer(cols = starts_with("true"), names_to = "True_Param", values_to = "True_Value") %>%
        mutate(Parameter = ifelse(str_detect(True_Param, "a"), "a", "b")) %>%
        distinct(item, Parameter, .keep_all = TRUE)  # Ensure unique item-parameter pairs

    # join the true values for plotting
    combined.long <- left_join(combined.long, true.values, by = c("item", "Parameter"))

    # plot the Bias differences with true values annotated
    p <- ggplot(combined.long, aes(x = Parameter, y = Bias_Diff, fill = Method)) +
        geom_col(position = position_dodge(width = 0.8)) +
        facet_wrap(~ item, scales = "fixed", labeller = label_both) +
        geom_text(aes(label = sprintf("True %s: %.2f", Parameter, True_Value)),
                  position = position_dodge(width = 0.8), check_overlap = TRUE,
                  vjust = -0.5, size = 2.5) +
        scale_fill_manual(values = c("blue", "red")) +
        labs(title = "Bias Comparison for Each Item: BL vs. EM", x = "Parameter", y = "Bias Difference") +
        theme_minimal() +
        theme(legend.position = "bottom",
              strip.background = element_blank(),
              strip.text.x = element_text(size = 8)) + 
        ylim(-0.5, 0.5)  # Set the same scale for all plots

    print(p)
}


# 3/9/23 wanted to plot rmse using *args, passing the entire list into function 

plot.rmse.differences.test <- function(metrics.list) {
    combined_metrics <- bind_rows(metrics.list, .id = "Method") %>%
        mutate(Item = factor(item), 
               Method = factor(Method, levels = names(metrics.list)))

    # Prepare the column names for RMSE based on the presence of "rmse."
    rmse_cols <- names(combined_metrics)[grepl("rmse", names(combined_metrics))]

    # Pivot to longer format for easier plotting
    combined_long <- combined_metrics %>%
        pivot_longer(cols = all_of(rmse_cols), names_to = "Parameter_Type", values_to = "Value") %>%
        mutate(
            Parameter = ifelse(str_detect(Parameter_Type, "a"), "a", "b"),
            Interaction = interaction(Item, Parameter, sep = ": ")
        )

    # Extract true values for annotation
    true_values <- combined_metrics %>%
        select(Item, starts_with("true")) %>%
        pivot_longer(cols = starts_with("true"), names_to = "True_Parameter", values_to = "True_Value") %>%
        mutate(Parameter = ifelse(str_detect(True_Parameter, "a"), "a", "b"))

    # Plot the RMSEs for all methods
    p <- ggplot(combined_long, aes(x = Interaction, y = Value, fill = Method)) +
        geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
        facet_wrap(~ Item, ncol = 1, scales = "free_y") +
        scale_fill_brewer(palette = "Set1") +
        labs(title = "RMSE Comparison for Each Item: Multiple Methods", 
             x = "Item and Parameter", y = "RMSE") +
        theme_minimal() +
        theme(
            legend.position = "bottom",
            strip.background = element_blank(),
            strip.text.x = element_text(size = 8),
            strip.text.y = element_text(hjust = 0),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
        ) +
        ylim(-0.5, 0.5) # Set fixed y-axis limits
    
    # Add the true values to the plot
    p <- p + geom_text(
        data = true_values,
        aes(label = sprintf("True %s: %.2f", Parameter, True_Value),
            y = ifelse(Value < 0, min(Value, na.rm = TRUE) * 1.1, max(Value, na.rm = TRUE) * 0.9),
            group = Interaction),
        position = position_dodge(width = 0.8),
        size = 2.5, check_overlap = TRUE,
        vjust = ifelse(Value < 0, 1, -1)
    )

    print(p)
}



