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


# 3/9/24 RMSE/Bias plot compares method for each item

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


# 3/9/24 Implemented icc curves using prob.data and dist.data, and metrics 

plot.characteristic.curves <- function(dist.data, probabilities.data, metrics){
    # Ensure the presence of necessary columns
    if (!"theta" %in% names(dist.data)) {
        stop("dist.data must contain a 'theta' column.")
    }
    
    # Ensure the number of rows match
    if (nrow(dist.data) != nrow(probabilities.data)) {
        stop("The number of rows in dist.data and probabilities must be the same.")
    }

    # print the distribution in dist.data to be aware
    message(paste("Please check to ensure that dist.type is correct:", names(all.distributions)))


    # combine the theta and probabilites into 1 dataframe
    combined <- cbind(dist.data, probabilities.data) %>%
        pivot_longer(cols = starts_with("P"), names_to = "item", values_to = "probability") %>%
        mutate(item = as.numeric(gsub("P", "", item)))

    # subset to only include the items corresponding to "b" parameter 
    message(paste("Fixing parameter a = 1.5"))
    items.to.plot <- cross.param[cross.param$a == 1.5, "b"]
    combined <- combined %>% filter(item %in% items.to.plot)
    
    calc_prob_2pl <- function(theta, a, b) {
        1 / (1 + exp(-a * (theta - b)))
    }
    
    curve_data <- metrics %>%
        rowwise() %>%
        do({
            data.frame(
                theta = dist.data$theta,
                probability_true = calc_prob_2pl(dist.data$theta, .$true.a, .$true.b),
                probability_est = calc_prob_2pl(dist.data$theta, .$est.a, .$est.b),
                item = as.factor(.$item)
            )
        }) %>%
        ungroup()

    p <- ggplot() +
        geom_line(data = combined, aes(x = theta, y = probability, group = item, color = item), alpha = 0.3) +
        geom_line(data = curve_data, aes(x = theta, y = probability_true, color = item), linewidth = 1.2) +
        geom_line(data = curve_data, aes(x = theta, y = probability_est, color = item), linewidth = 0.6, linetype = "dashed") +
        scale_color_manual(values = rainbow(length(unique(curve_data$item)))) +
        labs(title = "Item Characteristic Curves: True vs. Estimated", x = "Theta", y = "Probability") +
        theme_minimal() +
        theme(legend.position = "none")

    print(p)
    
} # end characteristic.curves


# 3/9/24 wanted to plot rmse using *args, passing the entire list into function 

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


# 3/19/24 

plot.rmse.a.diff <- function(metrics_list){
    library(purrr)
  rmse_data <- imap(metrics_list, function(df, name) {
    # Extract method and dentype from the name
    parts <- strsplit(name, "_")[[1]]
    method <- parts[2]
    dentype <- parts[4]
    
    # Select relevant columns and add method, dentype
    mutate(df, method = method, dentype = dentype, item = as.factor(item))
  }) %>% bind_rows()
  
  # Calculate differences in RMSE between true and estimated 'a'
  rmse_data <- rmse_data %>%
    group_by(item, method, dentype) %>%
    summarize(rmse_difference_a = mean(rmse.a), .groups = 'drop')
  
  # Plot
  ggplot(rmse_data, aes(x = item, y = rmse_difference_a, color = method, group = dentype)) +
    geom_line() +
    geom_point() +
    facet_wrap(~dentype, scales = "free_y") +
    labs(title = "RMSE differences for parameter 'a'",
         x = "Item",
         y = "RMSE Difference for 'a'") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom")
} # end plot.rmse





plot.rmse.bias.given.metrics.dfs <- function(metrics.BL, metrics.EM) {
  # Join metrics
  combined.metrics <- full_join(metrics.BL, metrics.EM, by = "item", suffix = c("_BL", "_EM"))
  
  # Convert to longer format
  combined.long <- combined.metrics %>%
    pivot_longer(cols = c("rmse.a_BL", "rmse.a_EM", "rmse.b_BL", "rmse.b_EM", "bias.a_BL", "bias.a_EM", "bias.b_BL", "bias.b_EM"),
                 names_to = "Parameter_Method", values_to = "Value") %>%
    mutate(Method = ifelse(str_detect(Parameter_Method, "_BL"), "BL", "EM"), # extract method
           Type = ifelse(str_detect(Parameter_Method, "rmse"), "RMSE", "Bias"), # extract type
           Parameter = ifelse(str_detect(Parameter_Method, "a_"), "a", "b")) %>% # extract parameter
    select(-Parameter_Method) # remove intermediate column

  # Plot the RMSE and Bias with true values annotated
  p <- ggplot(combined.long, aes(x = interaction(Parameter, Type, sep = " "), y = Value, fill = Method)) +
    geom_col(position = position_dodge(width = 0.8)) +
    facet_wrap(~ item, scales = "free", ncol = 4) +  # Using free scales if the range of values varies significantly
    geom_text(aes(label = sprintf("%s", Parameter)),
              position = position_dodge(width = 0.8), check_overlap = TRUE,
              vjust = 1.5, size = 3, angle = 90) +  # Annotations for 'a' and 'b'
    scale_fill_manual(values = c("blue", "red")) +
    labs(title = "RMSE and Bias for Each Item: BL vs. EM", x = "Parameter and Type", y = "Value") +
    theme_minimal() +
    theme(legend.position = "bottom",
          strip.background = element_blank(),
          strip.text.x = element_text(size = 8)) 

  print(p)
  
}



plot.rmse.bias.hist.metrics.list <- function(metrics_list) {
  # Extract information for each dataframe in the list
  metrics_info <- extract.info.from.list(metrics_list)
  
  # Combine all dataframes into one, ensuring 'Method' and 'Dentype' are included
  combined.metrics <- do.call(rbind, lapply(metrics_info, function(info) {
    cbind(info$data, Method = info$method, Dentype = info$dentype)
  }))

  # Make sure all dataframe elements have a uniform 'item' column
  combined.metrics$item <- factor(combined.metrics$item)

  # Pivot the data to a long format suitable for plotting
  combined.long <- combined.metrics %>%
    pivot_longer(cols = starts_with("rmse") | starts_with("bias"),
                 names_to = c("Metric_Type", "Parameter"),
                 names_sep = "\\.",
                 values_to = "Value") 

  # Add labels for RMSE and BIAS to the plot
  combined.long$Label <- paste(combined.long$Metric_Type, combined.long$Parameter, sep = "\n")

  # Plot the RMSE and Bias for each parameter-method-dentype combination
  p <- ggplot(combined.long, aes(x = interaction(Dentype, Parameter, sep = ": "), y = Value, fill = Method, label = Label)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.7)) +
    geom_text(aes(y = ifelse(Value > 0, Value + 0.02, Value - 0.02)), # Adjust text position based on value
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
          axis.ticks.y = element_line(color = "black", size = 1.2)) +
    scale_y_continuous(limits = c(NA, NA), oob = scales::rescale_none) # Use oob to handle out of bounds

  return(p)
}

