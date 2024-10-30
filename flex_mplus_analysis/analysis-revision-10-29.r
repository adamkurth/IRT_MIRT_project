# analysis.r
rm(list = ls()) 
library(dplyr) 
library(readr)
library(tidyr)

# working directory
cwd <- getwd()
data.path <- paste0(cwd, "/flex_mplus_analysis/data/", sep = "")

# specify paths
path.flex <- paste0(data.path, "results_fllexmirt_intercepts_slopes.csv")
path.mplus <- paste0(data.path, "results_mplus_intercepts_slopes.csv")
path.mirt <- paste0(data.path, "mirt_processed_MR_copy.csv")

# load data for flex/mplus
data.flex <- read.csv(path.flex, header = TRUE, sep = ",")
data.mplus <- read.csv(path.mplus, header = TRUE, sep = ",")
data.mirt <- read.csv(path.mirt, header = TRUE, sep = ",")

# rename columns
colnames(data.mirt) <- c("Parameter", "Item", "Intercept_Mean", "Intercept_Var", "Intercept_RMSD", "Intercept_Ratio", 
                        "Slope_Mean", "Slope_Var", "Slope_RMSD", "Slope_Ratio")

# calc.mean <- function(x){
#   n <- length(x)
#   mean.x <- sum(x) / n
#   return(mean.x)
# }

# Calculate Variance
calc.var <- function(x){
  R <- length(x)
  mean.x <- mean(x)
  var <- sum((x - mean.x)^2) / (R - 1)  # no sqrt
  return(var)
}

# Calculate Mean
calc.mean <- function(x){
  mean(x)
}

# Calculate RMSD
calc.rmsd <- function(x) {
  R <- length(x)
  mean.x <- mean(x)
  rmsd <- sqrt(sum((x - mean.x)^2) / R)  # RMSD calculation
  return(rmsd)
}

#  main function to process data
process.data <- function(data){
  results <- list()
  for (i in 1:20){
    for (param in c("intercept", "slope")) {
      col <- paste0(param, "_", i)
      col_data <- data[[col]]
      
      results[[col]] <- c(
        VAR = calc.var(col_data),
        RMSD = calc.rmsd(col_data),
        MEAN = calc.mean(col_data)
      )
    }
  }
  data.processed <- as.data.frame(do.call(rbind, results))
  colnames(data.processed) <- c('VAR', 'RMSD', 'MEAN')  # Reorder columns to match MIRT output
  return(data.processed)
}

# reorder data to intercept_1, intercept_2, ..., slope_19, slope_20
# reorder.data <- function(data) {
#   # confirmed that this is working correctly
#   intercepts <- data[seq(1, 39, 2), ] # select every other row starting from 1
#   slopes <- data[seq(2, 40, 2), ]     # select every other row starting from 2
#   return(rbind(intercepts, slopes))   # combine the two dataframes
# }

# Process flex and mplus data
data.flex.processed <- process.data(data.flex)
data.mplus.processed <- process.data(data.mplus)

print(summary(data.flex.processed))
print(summary(data.mplus.processed))
print(summary(data.mirt.processed))

# write to csv
write.csv(data.flex.processed, file = paste0(data.path, "flex_processed.csv"), row.names = TRUE)
write.csv(data.mplus.processed, file = paste0(data.path, "mplus_processed.csv"), row.names = TRUE)

# prep mirt data
data.mirt.processed <- data.frame(
  MEAN = c(data.mirt$Intercept_Mean, data.mirt$Slope_Mean),
  RMSD = c(data.mirt$Intercept_RMSD, data.mirt$Slope_RMSD),
  VAR = c(data.mirt$Intercept_Var, data.mirt$Slope_Var)
)

# make rownames uniform
rownames(data.flex.processed) <- paste0(rep(c("intercept_", "slope_"), each=20), rep(1:20, 2))
rownames(data.mplus.processed) <- paste0(rep(c("intercept_", "slope_"), each=20), rep(1:20, 2))
rownames(data.mirt.processed) <- paste0(rep(c("intercept_", "slope_"), each=20), rep(1:20, 2))

# CLARIFICATION:
# should create dataframes with 40 rows and 4 columns
# data.processed should look like this:
# 
#               SE    RMSD  MEAN  VAR
# intercept_1  0.1   0.2    0.3   0.4
# intercept_2  0.1   0.2    0.3   0.4
# ...
# slope_19     0.1   0.2    0.3   0.4
# slope_20     0.1   0.2    0.3   0.4 

# safely calculate ratio of two vectors
safe_ratio <- function(numerator, denominator) {
  if (length(numerator) != length(denominator)) {
    warning("Lengths of numerator and denominator do not match")
    return(rep(NA, max(length(numerator), length(denominator))))
  }
  ifelse(denominator == 0, NA, numerator / denominator)
}

# Calculate ratios
ratios <- tryCatch({
  data.frame(
    Parameter = rep(c("Intercept", "Slope"), 20),
    Item = rep(1:20, each = 2),
    Ratio_MEAN_Flex_MIRT = safe_ratio(data.flex.processed$MEAN, data.mirt.processed$MEAN),
    Ratio_MEAN_Flex_Mplus = safe_ratio(data.flex.processed$MEAN, data.mplus.processed$MEAN),
    Ratio_MEAN_MIRT_Mplus = safe_ratio(data.mirt.processed$MEAN, data.mplus.processed$MEAN),
    Ratio_RMSD_Flex_MIRT = safe_ratio(data.flex.processed$RMSD, data.mirt.processed$RMSD),
    Ratio_RMSD_Flex_Mplus = safe_ratio(data.flex.processed$RMSD, data.mplus.processed$RMSD),
    Ratio_RMSD_MIRT_Mplus = safe_ratio(data.mirt.processed$RMSD, data.mplus.processed$RMSD)
  )
  
}, error = function(e) {
  print(paste("Error in creating ratios data frame:", e$message))
  return(NULL)
})

# group by parameter and item
ratios <- ratios %>% group_by(Parameter) %>% arrange(Parameter, Item)


# # check if ratios were calculated successfully
# if (!is.null(ratios)) {
#   print("Ratios calculated successfully")
#   print(head(ratios))
# } else {
#   print("Failed to calculate ratios")
# }

# 1. After loading the data and before the data processing, add true parameter values:
true.parameters <- data.frame(
    Item = 1:20,
    intercept = c(
        # Add your 20 true intercept values here
        # Replace NA with actual values
        rep(NA, 20)
    ),
    slope = c(
        # Add your 20 true slope values here
        # Replace NA with actual values
        rep(NA, 20)
    )
)

# 2. Modify the final.table creation section to include true values:
final.table <- data.frame(
    Parameter = rep(c("Intercept", "Slope"), 20),
    Item = rep(1:20, each = 2),
    True_Value = c(true.parameters$intercept, true.parameters$slope),  # Add this line
    FlexMIRT_Mean = data.flex.processed$MEAN,
    MIRT_Mean = data.mirt.processed$MEAN,
    Mplus_Mean = data.mplus.processed$MEAN,
    Ratio_MEAN_Flex_MIRT = ratios$Ratio_MEAN_Flex_MIRT,
    Ratio_MEAN_Flex_Mplus = ratios$Ratio_MEAN_Flex_Mplus,
    Ratio_MEAN_MIRT_Mplus = ratios$Ratio_MEAN_MIRT_Mplus,
    FlexMIRT_RMSD = data.flex.processed$RMSD,
    MIRT_RMSD = data.mirt.processed$RMSD,
    Mplus_RMSD = data.mplus.processed$RMSD,
    FlexMIRT_VAR = data.flex.processed$VAR,
    MIRT_VAR = data.mirt.processed$VAR,
    Mplus_VAR = data.mplus.processed$VAR
)
# group by parameter and item
final.table <- final.table %>% group_by(Parameter) %>% arrange(Parameter, Item)
final.table <- as.data.frame(final.table)

# 3. Create separate tables for means, RMSD, and variance
# Add this after the final.table creation:

# Means table
means.table <- data.frame(
    Parameter = final.table$Parameter,
    Item = final.table$Item,
    True_Value = final.table$True_Value,
    FlexMIRT = final.table$FlexMIRT_Mean,
    MIRT = final.table$MIRT_Mean,
    Mplus = final.table$Mplus_Mean,
    Ratio_Flex_MIRT = final.table$Ratio_MEAN_Flex_MIRT,
    Ratio_Flex_Mplus = final.table$Ratio_MEAN_Flex_Mplus,
    Ratio_MIRT_Mplus = final.table$Ratio_MEAN_MIRT_Mplus
)

# RMSD table
rmsd.table <- data.frame(
    Parameter = final.table$Parameter,
    Item = final.table$Item,
    True_Value = final.table$True_Value,
    FlexMIRT = final.table$FlexMIRT_RMSD,
    MIRT = final.table$MIRT_RMSD,
    Mplus = final.table$Mplus_RMSD,
    Ratio_Flex_MIRT = final.table$FlexMIRT_RMSD / final.table$MIRT_RMSD,
    Ratio_Flex_Mplus = final.table$FlexMIRT_RMSD / final.table$Mplus_RMSD,
    Ratio_MIRT_Mplus = final.table$MIRT_RMSD / final.table$Mplus_RMSD
)

# Variance table
var.table <- data.frame(
    Parameter = final.table$Parameter,
    Item = final.table$Item,
    True_Value = final.table$True_Value,
    FlexMIRT = final.table$FlexMIRT_VAR,
    MIRT = final.table$MIRT_VAR,
    Mplus = final.table$Mplus_VAR,
    Ratio_Flex_MIRT = final.table$FlexMIRT_VAR / final.table$MIRT_VAR,
    Ratio_Flex_Mplus = final.table$FlexMIRT_VAR / final.table$Mplus_VAR,
    Ratio_MIRT_Mplus = final.table$MIRT_VAR / final.table$Mplus_VAR
)

# 4. Add these additional write.csv commands after the existing ones:
write.csv(means.table, file = paste0(data.path, "means_table.csv"), row.names = FALSE)
write.csv(rmsd.table, file = paste0(data.path, "rmsd_table.csv"), row.names = FALSE)
write.csv(var.table, file = paste0(data.path, "var_table.csv"), row.names = FALSE)

# 5. Add overall ratio calculations for each table:
calculate_overall_ratios <- function(table) {
    table %>%
        group_by(Parameter) %>%
        summarize(
            Overall_Ratio_Flex_MIRT = mean(Ratio_Flex_MIRT, na.rm = TRUE),
            Overall_Ratio_Flex_Mplus = mean(Ratio_Flex_Mplus, na.rm = TRUE),
            Overall_Ratio_MIRT_Mplus = mean(Ratio_MIRT_Mplus, na.rm = TRUE)
        )
    table <- table %>% group_by(Parameter) %>% arrange(Parameter, Item)
    table <- as.data.frame(table)
}

# Calculate overall ratios
means.overall <- calculate_overall_ratios(means.table)
rmsd.overall <- calculate_overall_ratios(rmsd.table)
var.overall <- calculate_overall_ratios(var.table)

# Write overall ratios
write.csv(means.overall, file = paste0(data.path, "means_overall_ratios.csv"), row.names = FALSE)
write.csv(rmsd.overall, file = paste0(data.path, "rmsd_overall_ratios.csv"), row.names = FALSE)
write.csv(var.overall, file = paste0(data.path, "variance_overall_ratios.csv"), row.names = FALSE)

# 6. Add these print statements at the end to verify the output:
print("Means Table (first few rows):")
print(head(means.table))
print("\nRMSD Table (first few rows):")
print(head(rmsd.table))
print("\nVariance Table (first few rows):")
print(head(var.table))
print("\nOverall Ratios:")
print("Means:")
print(means.overall)
print("RMSD:")
print(rmsd.overall)
print("Variance:")
print(var.overall)
