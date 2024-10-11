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
colnames(data.mirt) <- c("Parameter", "Item", "Intercept_Mean", "Intercept_Var", "Intercept_RMSD", "Intercept_Ratio", "Slope_Mean", "Slope_Var", "Slope_RMSD", "Slope_Ratio")

# Calculate SE (Standard Error)
calc.se <- function(x){
  R <- length(x)
  mean.x <- mean(x)
  se <- sqrt(sum((x - mean.x)^2) / (R - 1))  # sqrt for SE
  return(se)
}

calc.mean <- function(x){
  n <- length(x)
  mean.x <- sum(x) / n
  return(mean.x)
}

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
        SE = calc.se(col_data),
        VAR = calc.var(col_data),
        RMSD = calc.rmsd(col_data),
        MEAN = calc.mean(col_data)  # Use calc.mean function
      )
    }
  }
  data.processed <- as.data.frame(do.call(rbind, results))
  colnames(data.processed) <- c('SE', 'VAR', 'RMSD', 'MEAN')  # Reorder columns to match MIRT output
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
# data.flex.processed <- reorder.data(data.flex.processed)
# data.mplus.processed <- reorder.data(data.mplus.processed)

print(summary(data.flex.processed))
print(summary(data.mplus.processed))
print(summary(data.mirt.processed))

# write to csv
write.csv(data.flex.processed, file = paste0(data.path, "flex_processed.csv"), row.names = TRUE)
write.csv(data.mplus.processed, file = paste0(data.path, "mplus_processed.csv"), row.names = TRUE)

# prep mirt data
data.mirt.processed <- data.frame(
  SE = c(sqrt(data.mirt$Intercept_Var), sqrt(data.mirt$Slope_Var)),
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
# check if ratios were calculated successfully
if (!is.null(ratios)) {
  print("Ratios calculated successfully")
  print(head(ratios))
} else {
  print("Failed to calculate ratios")
}

final.table <- data.frame(
  Item = rep(1:20, each = 2),
  Parameter = rep(c("Intercept", "Slope"), 20),
  FlexMIRT_Mean = data.flex.processed$MEAN,
  MIRT_Mean = data.mirt.processed$MEAN,
  Mplus_Mean = data.mplus.processed$MEAN,
  Ratio_MEAN_Flex_MIRT = ratios$Ratio_MEAN_Flex_MIRT,
  Ratio_MEAN_Flex_Mplus = ratios$Ratio_MEAN_Flex_Mplus,
  Ratio_MEAN_MIRT_Mplus = ratios$Ratio_MEAN_MIRT_Mplus,
  FlexMIRT_SE = data.flex.processed$SE,
  MIRT_SE = data.mirt.processed$SE,
  Mplus_SE = data.mplus.processed$SE,
  FlexMIRT_RMSD = data.flex.processed$RMSD,
  MIRT_RMSD = data.mirt.processed$RMSD,
  Mplus_RMSD = data.mplus.processed$RMSD,
  FlexMIRT_VAR = data.flex.processed$VAR,
  MIRT_VAR = data.mirt.processed$VAR,
  Mplus_VAR = data.mplus.processed$VAR
)

# Group by Parameter
final.table <- final.table %>%
  group_by(Parameter) %>%
  arrange(Parameter, Item)
  
final.table <- as.data.frame(final.table)

# write final to csv
write.csv(final.table, file = paste0(data.path, "final_table.csv"), row.names = FALSE)
write.csv(ratios, file=paste0(data.path, "ratios.csv"), row.names = FALSE)

head(final.table, 20)
