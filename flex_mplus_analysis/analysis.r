# analysis.r

rm(list = ls()) 

library(dplyr) 
library(readr)

# working directory
cwd <- getwd()
data.path <- paste0(cwd, "/flex_mplus_analysis/data", sep = "")

# specify paths
path.flex <- paste0(data.path, "/results_fllexmirt_intercepts_slopes.csv")
path.mplus <- paste0(data.path, "/results_mplus_intercepts_slopes.csv")

# load data for flex/mplus
data.flex <- read.csv(path.flex, header = TRUE, sep = ",")
data.mplus <- read.csv(path.mplus, header = TRUE, sep = ",")

# flex intercept/slope data
data.flex.se.intercepts <- data.flex %>% select(starts_with("intercept_")) 
data.flex.se.slopes <- data.flex %>% select(starts_with("slope_")) 

# flex intercept/slope data
data.mplus.se.intercepts <- data.mplus %>% select(starts_with("intercept_")) 
data.mplus.se.slopes <- data.mplus %>% select(starts_with("slope_")) 

# Calculations 
# calculate SE for each pair (intercept/slope)
calc.se <- function(x){
    R <- length(x)
    mean.x <- mean(x)
    se <- sqrt(sum((x - mean.x)^2) / (R - 1))
    return(se)
}

# calculate avg SE for each pair (intercept/slope)
calc.avg.se <- function(x){
    R <- length(x)
    mean.x <- mean(x)
    se <- sqrt(sum((x - mean.x)^2) / (R - 1))
    # normalize to get average se
    se <- se / sqrt(R)
    return(se)
}

# calculate var for each pair (intercept/slope)
calc.var <- function(x){
    R <- length(x)
    mean.x <- mean(x)
    var <- sum((x - mean.x)^2) / (R - 1)
    return(var)
}

# calculate RMSD for each pair (intercept/slope)
calc.rmsd <- function(x, se.true){
    R <- length(x)
    # use calc.se for each val until R
    se.est <- sapply(x, function(est) calc.se(rep(est, R)))
    rmsd <- sqrt(sum((se.est - se.true)^2) / R)
    return(rmsd)
}

# Central calculation
# calculate the row of either intercept/slope data
calc.col.stats <- function(col){
    # call each calculation function
    se <- calc.se(col)
    avg.se <- calc.avg.se(col)
    var <- calc.var(col)
    rmsd <- calc.rmsd(col, se) 
    return(c(se=se, avg.se=avg.se, var=var, rmsd=rmsd))
}

calc.pair.stats <- function(intercept.col, row.col){
    # calculate 
}

# calculate stats for each column
intercept.flex.df <- as.data.frame(t(apply(data.flex.se.intercepts, 2, calc.col.stats)))
slope.flex.df <- as.data.frame(t(apply(data.flex.se.slopes, 2, calc.col.stats)))

intercept.mplus.df <- as.data.frame(t(apply(data.mplus.se.intercepts, 2, calc.col.stats)))
slope.mplus.df <- as.data.frame(t(apply(data.mplus.se.slopes, 2, calc.col.stats)))

# Rename rows for clarity
colnames(intercept.flex.df) <- c("SE", "Avg.SE", "Variance", "RMSD")
colnames(slope.flex.df) <- c("SE", "Avg.SE", "Variance", "RMSD")
colnames(intercept.mplus.df) <- c("SE", "Avg.SE", "Variance", "RMSD")
colnames(slope.mplus.df) <- c("SE", "Avg.SE", "Variance", "RMSD")

# add row names
rownames(intercept.flex.df) <- paste0("Intercept ", seq_len(nrow(intercept.flex.df)))
rownames(slope.flex.df) <- paste0("Slope ", seq_len(nrow(intercept.flex.df)))
rownames(intercept.mplus.df) <- paste0("Intercept ", seq_len(nrow(intercept.mplus.df)))
rownames(slope.mplus.df) <- paste0("Slope ", seq_len(nrow(slope.mplus.df)))

# write to csv
save.path <- paste0(data.path, "/tables/")
write.csv(intercept.flex.df, file = paste0(save.path, "intercept_flex.csv"), row.names = TRUE)
write.csv(slope.flex.df, file = paste0(save.path, "slope_flex.csv"), row.names = TRUE)
write.csv(intercept.mplus.df, file = paste0(save.path, "intercept_mplus.csv"), row.names = TRUE)
write.csv(slope.mplus.df, file = paste0(save.path, "slope_mplus.csv"), row.names = TRUE)