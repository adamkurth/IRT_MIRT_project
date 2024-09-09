rm(list = ls())
# Purpose: this script takes in Mplus data from ./data/Mplus20varn300estimatesflexlongline.dat and grabs slope/intercept values for each item
cwd <- getwd()

# read data
# path <- paste(cwd, "/data/Mplus20varn300estimatesflexlongline.txt", sep = "")
path <- paste(cwd, "/flex_mplus_analysis/data/Mplus20varn300estimatesflexlongline.txt", sep = "")
data <- read.table(path, header = FALSE, sep = " ")

# initialize vectors 
intercepts <- list()
slopes <- list()

# process each row in data 
for (i in 1:nrow(data)) {
    # initialize vectors for each row
    row.intercepts <- numeric()
    row.slopes <- numeric()

    # grab first intercept/slope pair 
    row.intercepts <- c(row.intercepts, data[i, 41]) #starts at 41 for mplus
    row.slopes <- c(row.slopes, data[i, 61]) #starts at 61 for mplus
    
    print(paste("Row", i, "- First intercept:", row.intercepts[1]))
    print(paste("Row", i, "- First slope:", row.slopes[1]))
    
    # collect all intercepts first
    for (j in 1:19){ # n-1 we only need the remaining 19
        row.intercepts <- c(row.intercepts, data[i, 41 + j])
        row.slopes <- c(row.slopes, data[i, 61 + j])
    }

    # add collected intercepts/slopes for this row
    intercepts[[i]] <- row.intercepts
    slopes[[i]] <- row.slopes

    print(paste("Row", i, "- Total intercepts:", length(row.intercepts)))
    print(paste("Row", i, "- Total slopes:", length(row.slopes)))
}

# create dataframe 
max.pairs <- 20 
results.df <- data.frame(matrix(NA, nrow = nrow(data), ncol = max.pairs * 2))
colnames(results.df) <- paste0(rep(c("intercept_", "slope_"), max.pairs), rep(1:max.pairs, each = 2))

# populate the data frame
for (i in 1:nrow(data)) {
  for (j in 1:length(intercepts[[i]])) {
    results.df[i, paste0("intercept_", j)] <- intercepts[[i]][j]
    results.df[i, paste0("slope_", j)] <- slopes[[i]][j]
  }
}

print(head(results.df))
cat("Number of columns:", ncol(results.df), "\n")
cat("Number of rows:", nrow(results.df), "\n")

# write to file
write.csv(results.df, file = paste(cwd, "/flex_mplus_analysis/data/results_mplus_intercepts_slopes.csv", sep = ""), row.names = TRUE)
