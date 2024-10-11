rm(list = ls())
# Purpose: this script takes in flexMIRT data from ./data/20varn300result.txt and grabs slope/intercept values for each item
cwd <- getwd()

# read data
path <- "/Users/adamkurth/Documents/vscode/code/IRT_MIRT_project/flex_mplus_analysis/data/20varn300result.txt"
# path <- paste(cwd, "/data/20varn300result.txt", sep = "")
data <- read.table(path, header = FALSE, sep = "\t")

# initialize vectors 
intercepts <- list()
slopes <- list()
 
# process each row in data
for (i in 1:nrow(data)){
  
  # check 47th column for 0, then 48th for 1
  if (data[i, 47] == 0 && data[i, 48] == 1) {
    
    # initialize vectors for this row
    row.intercepts <- numeric()
    row.slopes <- numeric()
    
    # grab first intercept/slope pair
    row.intercepts <- c(row.intercepts, data[i, 49])
    row.slopes <- c(row.slopes, data[i, 50])
    
    print(paste("Row", i, "- First intercept:", row.intercepts[1]))
    print(paste("Row", i, "- First slope:", row.slopes[1]))
    
    # initialize counters
    pair.count <- 1 # counter for intercept/slope pair

    for (j in seq(51, ncol(data), by = 2)){
      if (pair.count == 20) break  # Stop after collecting 20 pairs
      
      # Collect intercept and slope
      row.intercepts <- c(row.intercepts, data[i, j])
      row.slopes <- c(row.slopes, data[i, j+1])
      
      pair.count <- pair.count + 1
    }
        # add collected intercepts/slopes for this row
        intercepts[[i]] <- row.intercepts
        slopes[[i]] <- row.slopes

        print(paste("Row", i, "- Total intercepts:", length(row.intercepts)))
        print(paste("Row", i, "- Total slopes:", length(row.slopes)))
    } else {
        # # if conditions are not met add null to maintain row order
        # intercepts[[i]] <- NULL
        # slopes[[i]] <- NULL
    }
}

# create dataframe
max.pairs <- 20
results.df <- data.frame(matrix(NA, nrow = nrow(data), ncol = max.pairs * 2))
colnames(results.df) <- paste0(rep(c("intercept_", "slope_"), max.pairs), rep(1:max.pairs, each = 2))

# populate the data frame
for (i in 1:nrow(data)) {
    # if intercepts are not null, populate the dataframe
    if (!is.null(intercepts[[i]])) {
        # for each intercept/slope pair, populate the dataframe
        for (j in 1:length(intercepts[[i]])) {
            # populate the dataframe
            results.df[i, paste0("intercept_", j)] <- intercepts[[i]][j]
            results.df[i, paste0("slope_", j)] <- slopes[[i]][j]
        }
    }
}

print(head(results.df))
cat("Number of columns:", ncol(results.df), "\n")
cat("Number of rows:", nrow(results.df), "\n")

# write to file
write.csv(results.df, file = paste(cwd, "/data/results_fllexmirt_intercepts_slopes.csv", sep = ""), row.names = TRUE)
