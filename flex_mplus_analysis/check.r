# First, modify how we process the MIRT data. Replace the current MIRT data prep section with:

# Read MIRT data
data.mirt <- read.csv(path.mirt, header = TRUE, sep = ",")

# Verify the structure matches the file
print("First few rows of raw MIRT data:")
print(head(data.mirt))

# Create proper MIRT processed data frame
data.mirt.processed <- data.frame(
    # First 20 rows are intercepts, second 20 are slopes
    MEAN = c(
        data.mirt$intsemeans[data.mirt$Parameter == "Intercept"],  # First 20 intercepts
        data.mirt$slopesemeans[data.mirt$Parameter == "Intercept"] # Then 20 slopes
    ),
    RMSD = c(
        data.mirt$intRMSSD[data.mirt$Parameter == "Intercept"],
        data.mirt$slopeRMSSD[data.mirt$Parameter == "Intercept"]
    ),
    VAR = c(
        data.mirt$intsevar[data.mirt$Parameter == "Intercept"],
        data.mirt$slopesevar[data.mirt$Parameter == "Intercept"]
    )
)

# Set row names to match the pattern
rownames(data.mirt.processed) <- paste0(
    rep(c("intercept_", "slope_"), each=20),
    rep(1:20, 2)
)

# Verify the structure
print("First few rows of processed MIRT data:")
print(head(data.mirt.processed))