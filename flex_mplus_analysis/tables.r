#tables.r
rm(list=ls())
library(xtable)

# packages for table generation
packages <- c("flex", "mplus")
# packages <- c("flex", "mplus", "mirt")

data.flex.intercepts <- read.csv("flex_mplus_analysis/data/tables/intercept_flex.csv", header = TRUE, sep = ",")
data.flex.slopes <- read.csv("flex_mplus_analysis/data/tables/slope_flex.csv", header = TRUE, sep = ",")
data.mplus.intercepts <- read.csv("flex_mplus_analysis/data/tables/intercept_mplus.csv", header = TRUE, sep = ",")
data.mplus.slopes <- read.csv("flex_mplus_analysis/data/tables/slope_mplus.csv", header = TRUE, sep = ",")

create.package.table <- function(data.intercepts, data.slopes, package) {
    # Use the minimum number of rows between intercepts and slopes
    N <- min(nrow(data.intercepts), nrow(data.slopes))
    
    table <- data.frame(
        Parameter = c(rep("Intercept", N), rep("Slope", N)),
        Itemnumber = rep(1:N, 2),
        SE_Mean = c(data.intercepts$Avg.SE[1:N], data.slopes$Avg.SE[1:N]),
        SE_RMSD = c(data.intercepts$RMSD[1:N], data.slopes$RMSD[1:N]),
        SE_Var = c(data.intercepts$Variance[1:N], data.slopes$Variance[1:N]),
        SE_SE = c(data.intercepts$SE[1:N], data.slopes$SE[1:N])
    )

    colnames(table) <- c("Parameter", "Item Number", "SE Mean", "SE RMSD", "SE Var", "SE")
    
    # Create xtable object
    xtable <- xtable(table, digits = c(2, 4, 2, 4, 5, 4, 5))
    
    # Generate LaTeX code
    latex.code <- print(xtable,
        include.rownames = FALSE,
        caption.placement = "top",
        caption = paste(package, "SE Mean, RMSD, Var, and SE for Intercepts and Slopes"),
        label = paste0("tab:", tolower(package), "_se_stats"),
        table.placement = "H",
        print.results = FALSE
    )
    
    # Save LaTeX code to file
    writeLines(latex.code, paste0("flex_mplus_analysis/data/tables/", tolower(package), "_se_stats.tex"))
    
    cat("LaTeX table for", package, "has been generated and saved.\n")
}

# Create tables for each package
create.package.table(data.flex.intercepts, data.flex.slopes, "flex")
create.package.table(data.mplus.intercepts, data.mplus.slopes, "mplus")

cat("All LaTeX tables have been generated and saved.\n")

