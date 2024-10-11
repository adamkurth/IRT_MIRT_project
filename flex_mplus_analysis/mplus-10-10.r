# 10/10 for MPLUS
cwd <- getwd()
data.path <- "/Users/adamkurth/Documents/vscode/code/IRT_MIRT_project/Mplus20varn300estimatesflexlongline.dat"
data <- read.table(data.path, quote =  "/",, comment.char = " ")

intse <- data[,41:60]
sapply(intse, mean)
sapply(intse, var)

slopese <- data[,60:80]
sapply(slopese, mean)
sapply(slopese, var)