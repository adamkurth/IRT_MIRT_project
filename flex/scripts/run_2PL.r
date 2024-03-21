# navigate to the directory where the flexMIRT program setup is located
setwd("C:/Program Files/flexMIRT 3.6/")
print(paste('Current working directory:', getwd()))


execute <- function(filename, path) {
    # changes to directory of .exe file
    setwd("C:\\Program Files\\flexMIRT 3.6\\")
    print(paste('Current working directory:', getwd()))

    # Define the command to run the program
    full_path <- paste(path, filename, sep = "")
    # example: command <- "winflexmirt -r test.flexmirt"
    command <- paste("winflexmirt", "-r", full_path, sep=" ")
    print(paste('Command:', command))

    # Run the program using the system command
    system(command, wait = FALSE)
    print("Program has been executed.")
}

filename <- "2PL_calibration.flexmirt"
path <- "C:\\Users\\k33go\\Documents\\vscode\\IRT_MIRT_project\\" # path using R for windows

execute(filename, path)
