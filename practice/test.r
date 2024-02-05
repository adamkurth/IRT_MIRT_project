a <- c(1, 2, 3, 4)
b <- c(5, 6, 7, 8)

params <- data.frame(a, b)

sim <- function(params){
    print(params$a)
    print(params$b)
    return(df <- data.frame(params$a, params$b))
}


df <- sim(params)
print(df)
