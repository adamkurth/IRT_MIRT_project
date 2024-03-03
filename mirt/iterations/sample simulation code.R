rm(list=ls())

# True item parameters

item.pars <- expand.grid(b = c(-2.5, -1.25, 0, 1.25, 2.5), a = c(0.5, 1, 1.5, 2.5))
a <- item.pars$a
b <- item.pars$b

n.items <- nrow(item.pars)
n.persons <- 300
a.est <- b.est <- NULL

# Replications

for (r in 1:100){
  
  cat(r, "\n")
  
  # Person trait level, e.g., math skill
  
  theta <- rnorm(n = n.persons, mean = 0, sd = 1)
  
  # Simulate data with Monte Carlo experiments
  
  response.data <- matrix(NA, n.persons, n.items)
  
  for (i in 1:n.persons){
    for (j in 1:n.items){
      p <- 1 / (1 + exp(-(a[j] * (theta[i] - b[j]))))
      #p <- 1 / (1 + exp(-(a[j] * theta[i] + d[j])))
      u <- runif(n = 1, min = 0, max = 1)
      if (u < p){
        response.data[i, j] <- 1
      } else {
        response.data[i, j] <- 0
      }
    }
  }
  
  tmp <- NULL
  for (j in 1:n.items){
    tmp <- c(tmp, sprintf("I%d", j))
  }
  colnames(response.data) <- tmp
  
  # Calibrate item parameters using the mirt package, default settings
  
  mirt.out <- mirt::mirt(data = response.data, model = 1, itemtype = '2PL')
  
  # See the estimated item parameters
  #mirt::coef(mirt.out, simplify=T, IRTpars=T)
  #mirt::coef(mirt.out, simplify=T, IRTpars=F)
  
  # Calculate estimation accuracy -------------------------------------------
  
  a.est <- cbind(a.est, mirt::coef(mirt.out, simplify=T, IRTpars=T)$items[, 1])
  b.est <- cbind(b.est, mirt::coef(mirt.out, simplify=T, IRTpars=T)$items[, 2])
}

# Calculate MSE and bias

MSE.a <- MSE.b <- bias.a <- bias.b <- rep(NA, n.items)

for (i in 1:n.items){
  MSE.a[i] <- mean((a.est[i, ] - a[i]) ^ 2)
  MSE.b[i] <- mean((b.est[i, ] - b[i]) ^ 2)
  bias.a[i] <- mean(a.est[i, ] - a[i])
  bias.b[i] <- mean(b.est[i, ] - b[i])
}

round(data.frame(item.pars$b, bias.b, MSE.b, item.pars$a, bias.a, MSE.a), 3)

