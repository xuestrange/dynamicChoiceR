#clears environment
rm(list = ls())

#import data
data <- read.table("./data.txt", header = TRUE)

#find the path
I <- which(data[, 3] == 1)

new.data <- data[I,]

#parameters
beta <- .9
M <- 200

#log-likelihood function
LL <- function(theta1, theta2) {
  #utility function
  u <- matrix(0, M + 1, 2)
  
  for (x in 1:(M + 1)) {
    for (i in 1:2) {
      u[x, i] <-
        ifelse(x == 1, 0, ifelse(i == 2, -theta2, -theta1 * (x - 1)))
    }
  }
  
  #expected value function
  EV0 <- matrix(1, M + 1, 2)
  EV1 <- matrix(0, M + 1, 2)
  
  #nested fixed-point to calculate EV
  while (max(EV0 - EV1) > .00001)
  {
    EV0 <- EV1
    EV1[1, 1] <-
      log(exp(u[2, 1] + beta * EV0[2, 1]) +
            exp(u[2, 2] + beta * EV0[2, 2]))
    EV1[1, 2] <-
      log(exp(u[2, 1] + beta * EV0[2, 1]) +
            exp(u[2, 2] + beta * EV0[2, 2]))
    
    for (x in 2:M) {
      EV1[x, 1] <-
        log(exp(u[x + 1, 1] + beta * EV0[x + 1, 1]) +
              exp(u[x + 1, 2] + beta * EV0[x + 1, 2]))
      EV1[x, 2] <-
        log(exp(u[2, 1] + beta * EV0[2, 1]) +
              exp(u[2, 2] + beta * EV0[2, 2]))
    }
    
    EV1[M + 1, 1] <-
      log(exp(u[M + 1, 1] + beta * EV0[M + 1, 1]) +
            exp(u[M + 1, 2] + beta * EV0[M + 1, 2]))
    EV1[M + 1, 2] <-
      log(exp(u[2, 1] + beta * EV0[2, 1]) +
            exp(u[2, 2] + beta * EV0[2, 2]))
  }
  
  #probability of replacement
  prob.r <- numeric(M + 1)
  
  prob.r[1] <- .5
  for (x in 2:(M + 1)) {
    prob.r[x] = exp(u[2, 2] + beta * EV1[2, 2] -
                      max(u[x, 1] + beta * EV1[x, 1], u[2, 2] + beta * EV1[2, 2])) /
      ((exp(
        u[x, 1] + beta * EV1[x, 1] -
          max(u[x, 1] + beta * EV1[x, 1], u[2, 2] + beta * EV1[2, 2])
      )) +
        exp(u[2, 2] + beta * EV1[2, 2] -
              max(u[x, 1] + beta * EV1[x, 1], u[2, 2] + beta * EV1[2, 2])))
  }
  #probability of maintenance
  prob.m <- 1 - prob.r
  
  I.r <- which(new.data[, 5] == 1)
  I.m <- which(new.data[, 5] == 0)
  - sum(log(prob.m[new.data[I.m, 4] + 1])) - sum(log(prob.r[new.data[I.r, 4] + 1]))
}

library(stats4)

#max likelihood estimation
time.nfp <- system.time(res <- mle(LL,
                                   start = list(theta1 = 1, theta2 = 1)))
