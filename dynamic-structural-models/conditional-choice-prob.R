#clears environment
rm(list = ls())

#import data
data <- read.table("~/dropbox/data_test_assg2.txt", header = FALSE)

#parameters
beta <- .9
M <- max(data[, 4])

#find the path
I <- which(data[, 3] == 0)
new.data <- data[-I, ]


#CCP
CCP.r <- numeric(M + 1)
for (x in 1:(M + 1)) {
  CCP.r[x] <-
    length(which(new.data[which(new.data$V4 == (x - 1)), 5] == 1)) /
    length(which(new.data$V4 == (x - 1)))
}


#log-likelihood function
LL <- function(theta1, theta2) {
  #utility function
  u <- matrix(0, M + 1, 2)
  
  for (x in 1:(M + 1)) {
    for (i in 1:2) {
      u[x, i] <-
        ifelse(x == 1, 0, ifelse(i == 2,-theta2,-theta1 * (x - 1)))
    }
  }
  
  v <- matrix(0, M + 1, 2)
  
  for (x in 1:(M + 1)) {
    for (i in 1:2) {
      v[x, i] <-
        ifelse(i == 2,
               u[x, i] + beta * (u[2, 2] - log(CCP.r[2])),
               ifelse(i == 1 &&
                        x < M, u[x, i] + beta * (u[2, 2] - log(CCP.r[x + 1])),
                      u[x, i] + beta * (u[2, 2] - log(CCP.r[x]))))
    }
  }
  
  #probability of replacement
  prob.r <- numeric(M + 1)
  
  for (x in 2:(M + 1)) {
    prob.r[x] <- exp(v[x, 2]) / (exp(v[x, 2]) + exp(v[x, 1]))
  }
  
  #probability of maintenance
  prob.m <- 1 - prob.r
  
  I.r <- which(new.data[, 5] == 1)
  I.m <- which(new.data[, 5] == 0)
  - sum(log(prob.m[new.data[I.m, 4] + 1])) - sum(log(prob.r[new.data[I.r, 4] + 1]))
}


library(stats4)

#max likelihood estimation
system.time(result <- mle(LL,
              start = list(theta1 =0, theta2 = 0)))

