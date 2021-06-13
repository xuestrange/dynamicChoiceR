#clears environment
rm(list = ls())

#import data
data <-
  read.table("./data.txt",
             header = TRUE)

#find the path
I <- which(data[, 3] == 1)

new.data <- data[I, ]

#parameters
beta <- .9
M <- max(new.data[, 4]) + 20
N <- max(new.data[, 1]) + 1


#initial guesses
n.theta1 <- 0
n.theta2 <- 1
n.theta3 <- 1

n.pi0 <- 0.5

#initializing
theta1 <- 0
theta2 <- 0
theta3 <- 0

pi0 <- 0

#number of iterations for EM to converge

count <- 0


time <- system.time(while (abs(theta1 - n.theta1) > 0.00000000001 ||
                           abs(theta2 - n.theta2) > 0.00000000001 ||
                           abs(theta3 - n.theta3) > 0.00000000001) {
  theta1 <- n.theta1
  theta2 <- n.theta2
  theta3 <- n.theta3
  
  pi0 <- n.pi0
  
  
  #The EVs and choice probabilities for s=0
  #utility function
  u0 <- matrix(0, M + 1, 2)
  
  for (x in 1:(M + 1)) {
    for (i in 1:2) {
      u0[x, i] <-
        ifelse(i == 2, -theta2, -theta1 * (x - 1))
    }
  }
  
  #expected value function
  V.bar.0 <- numeric(M + 1) + 1
  n.V.bar.0 <- numeric(M + 1)
  #nested fixed-point to calculate EV
  while (max(abs(n.V.bar.0 - V.bar.0)) > .000001)
  {
    V.bar.0 <- n.V.bar.0
    for (x in 1:M) {
      n.V.bar.0[x] <-
        log(exp(u0[x, 1] + beta * V.bar.0[x + 1]) +
              exp(u0[x, 2] + beta * V.bar.0[2]))
    }
    
    n.V.bar.0[M + 1] <-
      log(exp(u0[M + 1, 1] + beta * V.bar.0[M + 1]) +
            exp(u0[M + 1, 2] + beta * V.bar.0[2]))
  }
  
  
  
  #The EVs and choice probabilities for s=1
  #utility function
  u1 <- matrix(0, M + 1, 2)
  
  for (x in 1:(M + 1)) {
    for (i in 1:2) {
      u1[x, i] <-
        ifelse(i == 2, -theta2 - theta3, -theta1 * (x - 1))
    }
  }
  
  #expected value function
  V.bar.1 <- numeric(M + 1) + 1
  n.V.bar.1 <- numeric(M + 1)
  #nested fixed-point to calculate EV
  while (max(abs(n.V.bar.1 - V.bar.1)) > .000001)
  {
    V.bar.1 <- n.V.bar.1
    for (x in 1:M) {
      n.V.bar.1[x] <-
        log(exp(u1[x, 1] + beta * V.bar.1[x + 1]) +
              exp(u1[x, 2] + beta * V.bar.1[2]))
    }
    
    n.V.bar.1[M + 1] <-
      log(exp(u1[M + 1, 1] + beta * V.bar.1[M + 1]) +
            exp(u1[M + 1, 2] + beta * V.bar.1[2]))
  }
  
  #probability of replacement
  prob.r0 <- numeric(M + 1)
  
  for (x in 2:(M)) {
    prob.r0[x] = exp(u0[x, 2] + beta * V.bar.0[2]) /
      (exp(u0[x, 2] + beta * V.bar.0[2])
       +
         exp(u0[x, 1] + beta * V.bar.0[x + 1]))
    
  }
  prob.r0[M + 1] = exp(u0[x, 2] + beta * V.bar.0[2]) /
    (exp(u0[M + 1, 2] + beta * V.bar.0[2])
     +
       exp(u0[M + 1, 1] + beta * V.bar.0[M + 1]))
  #probability of maintenance
  prob.m0 <- 1 - prob.r0
  
  
  #probability of replacement
  prob.r1 <- numeric(M + 1)
  
  for (x in 2:(M)) {
    prob.r1[x] = exp(u1[x, 2] + beta * V.bar.1[2]) /
      (exp(u1[x, 2] + beta * V.bar.1[2])
       +
         exp(u1[x, 1] + beta * V.bar.1[x + 1]))
    
  }
  prob.r1[M + 1] = exp(u1[M + 1, 2] + beta * V.bar.1[2]) /
    (exp(u1[M + 1, 2] + beta * V.bar.1[2])
     +
       exp(u1[M + 1, 1] + beta * V.bar.1[M + 1]))
  #probability of maintenance
  prob.m1 <- 1 - prob.r1
  
  #Expectation step
  rho.0 <- numeric(N)
  
  for (n in 1:N) {
    J <- which(new.data[, 1] == (n - 1))
    data.bus.n <- new.data[J,]
    rho.0[n] <-
      (pi0 * exp(sum(log(prob.r0[1 + data.bus.n[which(data.bus.n[, 5] == 1), 4]]))) *
         exp(sum(log(prob.m0[1 + data.bus.n[which(data.bus.n[, 5] == 0), 4]])))) /
      ((pi0 * exp(sum(log(
        prob.r0[1 + data.bus.n[which(data.bus.n[, 5] == 1), 4]]
      )))) *
        exp(sum(log(prob.m0[1 + data.bus.n[which(data.bus.n[, 5] == 0), 4]]))) +
        ((1 - pi0) * exp(sum(log(
          prob.r1[1 + data.bus.n[which(data.bus.n[, 5] == 1), 4]]
        ))) *
          exp(sum(log(
            prob.m1[1 + data.bus.n[which(data.bus.n[, 5] == 0), 4]]
          )))))
  }
  
  LL <- function(theta) {
    m.theta1 <- theta[1]
    m.theta2 <- theta[2]
    m.theta3 <- theta[3]
    
    #utility function
    m.u0 <- matrix(0, M + 1, 2)
    
    for (x in 1:(M + 1)) {
      for (i in 1:2) {
        m.u0[x, i] <-
          ifelse(i == 2,-m.theta2,-m.theta1 * (x - 1))
      }
    }
    
    
    #expected value function
    m.V.bar.0 <- numeric(M + 1) + 1
    n.m.V.bar.0 <- numeric(M + 1)
    #nested fixed-point to calculate EV
    while (max(abs(n.m.V.bar.0 - m.V.bar.0)) > .000001)
    {
      m.V.bar.0 <- n.m.V.bar.0
      for (x in 1:M) {
        n.m.V.bar.0[x] <-
          max(m.u0[x, 1] + beta * m.V.bar.0[x + 1], m.u0[x, 2] + beta * m.V.bar.0[2]) +
          log(exp(
            m.u0[x, 1] + beta * m.V.bar.0[x + 1] - max(m.u0[x, 1] + beta * m.V.bar.0[x + 1], m.u0[x, 2] + beta * m.V.bar.0[2])
          ) +
            exp(
              m.u0[x, 2] + beta * m.V.bar.0[2] - max(m.u0[x, 1] + beta * m.V.bar.0[x + 1], m.u0[x, 2] + beta * m.V.bar.0[2])
            ))
      }
      
      n.m.V.bar.0[M + 1] <-
        log(exp(m.u0[M + 1, 1] + beta * m.V.bar.0[M + 1]) +
              exp(m.u0[M + 1, 2] + beta * m.V.bar.0[2]))
    }
    
    
    
    #The EVs and choice probabilities for s=1
    #utility function
    m.u1 <- matrix(0, M + 1, 2)
    
    for (x in 1:(M + 1)) {
      for (i in 1:2) {
        m.u1[x, i] <-
          ifelse(i == 2,-m.theta2 - m.theta3,-m.theta1 * (x - 1))
      }
    }
    
    #expected value function
    m.V.bar.1 <- numeric(M + 1) + 1
    n.m.V.bar.1 <- numeric(M + 1)
    #nested fixed-point to calculate EV
    while (max(abs(n.m.V.bar.1 - m.V.bar.1)) > .000001)
    {
      m.V.bar.1 <- n.m.V.bar.1
      for (x in 1:M) {
        n.m.V.bar.1[x] <-
          max(m.u1[x, 1] + beta * m.V.bar.1[x + 1], m.u1[x, 2] + beta * m.V.bar.1[2]) +
          log(exp(
            m.u1[x, 1] + beta * m.V.bar.1[x + 1] - max(m.u1[x, 1] + beta * m.V.bar.1[x + 1], m.u1[x, 2] + beta * m.V.bar.1[2])
          ) +
            exp(
              m.u1[x, 2] + beta * m.V.bar.1[2] - max(m.u1[x, 1] + beta * m.V.bar.1[x + 1], m.u1[x, 2] + beta * m.V.bar.1[2])
            ))
      }
      
      n.m.V.bar.1[M + 1] <-
        log(exp(m.u1[M + 1, 1] + beta * m.V.bar.1[M + 1]) +
              exp(m.u1[M + 1, 2] + beta * m.V.bar.1[2]))
    }
    
    
    
    #probability of replacement
    m.prob.r0 <- numeric(M + 1)
    
    for (x in 2:(M + 1)) {
      m.prob.r0[x] = exp(m.u0[x, 2] + beta * m.V.bar.0[2]) /
        (exp(m.u0[x, 2] + beta * m.V.bar.0[2])
         +
           exp(m.u0[x, 1] + beta * m.V.bar.0[x + 1]))
      
    }
    m.prob.r0[M + 1] = exp(m.u0[M + 1, 2] + beta * m.V.bar.0[2]) /
      (exp(m.u0[M + 1, 2] + beta * m.V.bar.0[2])
       +
         exp(m.u0[M + 1, 1] + beta * m.V.bar.0[M + 1]))
    #probability of maintenance
    m.prob.m0 <- 1 - m.prob.r0
    
    
    #probability of replacement
    m.prob.r1 <- numeric(M + 1)
    
    for (x in 2:(M + 1)) {
      m.prob.r1[x] = exp(m.u1[x, 2] + beta * m.V.bar.1[2]) /
        (exp(m.u1[x, 2] + beta * m.V.bar.1[2])
         +
           exp(m.u1[x, 1] + beta * m.V.bar.1[x + 1]))
      
    }
    m.prob.r1[M + 1] = exp(m.u1[M + 1, 2] + beta * m.V.bar.1[2]) /
      (exp(m.u1[M + 1, 2] + beta * m.V.bar.1[2])
       +
         exp(m.u1[M + 1, 1] + beta * m.V.bar.1[M + 1]))
    #probability of maintenance
    m.prob.m1 <- 1 - m.prob.r1
    
    
    L1 <- numeric(N)
    
    for (nn in 1:N) {
      J <- which(new.data[, 1] == (nn - 1))
      data.bus.N <- new.data[J,]
      L1[nn] <-
        (rho.0[nn]) * (sum(log(m.prob.r0[1 + data.bus.N[which(data.bus.N[, 5] == 1), 4]])) +
                         sum(log(m.prob.m0[1 + data.bus.N[which(data.bus.N[, 5] == 0), 4]]))) +
        (1 - rho.0[nn]) * (sum(log(m.prob.r1[1 + data.bus.N[which(data.bus.N[, 5] == 1), 4]])) +
                             sum(log(m.prob.m1[1 + data.bus.N[which(data.bus.N[, 5] == 0), 4]])))
    }
    - sum(L1)
  }
  
  library(stats4)
  
  res <- optim(c(theta1, theta2, theta3),
               LL,
               hessian = T)
  
  
  n.pi0 <- sum(rho.0) / N
  
  count <- count + 1
  
  pars <- res$par
  
  value <- res$value
  
  print(c(count, pars, value, n.pi0))
  
  
  n.theta1 <- pars[1]
  n.theta2 <- pars[2]
  n.theta3 <- pars[3]
  if (abs(theta1 - n.theta1) > 0.00000000001 ||
      abs(theta2 - n.theta2) > 0.00000000001 ||
      abs(theta3 - n.theta3) > 0.00000000001) {
    vars <- diag(solve(res$hessian))
    
    se.theta1 <- sqrt(vars[1])
    se.theta2 <- sqrt(vars[2])
    se.theta3 <- sqrt(vars[3])
  }
})
