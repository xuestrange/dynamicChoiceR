#clears environment
rm(list = ls())

#inputs

N <- 100
Tee <- 100
M <- 200
beta <- .9
theta1 <- 1
theta2 <- 10
theta3 <- 5

#The EVs and choice probabilities for s=0
#utility function
u0 <- matrix(0, M + 1, 2)

for (x in 1:(M + 1)) {
  for (i in 1:2) {
    u0[x, i] <-
      ifelse(x == 1, 0, ifelse(i == 2, -theta2, -theta1 * (x - 1)))
  }
}

#expected value function
EV0.0 <- matrix(1, M + 1, 2)
EV1.0 <- matrix(0, M + 1, 2)

#nested fixed-point to calculate EV
while (max(EV0.0 - EV1.0) > .0000001)
{
  EV0.0 <- EV1.0
  EV1.0[1, 1] <-
    log(exp(u0[2, 1] + beta * EV0.0[2, 1]) +
          exp(u0[2, 2] + beta * EV0.0[2, 2]))
  EV1.0[1, 2] <-
    log(exp(u0[2, 1] + beta * EV0.0[2, 1]) +
          exp(u0[2, 2] + beta * EV0.0[2, 2]))
  
  for (x in 2:M) {
    EV1.0[x, 1] <-
      log(exp(u0[x + 1, 1] + beta * EV0.0[x + 1, 1]) +
            exp(u0[x + 1, 2] + beta * EV0.0[x + 1, 2]))
    EV1.0[x, 2] <-
      log(exp(u0[2, 1] + beta * EV0.0[2, 1]) +
            exp(u0[2, 2] + beta * EV0.0[2, 2]))
  }
  
  EV1.0[M + 1, 1] <-
    log(exp(u0[M + 1, 1] + beta * EV0.0[M + 1, 1]) +
          exp(u0[M + 1, 2] + beta * EV0.0[M + 1, 2]))
  EV1.0[M + 1, 2] <-
    log(exp(u0[2, 1] + beta * EV0.0[2, 1]) +
          exp(u0[2, 2] + beta * EV0.0[2, 2]))
}

#probability of replacement
prob.r0 <- numeric(M + 1)


for (x in 2:(M + 1)) {
  prob.r0[x] = exp(u0[2, 2] + beta * EV1.0[2, 2]) /
    ((exp(u0[x, 1] + beta * EV1.0[x, 1])) +
       exp(u0[2, 2] + beta * EV1.0[2, 2]))
}
#probability of maintenance
prob.m0 <- 1 - prob.r0

#The EVs and choice probabilities for s=1
#utility function
u1 <- matrix(0, M + 1, 2)

for (x in 1:(M + 1)) {
  for (i in 1:2) {
    u1[x, i] <-
      ifelse(x == 1, 0, ifelse(i == 2,-theta2 - theta3,-theta1 * (x - 1)))
  }
}

#expected value function
EV0.1 <- matrix(1, M + 1, 2)
EV1.1 <- matrix(0, M + 1, 2)

#nested fixed-point to calculate EV
while (max(EV0.1 - EV1.1) > .0000001)
{
  EV0.1 <- EV1.1
  EV1.1[1, 1] <-
    log(exp(u1[2, 1] + beta * EV0.1[2, 1]) +
          exp(u1[2, 2] + beta * EV0.1[2, 2]))
  EV1.1[1, 2] <-
    log(exp(u1[2, 1] + beta * EV0.1[2, 1]) +
          exp(u1[2, 2] + beta * EV0.1[2, 2]))
  
  for (x in 2:M) {
    EV1.1[x, 1] <-
      log(exp(u1[x + 1, 1] + beta * EV0.1[x + 1, 1]) +
            exp(u1[x + 1, 2] + beta * EV0.1[x + 1, 2]))
    EV1.1[x, 2] <-
      log(exp(u1[2, 1] + beta * EV0.1[2, 1]) +
            exp(u1[2, 2] + beta * EV0.1[2, 2]))
  }
  
  EV1.1[M + 1, 1] <-
    log(exp(u1[M + 1, 1] + beta * EV0.1[M + 1, 1]) +
          exp(u1[M + 1, 2] + beta * EV0.1[M + 1, 2]))
  EV1.1[M + 1, 2] <-
    log(exp(u1[2, 1] + beta * EV0.1[2, 1]) +
          exp(u1[2, 2] + beta * EV0.1[2, 2]))
}

#probability of replacement
prob.r1 <- numeric(M + 1)


for (x in 2:(M + 1)) {
  prob.r1[x] = exp(u1[2, 2] + beta * EV1.1[2, 2]) /
    ((exp(u1[x, 1] + beta * EV1.1[x, 1])) +
       exp(u1[2, 2] + beta * EV1.1[2, 2]))
}
#probability of maintenance
prob.m1 <- 1 - prob.r1


#generate data

new.dat <- matrix(0, 1, 7)

for (n in 1:N) {
  dat <- matrix(0, 2 * Tee, 7)
  
  #specify the bus number
  dat[, 1] = n - 1
  
  #specify the period
  for (l in 0:(Tee - 1)) {
    dat[2 * l + 1, 3] <- l
    dat[2 * l + 2, 3] <- l
  }
  
  
  #specify the action
  
  #bus type randomization
  rand <- as.matrix(runif(1, min = 0, max = 1))
  ifelse(rand < .5, I <- 0, I <- 1)
  
  #specify the bus type
  if (I == 1) {
    #specify the period
    for (j in 1:(2 * Tee)) {
      dat[j, 2] <- 1
    }
  }
  
  for (tee in 0:(Tee - 1)) {
    dat[2 * tee + 1, 4] <- 0
    dat[2 * tee + 2, 4] <- 1
  }
  
  dat[1, 6] <- (1 - I) * EV1.0[1, 1] + I * EV1.1[1, 1]
  dat[2, 6] <- (1 - I) * EV1.0[1, 2] + I * EV1.1[1, 2]
  
  rand1 <- as.matrix(runif(1, min = 0, max = 1))
  
  if (rand1 < (1 - I) * prob.r0[1] + I * prob.r1[1]) {
    dat[2, 7] <- 1
  } else{
    dat[1, 7] <- 1
  }
  
  dat[3, 5] <- 1
  dat[4, 5] <- 1
  
  dat[3, 6] <- (1 - I) * EV1.0[2, 1] + I * EV1.1[2, 1]
  dat[4, 6] <- (1 - I) * EV1.0[2, 2] + I * EV1.1[2, 2]
  
  rand2 <- as.matrix(runif(1, min = 0, max = 1))
  
  if (rand2 < (1 - I) * prob.r0[2] + I * prob.r1[2]) {
    dat[4, 7] <- 1
  } else{
    dat[3, 7] <- 1
  }
  
  for (t in 2:(Tee - 1)) {
    if (dat[2 * t, 7] == 1) {
      dat[2 * t + 1, 5] <- 1
      dat[2 * t + 2, 5] <- 1
      dat[2 * t + 1, 6] <- (1 - I) * EV1.0[2, 1] +
        I * EV1.1[2, 1]
      dat[2 * t + 2, 6] <- (1 - I) * EV1.0[2, 2] +
        I * EV1.1[2, 2]
      randd <- as.matrix(runif(1, min = 0, max = 1))
      if (randd < (1 - I) * prob.r0[dat[2 * t + 2, 5] + 1] +
          I * prob.r1[dat[2 * t + 2, 5] + 1]) {
        dat[2 * t + 2, 7] <- 1
      } else{
        dat[2 * t + 1, 7] <- 1
      }
    } else{
      if (dat[2 * t - 1, 5] < M) {
        dat[2 * t + 1, 5] <- dat[2 * t - 1, 5] + 1
        dat[2 * t + 2, 5] <- dat[2 * t, 5] + 1
        dat[2 * t + 1, 6] <-
          (1 - I) * EV1.0[dat[2 * t + 1, 5] + 1, 1] +
          I * EV1.1[dat[2 * t + 1, 5] + 1, 1]
        dat[2 * t + 2, 6] <-
          (1 - I) * EV1.0[dat[2 * t + 2, 5] + 1, 2] +
          I * EV1.1[dat[2 * t + 2, 5] + 1, 2]
        randd <- as.matrix(runif(1, min = 0, max = 1))
        if (randd < (1 - I) * prob.r0[dat[2 * t + 2, 5] + 1] +
            I * prob.r1[dat[2 * t + 2, 5] + 1]) {
          dat[2 * t + 2, 7] <- 1
        } else{
          dat[2 * t + 1, 7] <- 1
        }
      } else{
        dat[2 * t + 1, 5] <- dat[2 * t - 1, 5]
        dat[2 * t + 2, 5] <- dat[2 * t, 5]
        dat[2 * t + 1, 6] <-
          (1 - I) * EV1.0[dat[2 * t + 1, 5], 1] +
          I * EV1.1[dat[2 * t + 1, 5], 1]
        dat[2 * t + 2, 6] <-
          I * EV1.0[dat[2 * t + 2, 5], 2] +
          I * EV1.1[dat[2 * t + 2, 5], 2]
        randd <- as.matrix(runif(1, min = 0, max = 1))
        if (randd < (1 - I) * prob.r0[dat[2 * t + 2, 5] + 1] +
            I * prob.r1[dat[2 * t + 2, 5] + 1]) {
          dat[2 * t + 2, 7] <- 1
        } else{
          dat[2 * t + 1, 7] <- 1
        }
      }
    }
  }
  new.dat <- rbind(new.dat, dat)
}

#our data is ready!
new.dat <- new.dat[-1, ]


write.table(
  new.dat,
  "~/dropbox/data4.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)
