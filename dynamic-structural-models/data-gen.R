#clears environment
rm(list = ls())

#inputs

N <- 1000
Tee <- 200
M <- 6
beta <- .9
theta1 <- 2
theta2 <- 3

#initial guess of EV fnx

EV0 <- matrix(1, M + 1, 2)

EV1 <- matrix(0, M + 1, 2)

#nested fixed-point to calculate EV

while (max(abs(EV0 - EV1)) > .000001)
{
  EV0 <- EV1
  EV1[1, 1] <-
    log(exp(-theta1 * (1) + beta * EV0[2, 1]) + exp(-theta2 + beta * EV0[2, 2]))
  EV1[1, 2] <-
    log(exp(-theta1 * (1) + beta * EV0[2, 1]) + exp(-theta2 + beta * EV0[2, 2]))
  
  for (i in 2:M) {
    EV1[i, 1] <-
      log(exp(-theta1 * (i) + beta * EV0[i + 1, 1]) + exp(-theta2 + beta * EV0[i +
                                                                                 1, 2]))
    EV1[i, 2] <-
      log(exp(-theta1 * (1) + beta * EV0[2, 1]) + exp(-theta2 + beta * EV0[2, 2]))
  }
  
  EV1[M + 1, 1] <-
    log(exp(-theta1 * (M) + beta * EV0[M + 1, 1]) + exp(-theta2 + beta * EV0[M +
                                                                               1, 2]))
  EV1[M + 1, 2] <-
    log(exp(-theta1 * (1) + beta * EV0[2, 1]) + exp(-theta2 + beta * EV0[2, 2]))
}

#dynamic logit

prob <- c(1:(M + 1))

prob[1] <- .5

for (i in 2:M) {
  prob[i] = exp(-theta2 + beta * EV1[2, 2]) /
    ((exp(-theta1 * (i - 1) + beta * EV1[i, 1])) + exp(-theta2 + beta *
                                                         EV1[2, 2]))
}
prob[M + 1] = exp(-theta2 + beta * EV1[2, 2]) /
  ((exp(-theta1 * (M + 1) + beta * EV1[M + 1, 1])) + exp(-theta2 + beta *
                                                           EV1[2, 2]))

#generate data

newdat <- matrix(0, 1, 6)

for (nn in 1:N) {
  dat <- matrix(0, 2 * Tee, 6)
  
  #specify the bus number
  dat[, 1] = nn - 1
  
  #specify the period
  for (l in 0:(Tee - 1)) {
    dat[2 * l + 1, 2] <- l
    dat[2 * l + 2, 2] <- l
  }
  
  
  #specify the action
  for (m in 0:(Tee - 1)) {
    dat[2 * m + 1, 3] <- 0
    dat[2 * m + 2, 3] <- 1
  }
  
  dat[1, 5] <- EV1[1, 1]
  dat[2, 5] <- EV1[1, 2]
  
  rand1 <- as.matrix(runif(1, min = 0, max = 1))
  
  if (rand1 < prob[1]) {
    dat[2, 6] <- 1
  } else{
    dat[1, 6] <- 1
  }
  
  dat[3, 4] <- 1
  dat[4, 4] <- 1
  
  dat[3, 5] <- EV1[2, 1]
  dat[4, 5] <- EV1[2, 2]
  
  rand2 <- as.matrix(runif(1, min = 0, max = 1))
  
  if (rand2 < prob[2]) {
    dat[4, 6] <- 1
  } else{
    dat[3, 6] <- 1
  }
  
  for (n in 2:(Tee - 1)) {
    if (dat[2 * n, 6] == 1) {
      dat[2 * n + 1, 4] <- 1
      dat[2 * n + 2, 4] <- 1
      dat[2 * n + 1, 5] <- EV1[2, 1]
      dat[2 * n + 2, 5] <- EV1[2, 2]
      randd <- as.matrix(runif(1, min = 0, max = 1))
      if (randd < prob[dat[2 * n + 2, 4] + 1]) {
        dat[2 * n + 2, 6] <- 1
      } else{
        dat[2 * n + 1, 6] <- 1
      }
    } else{
      if (dat[2 * n - 1, 4] < M) {
        dat[2 * n + 1, 4] <- dat[2 * n - 1, 4] + 1
        dat[2 * n + 2, 4] <- dat[2 * n, 4] + 1
        dat[2 * n + 1, 5] <- EV1[dat[2 * n + 1, 4] + 1, 1]
        dat[2 * n + 2, 5] <- EV1[dat[2 * n + 2, 4] + 1, 2]
        randd <- as.matrix(runif(1, min = 0, max = 1))
        if (randd < prob[dat[2 * n + 2, 4] + 1]) {
          dat[2 * n + 2, 6] <- 1
        } else{
          dat[2 * n + 1, 6] <- 1
        }
      } else{
        dat[2 * n + 1, 4] <- dat[2 * n - 1, 4]
        dat[2 * n + 2, 4] <- dat[2 * n, 4]
        dat[2 * n + 1, 5] <- EV1[dat[2 * n + 1, 4], 1]
        dat[2 * n + 2, 5] <- EV1[dat[2 * n + 2, 4], 2]
        randd <- as.matrix(runif(1, min = 0, max = 1))
        if (randd < prob[dat[2 * n + 2, 4] + 1]) {
          dat[2 * n + 2, 6] <- 1
        } else{
          dat[2 * n + 1, 6] <- 1
        }
      }
    }
  }
  newdat <- rbind(newdat, dat)
}

#our data is ready!
newdat <- newdat[-1, ]

#the histogram of choice probabilities at each mileage level
barplot(prob,
        ylim = c(0, 1),
        xlab = "Mileage",
        ylab = "Prob of Replacement")


#the observed probability distribution of mileage
hist(
  newdat[, 4],
  freq = FALSE,
  breaks = seq(
    from = -.5,
    to = (M + .5),
    by = 1
  ),
  ylim = c(0, 1),
  main = "Observed PD of Mileage",
  xlab = "Mileage"
)

write.table(newdat,
            "./data.txt",
            sep = "\t",
            row.names = FALSE)
