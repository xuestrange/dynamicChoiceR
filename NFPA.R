# Title     : Nested Fixed Point Algorithm
# Created by: XUE
# Created on: 5/26/2021

library(msm)

lower <- 0
upper <- 15000
mu <- 6000
sigma <- 4000
p.x0 <- ptnorm(5000, mean = mu, sd = sigma, lower = lower, upper = upper)
p.x1 <- ptnorm(10000, mean = mu, sd = sigma, lower = lower, upper = upper) - p.x0
p.x2 <- 1 - p.x0 - p.x1
p <- c(p.x0, p.x1, p.x2)
RC <- 20
theta1.1 <- 0.5
theta1.2 <- 0.01
beta <- 0.75

myopic_costs <- function(S, MF, params) {
    rc <- params[1]
    maintain_cost <- rep(NA, S)
    replace_cost <- rep(NA, S)
    for (s in 1 : S) {
        maintain_cost[s] <- MF(s, params[-1])
        replace_cost[s] <- rc
    }
    return(cbind(maintain_cost, replace_cost))
}

maintain_cost.linear <- function(s, params) {
    return(s * params[1])
}

choice_prob <- function(cost_array) {
    # minus minmum of each row
    cost <- cost_array - apply(cost_array, 1, min)
    util <- exp(-cost)
    pchoice <- util / rowSums(util)
    return(pchoice)
}

contraction_mapping <- function(S, p, MF, params, beta = 0.75, threshold = 1e-6, suppr_output = FALSE) {
    m.st <- matrix(0, S, S)
    c.lp <- length(p)
    for (i in 1 : S) {
        for (j in 1 : c.lp) {
            if ((i + j - 1) < S) {
                m.st[i, i + j - 1] <- p[j]
            }else if ((i + j - 1) == S) {
                m.st[i, i + j - 1] <- sum(p[j : c.lp])
            }
        }
    }
    m.rep <- cbind(1, matrix(0, S, S - 1))
    k <- 0
    EV <- matrix(0, S, 2)
    EV.myopic <- EV.new <- myopic_costs(S, MF, params, p)
    achieved <- TRUE
    while (max(abs(EV.new - EV)) > threshold) {
        EV <- EV.new
        pchoice <- choice_prob(EV)
        expected_cost <- rowSums(pchoice * EV)
        future_utility_maintain <- m.st %*% expected_cost
        future_utility_repl <- m.rep %*% expected_cost
        future_utility <- cbind(future_utility_maintain, future_utility_repl)
        EV.new <- EV.myopic + beta * future_utility
        k <- k + 1
        if (k == 1000) {
            achieved <- FALSE
        }
    }
    if (!suppr_output) {
        if (achieved) {
            cat("========================================")
            cat("convergence achieved in ", k, "iterations")
            cat("========================================")
        }else {
            cat("========================================")
            cat("CM could not converge!")
            cat("========================================")
        }
    }
    return(list(CP_forward = choice_prob(EV.new), CP_myopic = choice_prob(EV.myopic)))
}

dynamicLogit <- function(params, data, S, p, MF) {
    endog <- data$choice
    exog <- data$state
    N <- length(endog)
    S <- max(exdog) * 2
    m.st <- matrix(0, S, N)
    for (s in 0 : (S - 1)) {
        m.st[s + 1,] <- (exog == s) * 1
    }
    m.d <- rbind(t(1 - endog), endog)
    util <- contraction_mapping(S = S, p = p, MF = MF, params = params, beta = 0.75, suppr_output = TRUE)
    pchoice <- util$CP_forword
    logprob <- log(t(m.st) %*% pchoice)
    return(-sum(logprob * t(m.d)))
}

bounds <- c(1e-6, Inf)
npars <- 2
lin_fit <- optim(par = rep(.1, npars),
                 fn = dynamicLogit,
                 method = "L-BFGS-B",
                 lower = bounds[1],
                 upper = bounds[2],
                 data = data,
                 S = S,
                 p = p,
                 MF = maintain_cost.linear,
                 control = list(fnscale = 1)
)