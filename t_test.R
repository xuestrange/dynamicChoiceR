# Title     : TODO
# Objective : TODO
# Created by: XUE
# Created on: 6/14/2021

t_test <- function (alpha, n, fit){
    t_statistic <- fit$par/(sqrt(diag(solve(fit$hessian)))/n)
    t_value <- qt(alpha, n - 2)
    print(paste("==================Hypothesis Test==================="))
    print(paste("alpha level is: ", alpha))
    for (i in 1:length(t_statistic)){
        if(abs(t_statistic[i]) > abs(t_value)){
            print(paste("   ",fit$par[i], ": significant"))
        }else{
            print(paste("   ",fit$par[i], ": insignificant"))
        }
    }

}
