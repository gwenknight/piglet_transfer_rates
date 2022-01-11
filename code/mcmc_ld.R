### Different MCMC: LaplacesDemon

library(LaplacesDemon)
source("code/piglet_mrsa_functions.R")

# data(demonsnacks)
# N <- nrow(demonsnacks)
# y <- log(demonsnacks$Calories)
# X <- cbind(1, as.matrix(log(demonsnacks[,c(1,4,10)]+1)))
# J <- ncol(X)
# for (j in 2:J) {X[,j] <- CenterScale(X[,j])}
# mon.names <- "LP"
# parm.names <- as.parm.names(list(beta=rep(0,J), sigma=0))
# pos.beta <- grep("beta", parm.names)
# pos.sigma <- grep("sigma", parm.names)
# PGF <- function(Data) {
#   beta <- rnorm(Data$J)
#   sigma <- runif(1)
#   return(c(beta, sigma))
# }
# MyData <- list(J=J, PGF=PGF, X=X, mon.names=mon.names,
#                parm.names=parm.names, pos.beta=pos.beta, pos.sigma=pos.sigma, y=y)
# 
# 
# Model <- function(parm, Data){
#   ### Parameters
#   beta <- parm[Data$pos.beta]
#   sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
#   parm[Data$pos.sigma] <- sigma
#   ### Log-Prior
#   beta.prior <- dnormv(beta, 0, 1000, log=TRUE)
#   sigma.prior <- dhalfcauchy(sigma, 25, log=TRUE)
#   ### Log-Likelihood
#   mu <- tcrossprod(beta, Data$X)
#   LL <- sum(dnorm(Data$y, mu, sigma, log=TRUE))
#   ### Log-Posterior
#   LP <- LL + sum(beta.prior) + sigma.prior
#   Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
#                    yhat=rnorm(length(mu), mu, sigma), parm=parm)
#   return(Modelout)
# }
# 
# Initial.Values <- c(rep(0,J), 1)
# 
# Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#                      Covar=NULL, Iterations=1000, Status=100, Thinning=1,
#                      Algorithm="AFSS", Specs=list(A=500, B=NULL, m=100, n=0, w=1))

###### PIGLET
### Likelihood 
dist_like <- read.csv("data/seen_predicted.csv")[,-1]
PigData$dist_like <- dist_like
data <- read.csv("data/data_to_fit.csv")[,-1] 
data_6 <- data %>% filter(variable %in% c("V2","V5","V6","V7","V8","V10"))
PigData$data_6 <- data_6
PigData$mon.names <- "LP"

Initial.Values = c(mu2 = 0.188878541032485, mu5 = 0.191811127529857, mu6 = 1.06342477430845,
                                mu7 = 0.650240103957843, mu8 = 1.13092685382576, mu10 = 1.36327733057528,
                                gamma2 = 0.299144953920434, gamma5 = 1.10429707130617, gamma6 = 1.15332183481689,
                                gamma7 = 0.586259749844624, gamma8 = 1.42781770425473, gamma10 = 0.669795153200574,
                                f2 = 0.0752560231516298, f5 = 0.149273263311754, f6 = 0.104081458666405,
                                f7 = 0.00298936941316733, f8 = 0.29765944219356, f10 = 0.30328253933773,
                                grow = 1.72361721520587)
PigData$parm <- names(Initial.Values)

Initial.Values = c(mu2 = 0.188878541032485, mu5 = 0.191811127529857, mu6 = 1.06342477430845,
                   mu7 = 0.650240103957843, mu8 = 1.13092685382576, mu10 = 1.36327733057528,
                   gamma2 = 0.299144953920434, gamma5 = 1.10429707130617, gamma6 = 1.15332183481689,
                   gamma7 = 0.586259749844624, gamma8 = 1.42781770425473, gamma10 = 0.669795153200574,
                   f2 = 0.0752560231516298, f5 = 0.149273263311754, f6 = 0.104081458666405,
                   f7 = 0.00298936941316733, f8 = 0.29765944219356, f10 = 0.30328253933773,
                   grow = 1.72361721520587)

Initial.Values = c(mu2 = 0.1, mu5 = 0.1, mu6 = 0.1,
                   mu7 = 0.1, mu8 = 0.1, mu10 = 0.1,
                   gamma2 = 0.1, gamma5 = 0.1, gamma6 = 0.1,
                   gamma7 = 0.1, gamma8 = 0.1, gamma10 = 0.1,
                   f2 = 0, f5 = 0, f6 = 0,
                   f7 = 0, f8 = 0, f10 = 0,
                   grow = 0.1)

PigData$parm <- names(Initial.Values)

Fit <- LaplacesDemon(run_sim_model_ld, Data=PigData, Initial.Values,
                     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
                     Algorithm="AFSS", Specs=list(A=500, B=NULL, m=100, n=0, w=1))

##### No fitness costs and all same 
Initial.Values = c(mu = 0.1, gamma = 0.1, grow = 0.1)

PigData$parm <- names(Initial.Values)

Fit <- LaplacesDemon(run_sim_model_ld, Data=PigData, Initial.Values,
                     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
                     Algorithm="AFSS", Specs=list(A=500, B=NULL, m=100, n=0, w=1))
?LaplacesDemon
