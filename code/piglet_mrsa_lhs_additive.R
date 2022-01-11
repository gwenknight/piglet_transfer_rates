##### LHS search with Q faster code
#### LHS search 

library(data.table)
library(janitor)
library(foreach)
library(doParallel)
library(tidyverse)
library(lhs)
library(tmvtnorm)
library(tidyverse)
library(here)
library(coda)
library(scales)
library(reshape2)
library(dplyr)
library(Rfast)


################### Data 
data <- read.csv("data/data_to_fit.csv")[,-1] 
data_6 <- data %>% filter(variable %in% c("V2","V5","V6","V7","V8","V10"))

### Likelihood 
dist_like <- read.csv("data/seen_predicted.csv")[,-1]

### Total number of bugs  time series
totals <- read.csv("data/totals_bug.csv")[,-1] %>% select(time,value, parent, lim) %>% 
  pivot_wider(names_from = lim) %>% mutate(weight = 1/(max - min))


numCores <- parallel::detectCores() - 2
numCores

### Functions
source("code/piglet_mrsa_functions.R")
source("code/mcmcmh.r") # had to change line172

#Parameter ranges for LHS;

## Theta that works: 
init.theta = c(mu2 = 0.188878541032485, mu5 = 0.191811127529857, mu6 = 1.06342477430845,
               mu7 = 0.650240103957843, mu8 = 1.13092685382576, mu10 = 1.36327733057528,
               gamma2 = 0.299144953920434, gamma5 = 1.10429707130617, gamma6 = 1.15332183481689,
               gamma7 = 0.586259749844624, gamma8 = 1.42781770425473, gamma10 = 0.669795153200574,
               f2 = 0.0752560231516298, f5 = 0.149273263311754, f6 = 0.104081458666405,
               f7 = 0.00298936941316733, f8 = 0.29765944219356, f10 = 0.30328253933773,
               grow = 1.72361721520587)

mu_mean = - diff(range(init.theta[c("mu2","mu5","mu6","mu7","mu8","mu10")]))
gamma_mean =  - diff(range(init.theta[c("gamma2","gamma5","gamma6","gamma7","gamma8","gamma10")]))
f_mean =  - diff(range(init.theta[c("f2","f5","f6","f7","f8","f10")]))

mu_sd =  diff(range(init.theta[c("mu2","mu5","mu6","mu7","mu8","mu10")]))
gamma_sd =  diff(range(init.theta[c("gamma2","gamma5","gamma6","gamma7","gamma8","gamma10")]))
f_sd =  diff(range(init.theta[c("f2","f5","f6","f7","f8","f10")]))


#Parameter ranges for LHS;

param_ranges <- as.data.frame(cbind(c(rep(mu_mean,10),
                                      rep(gamma_mean,10),
                                      rep(f_mean,10),
                                      init.theta["grow"]),
                                    c(rep(mu_sd,10),
                                      rep(gamma_sd,10),
                                      rep(f_sd,10),
                                      init.theta["grow"]/3)))


# From the ?lhs example page
# transform a Latin hypercube
nsamples <- 50000
nparameters <- 31
X <- randomLHS(nsamples, nparameters) # first = number of samples, 1000 maybe? second = number of parameters, here 2
Y <- matrix(0, nrow=nsamples, ncol=nparameters)
# Assume parameters normally arranged over the range (could assume normal etc if have evidence...)
for(ii in 1:nparameters){
  Y[,ii] <- qunif(X[,ii], min = param_ranges[ii,1], max = param_ranges[ii,2])
}
write.csv(Y, "fits/04_01/paraset.csv")


#### Parallel

registerDoParallel(numCores)

foreach (ii=1:nsamples) %dopar% {
  #for(ii in 1:nsamples){
  
  #print(c("Samples number: ", ii))
  
  # gamma_here <- Y[ii,1:10]
  # mu_here <- Y[ii,11:20]
  # growth_here = Y[ii,21:30]
  # grate_here = Y[ii,31]
  
  parameters <- c(mu2 = Y[ii,2],mu5 = Y[ii,5],mu6 = Y[ii,6],
                  mu7 = Y[ii,7],mu8 = Y[ii,8],mu10 = Y[ii,10],
                  gamma2 = Y[ii,12],gamma5 = Y[ii,15],gamma6 = Y[ii,16],
                  gamma7 = Y[ii,17],gamma8 = Y[ii,18],gamma10 = Y[ii,20],
                  f2 = Y[ii,22],f5 = Y[ii,25],f6 = Y[ii,26],
                  f7 = Y[ii,27],f8 = Y[ii,28],f10 = Y[ii,30],
                  grow = Y[ii,31])
  
  ## Time step = 1 hr
  #tsteps = 16 * 24
  
  ## Run
  #out <- run_sim(tsteps, c(mu_here, gamma_here, growth_here, grate_here))
  #out <- run_sim(tsteps, parameters)
  
  out <- run_sim_logPosterior(parameters)
  if(abs(out)<10000){write.csv(out,paste0("fits/22_12/",ii,".csv"))} # don't keep infinite values
  if(ii %% 1000 == 0){write.csv(out,paste0("fits/04_01/not_",ii,".csv"))} # check working
}

stopImplicitCluster()
