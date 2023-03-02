##### LHS search with Q faster code
#### LHS search 
# The third scenario allowed each of the six elements to have the different fitness cost values but the same gain and loss rates

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
library(prodlim)
library(ggforce)
library(ggridges)


################### Data 
data <- read.csv("data/data_to_fit.csv")[,-1] 
data_6 <- data %>% filter(variable %in% c("V2","V5","V6","V7","V8","V10"))

### Likelihood 
dist_like <- read.csv("data/seen_predicted.csv")[,-1]

### Total number of bugs  time series
totals <- read.csv("data/totals_bug.csv")[,-1] %>% select(time,value, parent, lim) %>% 
  pivot_wider(names_from = lim) %>% mutate(weight = 1/(max - min))

### Functions
source("code/piglet_mrsa_functions.R")

numCores <- parallel::detectCores() - 5
numCores

#Parameter ranges for LHS;

# parameters for Scenario 3
Initial.Values = c(mu = 0.01,
                   gamma = 0.00000001,
                   f2 = 0.000001, f5 = 0.000001, f6 = 0.000001, 
                   f7 = 0.000001, f8 = 0.000001, f10 = 0.000001, 
                   grow = 0.17, 
                   rel_fit = 0.99)
run_sim_logPosterior(Initial.Values)

param_ranges <- as.data.frame(cbind(Initial.Values - Initial.Values/10,
                                    Initial.Values + Initial.Values/10))
colnames(param_ranges) <- c("min","max")

# From the ?lhs example page
# transform a Latin hypercube
nsamples <- 5000
nparameters <- nrow(param_ranges)
X <- randomLHS(nsamples, nparameters) # first = number of samples, 1000 maybe? second = number of parameters, here 2
Y <- matrix(0, nrow=nsamples, ncol=nparameters)
# Assume parameters uniformly arranged over the range (could assume normal etc if have evidence...)
for(ii in 1:nparameters){
  Y[,ii] <- qunif(X[,ii], min = param_ranges[ii,1], max = param_ranges[ii,2])
}
write.csv(Y, "fits/scn3/lhs10/paraset.csv")


#### Parallel

registerDoParallel(numCores)

foreach (ii=1:nsamples) %dopar% {
  #for(ii in 1:nsamples){
  
  #print(c("Samples number: ", ii))
  
  parameters = c(mu = Y[ii,1],
                     gamma = Y[ii,2],
                     f2 = Y[ii,3], f5 = Y[ii,4], f6 = Y[ii,5], 
                     f7 = Y[ii,6], f8 = Y[ii,7], f10 = Y[ii,8], 
                     grow = Y[ii,9], 
                     rel_fit = Y[ii,10])
  
  ## Time step = 1 hr
  #tsteps = 16 * 24
  
  ## Run
  #out <- run_sim(tsteps, c(mu_here, gamma_here, growth_here, grate_here))
  #out <- run_sim(tsteps, parameters)
  
  out <- run_sim_logPosterior(parameters)
  if(abs(out)<10000){write.csv(out,paste0("fits/scn3/lhs10/",ii,".csv"))} # don't keep infinite values
}

stopImplicitCluster()


##### Which worked? 
setwd("fits/scn3/lhs10/")

which_para_csv <- list.files() 
which_para <- as.numeric(sub("\\..*", "",which_para_csv))

work <-list.files(pattern = "*.csv") %>% 
  map_df(~read_csv(.))
max_ll <- which_para[which.max(work$x)]



parameters = c(mu = Y[max_ll,1],
               gamma = Y[max_ll,2],
               f2 = Y[max_ll,3], f5 = Y[max_ll,4], f6 = Y[max_ll,5], 
               f7 = Y[max_ll,6], f8 = Y[max_ll,7], f10 = Y[max_ll,8], 
               grow = Y[max_ll,9], 
               rel_fit = Y[max_ll,10])

run_sim_logPosterior(parameters)
run_sim_logPosterior(Initial.Values)

# success rate
100 * length(which_para) / nsamples ## 69% Too wide => suggest wider 

# Plot output
setwd(here())
out <- piglet_mrsa_movement(tsteps, parameters, ini$bacteria, ini$difference_list)
plot_circles(ini$bacteria, out$all_results,"scn3_lhs10")


#######**************** Second LHS **********************##############
#######
param_ranges <- as.data.frame(cbind(Initial.Values - Initial.Values/5,
                                    Initial.Values + Initial.Values/5))
colnames(param_ranges) <- c("min","max")

# From the ?lhs example page
# transform a Latin hypercube
nsamples <- 5000
nparameters <- nrow(param_ranges)
X <- randomLHS(nsamples, nparameters) # first = number of samples, 1000 maybe? second = number of parameters, here 2
Y <- matrix(0, nrow=nsamples, ncol=nparameters)
# Assume parameters uniformly arranged over the range (could assume normal etc if have evidence...)
for(ii in 1:nparameters){
  Y[,ii] <- qunif(X[,ii], min = param_ranges[ii,1], max = param_ranges[ii,2])
}
write.csv(Y, "fits/scn3/lhs5/paraset.csv")


#### Parallel

registerDoParallel(numCores)

foreach (ii=1:nsamples) %dopar% {
  #for(ii in 1:nsamples){
  
  #print(c("Samples number: ", ii))

  parameters = c(mu = Y[ii,1],
                 gamma = Y[ii,2],
                 f2 = Y[ii,3], f5 = Y[ii,4], f6 = Y[ii,5], 
                 f7 = Y[ii,6], f8 = Y[ii,7], f10 = Y[ii,8], 
                 grow = Y[ii,9], 
                 rel_fit = Y[ii,10])
  
  ## Time step = 1 hr
  #tsteps = 16 * 24
  
  ## Run
  #out <- run_sim(tsteps, c(mu_here, gamma_here, growth_here, grate_here))
  #out <- run_sim(tsteps, parameters)
  
  out <- run_sim_logPosterior(parameters)
  if(abs(out)<10000){write.csv(out,paste0("fits/scn3/lhs5/",ii,".csv"))} # don't keep infinite values
}

stopImplicitCluster()


##### Which worked? 
setwd("fits/scn3/lhs5/")

which_para_csv <- list.files() 
which_para <- as.numeric(sub("\\..*", "",which_para_csv))

work <-list.files(pattern = "*.csv") %>% 
  map_df(~read_csv(.))
max_ll <- which_para[which.max(work$x)]

parameters = c(mu = Y[max_ll,1],
               gamma = Y[max_ll,2],
               f2 = Y[max_ll,3], f5 = Y[max_ll,4], f6 = Y[max_ll,5], 
               f7 = Y[max_ll,6], f8 = Y[max_ll,7], f10 = Y[max_ll,8], 
               grow = Y[max_ll,9], 
               rel_fit = Y[max_ll,10])

run_sim_logPosterior(parameters)
run_sim_logPosterior(Initial.Values)

# success rate
100 * length(which_para) / nsamples ## 35% Too wide => suggest wider 

# Plot output
setwd(here())
out <- piglet_mrsa_movement(tsteps, parameters, ini$bacteria, ini$difference_list)
plot_circles(ini$bacteria, out$all_results,"scn3_lhs5")

#######**************** Third LHS **********************##############
#######
param_ranges <- as.data.frame(cbind(Initial.Values - Initial.Values/3,
                                    Initial.Values + Initial.Values/3))
colnames(param_ranges) <- c("min","max")

# From the ?lhs example page
# transform a Latin hypercube
nsamples <- 5000
nparameters <- nrow(param_ranges)
X <- randomLHS(nsamples, nparameters) # first = number of samples, 1000 maybe? second = number of parameters, here 2
Y <- matrix(0, nrow=nsamples, ncol=nparameters)
# Assume parameters uniformly arranged over the range (could assume normal etc if have evidence...)
for(ii in 1:nparameters){
  Y[,ii] <- qunif(X[,ii], min = param_ranges[ii,1], max = param_ranges[ii,2])
}
write.csv(Y, "fits/scn3/lhs3/paraset.csv")


#### Parallel

registerDoParallel(numCores)

foreach (ii=1:nsamples) %dopar% {
  #for(ii in 1:nsamples){
  
  #print(c("Samples number: ", ii))
  
  parameters = c(mu = Y[ii,1],
                 gamma = Y[ii,2],
                 f2 = Y[ii,3], f5 = Y[ii,4], f6 = Y[ii,5], 
                 f7 = Y[ii,6], f8 = Y[ii,7], f10 = Y[ii,8], 
                 grow = Y[ii,9], 
                 rel_fit = Y[ii,10])
  
  ## Time step = 1 hr
  #tsteps = 16 * 24
  
  ## Run
  #out <- run_sim(tsteps, c(mu_here, gamma_here, growth_here, grate_here))
  #out <- run_sim(tsteps, parameters)
  
  out <- run_sim_logPosterior(parameters)
  if(abs(out)<10000){write.csv(out,paste0("fits/scn3/lhs3/",ii,".csv"))} # don't keep infinite values
}

stopImplicitCluster()


##### Which worked? 
setwd("fits/scn3/lhs3/")

which_para_csv <- list.files() 
which_para <- as.numeric(sub("\\..*", "",which_para_csv))

work <-list.files(pattern = "*.csv") %>% 
  map_df(~read_csv(.))
max_ll <- which_para[which.max(work$x)]

parameters = c(mu = Y[max_ll,1],
               gamma = Y[max_ll,2],
               f2 = Y[max_ll,3], f5 = Y[max_ll,4], f6 = Y[max_ll,5], 
               f7 = Y[max_ll,6], f8 = Y[max_ll,7], f10 = Y[max_ll,8], 
               grow = Y[max_ll,9], 
               rel_fit = Y[max_ll,10])
run_sim_logPosterior(parameters)
run_sim_logPosterior(Initial.Values)

# success rate
100 * length(which_para) / nsamples ## 35% Too wide => suggest wider 

# Plot output
setwd(here())
out <- piglet_mrsa_movement(tsteps, parameters, ini$bacteria, ini$difference_list)
plot_circles(ini$bacteria, out$all_results,"scn3_lhs3")

