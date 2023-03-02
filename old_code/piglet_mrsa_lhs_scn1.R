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
library(prodlim)


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

numCores <- parallel::detectCores() - 4
numCores

#Parameter ranges for LHS;

# parameters for Scenario 1
Initial.Values = c(mu = 0.01,
                   gamma = 0.00000001,
                   f = 0.000001, 
                   grow = 0.17, 
                   rel_fit = 0.99)


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
write.csv(Y, "fits/scn1/lhs10/paraset.csv")


#### Parallel

registerDoParallel(numCores)

foreach (ii=1:nsamples) %dopar% {
  #for(ii in 1:nsamples){
  
  #print(c("Samples number: ", ii))
  
  # gamma_here <- Y[ii,1:10]
  # mu_here <- Y[ii,11:20]
  # growth_here = Y[ii,21:30]
  # grate_here = Y[ii,31]
  
  parameters <-  c(mu = Y[ii,1],
                   gamma = Y[ii,2],
                   f = Y[ii,3], 
                   grow = Y[ii,4], 
                   rel_fit = Y[ii,5])
  
  ## Time step = 1 hr
  #tsteps = 16 * 24
  
  ## Run
  #out <- run_sim(tsteps, c(mu_here, gamma_here, growth_here, grate_here))
  #out <- run_sim(tsteps, parameters)
  
  out <- run_sim_logPosterior(parameters)
  if(abs(out)<10000){write.csv(out,paste0("fits/scn1/lhs10/",ii,".csv"))} # don't keep infinite values
}

stopImplicitCluster()


##### Which worked? 
setwd("fits/scn1/lhs10/")

which_para_csv <- list.files() 
which_para <- as.numeric(sub("\\..*", "",which_para_csv))

work <-list.files(pattern = "*.csv") %>% 
  map_df(~read_csv(.))
max_ll <- which.max(work$x)

parameters <-  c(mu = Y[max_ll,1],
                 gamma = Y[max_ll,2],
                 f = Y[max_ll,3], 
                 grow = Y[max_ll,4], 
                 rel_fit = Y[max_ll,5])
run_sim_logPosterior(parameters)
run_sim_logPosterior(Initial.Values)

# success rate
100 * length(which_para) / nsamples ## 69% Too narrow => suggest wider 

# Plot output
out <- piglet_mrsa_movement(tsteps, parameters, ini$bacteria, ini$difference_list)
plot_circles(ini$bacteria, out$all_results,"scn1_lhs10")


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
write.csv(Y, "fits/scn1/lhs5/paraset.csv")


#### Parallel

registerDoParallel(numCores)

foreach (ii=1:nsamples) %dopar% {
  #for(ii in 1:nsamples){
  
  #print(c("Samples number: ", ii))
  
  # gamma_here <- Y[ii,1:10]
  # mu_here <- Y[ii,11:20]
  # growth_here = Y[ii,21:30]
  # grate_here = Y[ii,31]
  
  parameters <-  c(mu = Y[ii,1],
                   gamma = Y[ii,2],
                   f = Y[ii,3], 
                   grow = Y[ii,4], 
                   rel_fit = Y[ii,5])
  
  ## Time step = 1 hr
  #tsteps = 16 * 24
  
  ## Run
  #out <- run_sim(tsteps, c(mu_here, gamma_here, growth_here, grate_here))
  #out <- run_sim(tsteps, parameters)
  
  out <- run_sim_logPosterior(parameters)
  if(abs(out)<10000){write.csv(out,paste0("fits/scn1/lhs5/",ii,".csv"))} # don't keep infinite values
}

stopImplicitCluster()


##### Which worked? 
setwd("fits/scn1/lhs5/")

which_para_csv <- list.files() 
which_para <- as.numeric(sub("\\..*", "",which_para_csv))

work <-list.files(pattern = "*.csv") %>% 
  map_df(~read_csv(.))
max_ll <- which.max(work$x)

parameters <-  c(mu = Y[max_ll,1],
                 gamma = Y[max_ll,2],
                 f = Y[max_ll,3], 
                 grow = Y[max_ll,4], 
                 rel_fit = Y[max_ll,5])
run_sim_logPosterior(parameters)
run_sim_logPosterior(Initial.Values)

# success rate
100 * length(which_para) / nsamples ## 35% Too wide => suggest wider 

# Plot output
out <- piglet_mrsa_movement(tsteps, parameters, ini$bacteria, ini$difference_list)
setwd("../..")
plot_circles(ini$bacteria, out$all_results,"scn1_lhs5")

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
write.csv(Y, "fits/scn1/lhs3/paraset.csv")


#### Parallel

registerDoParallel(numCores)

foreach (ii=1:nsamples) %dopar% {
  #for(ii in 1:nsamples){
  
  #print(c("Samples number: ", ii))
  
  # gamma_here <- Y[ii,1:10]
  # mu_here <- Y[ii,11:20]
  # growth_here = Y[ii,21:30]
  # grate_here = Y[ii,31]
  
  parameters <-  c(mu = Y[ii,1],
                   gamma = Y[ii,2],
                   f = Y[ii,3], 
                   grow = Y[ii,4], 
                   rel_fit = Y[ii,5])
  
  ## Time step = 1 hr
  #tsteps = 16 * 24
  
  ## Run
  #out <- run_sim(tsteps, c(mu_here, gamma_here, growth_here, grate_here))
  #out <- run_sim(tsteps, parameters)
  
  out <- run_sim_logPosterior(parameters)
  if(abs(out)<10000){write.csv(out,paste0("fits/scn1/lhs3/",ii,".csv"))} # don't keep infinite values
}

stopImplicitCluster()


##### Which worked? 
setwd("fits/scn1/lhs3/")

which_para_csv <- list.files() 
which_para <- as.numeric(sub("\\..*", "",which_para_csv))

work <-list.files(pattern = "*.csv") %>% 
  map_df(~read_csv(.))
max_ll <- which.max(work$x)

parameters <-  c(mu = Y[max_ll,1],
                 gamma = Y[max_ll,2],
                 f = Y[max_ll,3], 
                 grow = Y[max_ll,4], 
                 rel_fit = Y[max_ll,5])
run_sim_logPosterior(parameters)
run_sim_logPosterior(Initial.Values)

# success rate
100 * length(which_para) / nsamples ## 20% super

# Plot output
out <- piglet_mrsa_movement(tsteps, parameters, ini$bacteria, ini$difference_list)
setwd("../..")
plot_circles(ini$bacteria, out$all_results,"scn1_lhs3")


#####******* Grab all parameters ******#######
setwd("fits/scn1/lhs3/")
which_para_csv <- list.files() 
which_para <- as.numeric(sub("\\..*", "",which_para_csv))

work <-list.files(pattern = "*.csv") %>% 
  map_df(~read_csv(.))
max_ll <- which.max(work$x)
Y <- read.csv("paraset.csv")[,-1]
parameters <-  c(mu = Y[max_ll,1],
                 gamma = Y[max_ll,2],
                 f = Y[max_ll,3], 
                 grow = Y[max_ll,4], 
                 rel_fit = Y[max_ll,5])
run_sim_logPosterior(parameters)
run_sim_logPosterior(Initial.Values)


