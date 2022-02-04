##### LHS search with Q faster code
#### LHS search 
# The final scenario allowed all six elements to have their own gain, loss and fitness cost values. 
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
library(patchwork)


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

# parameters for Scenario 4
Initial.Values = c(mu2 = 0.0000265289370464156, mu5 = 0.000000130932829264206, mu6 = 0.01, 
                   mu7 = 0.00000481693463586813, mu8 = 0.01, mu10 = 0.005, 
                   gamma2 = 0.00001, gamma5 = 0.000001, gamma6 = 0.00000000001, 
                   gamma7 = 0.000001, gamma8 = 0.00000001, gamma10 = 0.000000000001, 
                   f2 = -0.3, f5 = -0.3, f6 = 0.9, 
                   f7 = -0.2, f8 = 0.3, f10 = 0.8, 
                   grow = 0.15, rel_fit = 1)
run_sim_logPosterior(Initial.Values)

param_ranges <- as.data.frame(cbind(Initial.Values - Initial.Values/10,
                                    Initial.Values + Initial.Values/10))
colnames(param_ranges) <- c("min","max")
param_ranges <- abs(param_ranges)

# From the ?lhs example page
# transform a Latin hypercube
nsamples <- 5000
nparameters <- nrow(param_ranges)
X <- randomLHS(nsamples, nparameters) # first = number of samples, 1000 maybe? second = number of parameters, here 2
Y <- matrix(0, nrow=nsamples, ncol=nparameters)
# Assume parameters uniformly arranged over the range (could assume normal etc if have evidence...)
for(ii in 1:nparameters){
  print(ii)
  Y[,ii] <- qunif(X[,ii], min = param_ranges[ii,1], max = param_ranges[ii,2])
}
Y[,c(13,14,16)] <- -Y[,c(13,14,16)]
write.csv(Y, "fits/scn4/lhs10/paraset.csv")


#### Parallel

registerDoParallel(numCores)

foreach (ii=1:nsamples) %dopar% {
  #for(ii in 1:nsamples){
  
  #print(c("Samples number: ", ii))
  
  parameters = c(mu2 = Y[ii,1],mu5 = Y[ii,2],mu6 = Y[ii,3],
                 mu7 = Y[ii,4],mu8 = Y[ii,5],mu10 = Y[ii,6],
                 gamma2 = Y[ii,7],gamma5 = Y[ii,8],gamma6 = Y[ii,9],
                 gamma7 = Y[ii,10],gamma8 = Y[ii,11],gamma10 = Y[ii,12],
                 f2 = Y[ii,13],f5 = Y[ii,14],f6 = Y[ii,15],
                 f7 = Y[ii,16],f8 = Y[ii,17],f10 = Y[ii,18],
                 grow = Y[ii,19], 
                 rel_fit = Y[ii,20])
  
  
  ## Time step = 1 hr
  #tsteps = 16 * 24
  
  ## Run
  #out <- run_sim(tsteps, c(mu_here, gamma_here, growth_here, grate_here))
  #out <- run_sim(tsteps, parameters)
  
  out <- run_sim_logPosterior(parameters)
  if(abs(out)<10000){write.csv(out,paste0("fits/scn4/lhs10/",ii,".csv"))} # don't keep infinite values
}

stopImplicitCluster()


##### Which worked? 
setwd("fits/scn4/lhs10/")

which_para_csv10 <- list.files() 
which_para10 <- as.numeric(sub("\\..*", "",which_para_csv10[-length(which_para_csv10)]))

work10 <-list.files(pattern = "*.csv") %>% 
  map_df(~read_csv(.))
max_ll <- which_para10[which.max(work10$x)]

parameters = c(mu2 = Y[max_ll,1],mu5 = Y[max_ll,2],mu6 = Y[max_ll,3],
               mu7 = Y[max_ll,4],mu8 = Y[max_ll,5],mu10 = Y[max_ll,6],
               gamma2 = Y[max_ll,7],gamma5 = Y[max_ll,8],gamma6 = Y[max_ll,9],
               gamma7 = Y[max_ll,10],gamma8 = Y[max_ll,11],gamma10 = Y[max_ll,12],
               f2 = Y[max_ll,13],f5 = Y[max_ll,14],f6 = Y[max_ll,15],
               f7 = Y[max_ll,16],f8 = Y[max_ll,17],f10 = Y[max_ll,18],
               grow = Y[max_ll,19], 
               rel_fit = Y[max_ll,20])

run_sim_logPosterior(parameters)
run_sim_logPosterior(Initial.Values)

# success rate
100 * length(which_para) / nsamples ## 69% Too wide => suggest wider 

# Plot output
setwd(here())
out <- piglet_mrsa_movement(tsteps, parameters, ini$bacteria, ini$difference_list)
plot_circles(ini$bacteria, out$all_results,"scn4_lhs10")


#######**************** Second LHS **********************##############
#######
param_ranges <- as.data.frame(cbind(Initial.Values - Initial.Values/5,
                                    Initial.Values + Initial.Values/5))
colnames(param_ranges) <- c("min","max")
param_ranges <- abs(param_ranges)

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
Y[,c(13,14,16)] <- -Y[,c(13,14,16)]
write.csv(Y, "fits/scn4/lhs5/paraset.csv")


#### Parallel

registerDoParallel(numCores)

foreach (ii=1:nsamples) %dopar% {
  #for(ii in 1:nsamples){
  
  #print(c("Samples number: ", ii))
  
  parameters = c(mu2 = Y[ii,1],mu5 = Y[ii,2],mu6 = Y[ii,3],
                 mu7 = Y[ii,4],mu8 = Y[ii,5],mu10 = Y[ii,6],
                 gamma2 = Y[ii,7],gamma5 = Y[ii,8],gamma6 = Y[ii,9],
                 gamma7 = Y[ii,10],gamma8 = Y[ii,11],gamma10 = Y[ii,12],
                 f2 = Y[ii,13],f5 = Y[ii,14],f6 = Y[ii,15],
                 f7 = Y[ii,16],f8 = Y[ii,17],f10 = Y[ii,18],
                 grow = Y[ii,19], 
                 rel_fit = Y[ii,20])
  
  ## Time step = 1 hr
  #tsteps = 16 * 24
  
  ## Run
  #out <- run_sim(tsteps, c(mu_here, gamma_here, growth_here, grate_here))
  #out <- run_sim(tsteps, parameters)
  
  out <- run_sim_logPosterior(parameters)
  if(abs(out)<10000){write.csv(out,paste0("fits/scn4/lhs5/",ii,".csv"))} # don't keep infinite values
}

stopImplicitCluster()


##### Which worked? 
setwd("fits/scn4/lhs5/")

which_para_csv5 <- list.files() 
which_para5 <- as.numeric(sub("\\..*", "",which_para_csv5))

which_para5 <- as.numeric(sub(" .*","",sub("\\..*", "",which_para_csv5)))[1:2840]
work5 <-list.files(pattern = "*.csv") %>% 
  map_df(~read_csv(.))
max_ll <- which_para5[which.max(work5$x)]

parameters = c(mu2 = Y[max_ll,1],mu5 = Y[max_ll,2],mu6 = Y[max_ll,3],
               mu7 = Y[max_ll,4],mu8 = Y[max_ll,5],mu10 = Y[max_ll,6],
               gamma2 = Y[max_ll,7],gamma5 = Y[max_ll,8],gamma6 = Y[max_ll,9],
               gamma7 = Y[max_ll,10],gamma8 = Y[max_ll,11],gamma10 = Y[max_ll,12],
               f2 = Y[max_ll,13],f5 = Y[max_ll,14],f6 = Y[max_ll,15],
               f7 = Y[max_ll,16],f8 = Y[max_ll,17],f10 = Y[max_ll,18],
               grow = Y[max_ll,19], 
               rel_fit = Y[max_ll,20])

run_sim_logPosterior(parameters)
run_sim_logPosterior(Initial.Values)

# success rate
100 * length(which_para) / nsamples ## 35% Too wide => suggest wider 

# Plot output
setwd(here())
out <- piglet_mrsa_movement(tsteps, parameters, ini$bacteria, ini$difference_list)
plot_circles(ini$bacteria, out$all_results,"scn4_lhs5")

#######**************** Third LHS **********************##############
#######
param_ranges <- as.data.frame(cbind(Initial.Values - Initial.Values/3,
                                    Initial.Values + Initial.Values/3))
colnames(param_ranges) <- c("min","max")
param_ranges <- abs(param_ranges)

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
Y[,c(13,14,16)] <- -Y[,c(13,14,16)]
write.csv(Y, "fits/scn4/lhs3/paraset.csv")


#### Parallel

registerDoParallel(numCores)

foreach (ii=1:nsamples) %dopar% {
  #for(ii in 1:nsamples){
  
  #print(c("Samples number: ", ii))
  
  parameters = c(mu2 = Y[ii,1],mu5 = Y[ii,2],mu6 = Y[ii,3],
                 mu7 = Y[ii,4],mu8 = Y[ii,5],mu10 = Y[ii,6],
                 gamma2 = Y[ii,7],gamma5 = Y[ii,8],gamma6 = Y[ii,9],
                 gamma7 = Y[ii,10],gamma8 = Y[ii,11],gamma10 = Y[ii,12],
                 f2 = Y[ii,13],f5 = Y[ii,14],f6 = Y[ii,15],
                 f7 = Y[ii,16],f8 = Y[ii,17],f10 = Y[ii,18],
                 grow = Y[ii,19], 
                 rel_fit = Y[ii,20])
  
  ## Time step = 1 hr
  #tsteps = 16 * 24
  
  ## Run
  #out <- run_sim(tsteps, c(mu_here, gamma_here, growth_here, grate_here))
  #out <- run_sim(tsteps, parameters)
  
  out <- run_sim_logPosterior(parameters)
  if(abs(out)<10000){write.csv(out,paste0("fits/scn4/lhs3/",ii,".csv"))} # don't keep infinite values
}

stopImplicitCluster()


##### Which worked? 
setwd("fits/scn4/lhs3/")

which_para_csv3 <- list.files() 
which_para3 <- as.numeric(sub("\\..*", "",which_para_csv3))

work3 <-list.files(pattern = "*.csv") %>% 
  map_df(~read_csv(.))
max_ll <- which_para3[which.max(work3$x)]
Y <- read.csv("paraset.csv")[,-1]

parameters = c(mu2 = Y[max_ll,1],mu5 = Y[max_ll,2],mu6 = Y[max_ll,3],
               mu7 = Y[max_ll,4],mu8 = Y[max_ll,5],mu10 = Y[max_ll,6],
               gamma2 = Y[max_ll,7],gamma5 = Y[max_ll,8],gamma6 = Y[max_ll,9],
               gamma7 = Y[max_ll,10],gamma8 = Y[max_ll,11],gamma10 = Y[max_ll,12],
               f2 = Y[max_ll,13],f5 = Y[max_ll,14],f6 = Y[max_ll,15],
               f7 = Y[max_ll,16],f8 = Y[max_ll,17],f10 = Y[max_ll,18],
               grow = Y[max_ll,19], 
               rel_fit = Y[max_ll,20])
run_sim_logPosterior(parameters)
run_sim_logPosterior(Initial.Values)

# success rate
100 * length(which_para) / nsamples ## 35% Too wide => suggest wider 

# Plot output
setwd(here())
out <- piglet_mrsa_movement(tsteps, parameters, ini$bacteria, ini$difference_list)
plot_circles(ini$bacteria, out$all_results,"scn4_lhs3")
plot_time_series(out,"scn4_lhs3_ts")

# Correlation
worked_para <- as.data.frame(cbind(Y[which_para3,],work3[1:length(which_para3),"x"]))
min(worked_para$x, na.rm = TRUE)
quantile(worked_para$x, prob = 0.75, na.rm = TRUE)
max(worked_para$x, na.rm = TRUE)
plot(worked_para$x)
top_para <- worked_para %>% filter(x > -738) %>% select(-x)
write.csv(cov(top_para), "fits/scn4_cov_lhs3.csv")

sd <- as.numeric(unlist(lapply(top_para, sd, 2)))
write.csv(sd, "fits/scn4_sd_lhs3.csv")

#######**************** Fourth LHS **********************##############
#######
param_ranges <- as.data.frame(cbind(Initial.Values/10,
                                    10 * Initial.Values))
colnames(param_ranges) <- c("min","max")
param_ranges <- abs(param_ranges)

# From the ?lhs example page
# transform a Latin hypercube
nsamples <- 100000
nparameters <- nrow(param_ranges)
X <- randomLHS(nsamples, nparameters) # first = number of samples, 1000 maybe? second = number of parameters, here 2
Y <- matrix(0, nrow=nsamples, ncol=nparameters)
# Assume parameters uniformly arranged over the range (could assume normal etc if have evidence...)
for(ii in 1:nparameters){
  Y[,ii] <- qunif(X[,ii], min = param_ranges[ii,1], max = param_ranges[ii,2])
}
Y[,c(13,14,16)] <- -Y[,c(13,14,16)]
write.csv(Y, "fits/scn4/lhs1000/paraset.csv")


#### Parallel

registerDoParallel(numCores)

foreach (ii=1:nsamples) %dopar% {
  #for(ii in 1:nsamples){
  
  #print(c("Samples number: ", ii))
  
  parameters = c(mu2 = Y[ii,1],mu5 = Y[ii,2],mu6 = Y[ii,3],
                 mu7 = Y[ii,4],mu8 = Y[ii,5],mu10 = Y[ii,6],
                 gamma2 = Y[ii,7],gamma5 = Y[ii,8],gamma6 = Y[ii,9],
                 gamma7 = Y[ii,10],gamma8 = Y[ii,11],gamma10 = Y[ii,12],
                 f2 = Y[ii,13],f5 = Y[ii,14],f6 = Y[ii,15],
                 f7 = Y[ii,16],f8 = Y[ii,17],f10 = Y[ii,18],
                 grow = Y[ii,19], 
                 rel_fit = Y[ii,20])
  
  ## Time step = 1 hr
  #tsteps = 16 * 24
  
  ## Run
  #out <- run_sim(tsteps, c(mu_here, gamma_here, growth_here, grate_here))
  #out <- run_sim(tsteps, parameters)
  
  out <- run_sim_logPosterior(parameters)
  if(abs(out)<10000){write.csv(out,paste0("fits/scn4/lhs1000/",ii,".csv"))} # don't keep infinite values
}

stopImplicitCluster()


##### Which worked? 
setwd("fits/scn4/lhs1000/")

which_para_csv1000 <- list.files() 
which_para1000 <- as.numeric(sub("\\..*", "",which_para_csv1000))

work1000 <-list.files(pattern = "*.csv") %>% 
  map_df(~read_csv(.))
max_ll <- which_para1000[which.max(work1000$x)]

parameters = c(mu2 = Y[max_ll,1],mu5 = Y[max_ll,2],mu6 = Y[max_ll,3],
               mu7 = Y[max_ll,4],mu8 = Y[max_ll,5],mu10 = Y[max_ll,6],
               gamma2 = Y[max_ll,7],gamma5 = Y[max_ll,8],gamma6 = Y[max_ll,9],
               gamma7 = Y[max_ll,10],gamma8 = Y[max_ll,11],gamma10 = Y[max_ll,12],
               f2 = Y[max_ll,13],f5 = Y[max_ll,14],f6 = Y[max_ll,15],
               f7 = Y[max_ll,16],f8 = Y[max_ll,17],f10 = Y[max_ll,18],
               grow = Y[max_ll,19], 
               rel_fit = Y[max_ll,20])
run_sim_logPosterior(parameters)
run_sim_logPosterior(Initial.Values)

### All worse likelihood than initial values.... 
## Can easily improve e.g. for P3 to 
# Initial.Values = c(mu2 = 0.00023631663994558, mu5 = 7.96232764582261e-07, mu6 = 0.0362005899163667, 
# mu7 = 3.84139033003573e-05, mu8 = 0.0516027884447939, mu10 = 0.0202872958267813, 
# gamma2 = 2.36264639681026e-05, gamma5 = 1.14560454890954e-06, 
# gamma6 = 3.33363775439206e-10, gamma7 = 6.92065123035721e-07, 
# gamma8 = 5.05865587238981e-06, gamma10 = 5.34202975235121e-12, 
# f2 = -2.10459304661928, f5 = -1.86778149978006, f6 = 0.673853238158728, 
# f7 = -0.334461607733397, f8 = 2.00273114386821, f10 = 0.706988436254512, 
# grow = 0.0436026095794833, rel_fit = 1.34916656665852)

# success rate
100 * length(which_para) / nsamples ## 1%!

# Plot output
setwd(here())
out <- piglet_mrsa_movement(tsteps, parameters, ini$bacteria, ini$difference_list)
plot_circles(ini$bacteria, out$all_results,"scn4_lhs1000")

##### Analyse output
para <- rbind(read_csv("fits/scn4/lhs10/paraset.csv") %>% mutate("level" = 10),
              read_csv("fits/scn4/lhs10/paraset.csv") %>% mutate("level" = 5),
              read_csv("fits/scn4/lhs10/paraset.csv") %>% mutate("level" = 3)) %>%
  mutate(ll = -Inf)

worked = c(which_para10, which_para5 + 5000, which_para3 + 2*5000)


which_para10 <- as.numeric(sub("\\..*", "",which_para_csv10))

work10 <-list.files(pattern = "*.csv") %>% 
  map_df(~read_csv(.))
max_ll <- which_para10[which.max(work10$x)]

