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
library(ggforce)
library(ggridges)
library(patchwork)
theme_set(theme_bw(base_size = 11))


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

# parameters for Scenario 2
Initial.Values = c(mu_phage = 0.01, mu_plasmid = 0.01, 
                   gamma_phage = 0.00000001,gamma_plasmid = 0.00000001,
                   f_phage = 0.000001, f_plasmid = 0.00000001,
                   grow = 0.17, 
                   rel_fit = 0.99)
run_sim_logPosterior(Initial.Values) # check 


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
write.csv(Y, "fits/scn2/lhs10/paraset.csv")


#### Parallel

registerDoParallel(numCores)

foreach (ii=1:nsamples) %dopar% {
  #for(ii in 1:nsamples){
  
  #print(c("Samples number: ", ii))
  
  
  parameters <-  c(mu_phage = Y[ii,1],
                   mu_plasmid = Y[ii,2],
                   gamma_phage = Y[ii,3],
                   gamma_plasmid = Y[ii,4],
                   f_phage = Y[ii,5], 
                   f_plasmid = Y[ii,6], 
                   grow = Y[ii,7], 
                   rel_fit = Y[ii,8])
  
  ## Time step = 1 hr
  #tsteps = 16 * 24
  
  ## Run
  #out <- run_sim(tsteps, c(mu_here, gamma_here, growth_here, grate_here))
  #out <- run_sim(tsteps, parameters)
  
  out <- run_sim_logPosterior(parameters)
  if(abs(out)<10000){write.csv(out,paste0("fits/scn2/lhs10/",ii,".csv"))} # don't keep infinite values
}

stopImplicitCluster()


##### Which worked? 
setwd("fits/scn2/lhs10/")

which_para_csv10 <- list.files() 
which_para10 <- as.numeric(sub("\\..*", "",which_para_csv10[-length(which_para_csv10)]))

work10 <-list.files(pattern = "*.csv") %>% 
  map_df(~read_csv(.))
max_ll <- which_para10[which.max(work10$x)]

parameters <-  c(mu_phage = Y[max_ll,1],
                 mu_plasmid = Y[max_ll,2],
                 gamma_phage = Y[max_ll,3],
                 gamma_plasmid = Y[max_ll,4],
                 f_phage = Y[max_ll,5], 
                 f_plasmid = Y[max_ll,6], 
                 grow = Y[max_ll,7], 
                 rel_fit = Y[max_ll,8])

run_sim_logPosterior(parameters)
run_sim_logPosterior(Initial.Values)

# success rate
100 * length(which_para) / nsamples ## 

# Plot output
setwd(here())
out <- piglet_mrsa_movement(tsteps, parameters, ini$bacteria, ini$difference_list)
plot_circles(ini$bacteria, out$all_results,"scn2_lhs10")


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
write.csv(Y, "fits/scn2/lhs5/paraset.csv")


#### Parallel

registerDoParallel(numCores)

foreach (ii=1:nsamples) %dopar% {
  #for(ii in 1:nsamples){
  
  #print(c("Samples number: ", ii))

  parameters <-  c(mu_phage = Y[ii,1],
                   mu_plasmid = Y[ii,2],
                   gamma_phage = Y[ii,3],
                   gamma_plasmid = Y[ii,4],
                   f_phage = Y[ii,5], 
                   f_plasmid = Y[ii,6], 
                   grow = Y[ii,7], 
                   rel_fit = Y[ii,8])
  
  
  ## Time step = 1 hr
  #tsteps = 16 * 24
  
  ## Run
  #out <- run_sim(tsteps, c(mu_here, gamma_here, growth_here, grate_here))
  #out <- run_sim(tsteps, parameters)
  
  out <- run_sim_logPosterior(parameters)
  if(abs(out)<10000){write.csv(out,paste0("fits/scn2/lhs5/",ii,".csv"))} # don't keep infinite values
}

stopImplicitCluster()


##### Which worked? 
setwd("fits/scn2/lhs5/")

which_para_csv5 <- list.files() 
which_para5 <- as.numeric(sub("\\..*", "",which_para_csv5[-length(which_para_csv5)]))

work5 <-list.files(pattern = "*.csv") %>% 
  map_df(~read_csv(.))
max_ll <- which_para5[which.max(work5$x)]

parameters <-  c(mu_phage = Y[max_ll,1],
                 mu_plasmid = Y[max_ll,2],
                 gamma_phage = Y[max_ll,3],
                 gamma_plasmid = Y[max_ll,4],
                 f_phage = Y[max_ll,5], 
                 f_plasmid = Y[max_ll,6], 
                 grow = Y[max_ll,7], 
                 rel_fit = Y[max_ll,8])

run_sim_logPosterior(parameters)
run_sim_logPosterior(Initial.Values)

# success rate
100 * length(which_para) / nsamples ## 

# Plot output
setwd(here())
out <- piglet_mrsa_movement(tsteps, parameters, ini$bacteria, ini$difference_list)
plot_circles(ini$bacteria, out$all_results,"scn2_lhs5")

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
write.csv(Y, "fits/scn2/lhs3/paraset.csv")


#### Parallel

registerDoParallel(numCores)

foreach (ii=1:nsamples) %dopar% {
  #for(ii in 1:nsamples){
  
  #print(c("Samples number: ", ii))
  
  parameters <-  c(mu_phage = Y[ii,1],
                   mu_plasmid = Y[ii,2],
                   gamma_phage = Y[ii,3],
                   gamma_plasmid = Y[ii,4],
                   f_phage = Y[ii,5], 
                   f_plasmid = Y[ii,6], 
                   grow = Y[ii,7], 
                   rel_fit = Y[ii,8])
  
  
  ## Time step = 1 hr
  #tsteps = 16 * 24
  
  ## Run
  #out <- run_sim(tsteps, c(mu_here, gamma_here, growth_here, grate_here))
  #out <- run_sim(tsteps, parameters)
  
  out <- run_sim_logPosterior(parameters)
  if(abs(out)<10000){write.csv(out,paste0("fits/scn2/lhs3/",ii,".csv"))} # don't keep infinite values
}

stopImplicitCluster()


##### Which worked? 
setwd("fits/scn2/lhs3/")

which_para_csv3 <- list.files() 
which_para3 <- as.numeric(sub("\\..*", "",which_para_csv3[-length(which_para_csv3)]))

work3 <-list.files(pattern = "*.csv") %>% 
  map_df(~read_csv(.))
max_ll <- which_para3[which.max(work3$x)]

parameters <-  c(mu_phage = Y[max_ll,1],
                 mu_plasmid = Y[max_ll,2],
                 gamma_phage = Y[max_ll,3],
                 gamma_plasmid = Y[max_ll,4],
                 f_phage = Y[max_ll,5], 
                 f_plasmid = Y[max_ll,6], 
                 grow = Y[max_ll,7], 
                 rel_fit = Y[max_ll,8])

run_sim_logPosterior(parameters)
run_sim_logPosterior(Initial.Values)

# success rate
100 * length(which_para) / nsamples ## 

# Plot output
setwd(here())
out <- piglet_mrsa_movement(tsteps, parameters, ini$bacteria, ini$difference_list)
plot_circles(ini$bacteria, out$all_results,"scn2_lhs3")

#######**************** Fourth LHS **********************##############
#######
param_ranges <- as.data.frame(cbind(Initial.Values/10,10 * Initial.Values))
colnames(param_ranges) <- c("min","max")

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
write.csv(Y, "fits/scn2/lhs1000/paraset.csv")


#### Parallel

registerDoParallel(numCores)

foreach (ii=1:nsamples) %dopar% {
  #for(ii in 1:nsamples){
  
  #print(c("Samples number: ", ii))
  
  parameters <-  c(mu_phage = Y[ii,1],
                   mu_plasmid = Y[ii,2],
                   gamma_phage = Y[ii,3],
                   gamma_plasmid = Y[ii,4],
                   f_phage = Y[ii,5], 
                   f_plasmid = Y[ii,6], 
                   grow = Y[ii,7], 
                   rel_fit = Y[ii,8])
  
  
  ## Time step = 1 hr
  #tsteps = 16 * 24
  
  ## Run
  #out <- run_sim(tsteps, c(mu_here, gamma_here, growth_here, grate_here))
  #out <- run_sim(tsteps, parameters)
  
  out <- run_sim_logPosterior(parameters)
  if(abs(out)<10000){write.csv(out,paste0("fits/scn2/lhs1000/",ii,".csv"))} # don't keep infinite values
}

stopImplicitCluster()


##### Which worked? 
setwd("fits/scn2/lhs1000/")

which_para_csv1000 <- list.files() 
which_para1000 <- as.numeric(sub("\\..*", "",which_para_csv1000[-length(which_para_csv1000)]))

work1000 <-list.files(pattern = "*.csv") %>% 
  map_df(~read_csv(.))
max_ll <- which_para1000[which.max(work1000$x)]

parameters <-  c(mu_phage = Y[max_ll,1],
                 mu_plasmid = Y[max_ll,2],
                 gamma_phage = Y[max_ll,3],
                 gamma_plasmid = Y[max_ll,4],
                 f_phage = Y[max_ll,5], 
                 f_plasmid = Y[max_ll,6], 
                 grow = Y[max_ll,7], 
                 rel_fit = Y[max_ll,8])

run_sim_logPosterior(parameters)
run_sim_logPosterior(Initial.Values)

# success rate
100 * length(which_para) / nsamples ## 

# Plot output
setwd(here())
out <- piglet_mrsa_movement(tsteps, parameters, ini$bacteria, ini$difference_list)
plot_circles(ini$bacteria, out$all_results,"scn2_lhs1000")



#####********* Analyse output ********########
para <- rbind(read_csv("fits/scn2/lhs10/paraset.csv") %>% mutate("level" = 10),
              read_csv("fits/scn2/lhs5/paraset.csv") %>% mutate("level" = 5),
              read_csv("fits/scn2/lhs3/paraset.csv") %>% mutate("level" = 3)) %>%
  mutate(ll = -Inf)

para[which_para10,"ll"] <- work10[1:length(which_para10),"x"] # work10 also has paraset read in at end
para[which_para5+5000,"ll"] <- work5[1:length(which_para5),"x"] # 
para[which_para3+2*5000,"ll"] <- work3[1:length(which_para3),"x"] # 

colnames(para) <- c("pa","mu_phage","mu_plasmid","gamma_phage","gamma_plasmid","f_phage","f_plasmid","grow","rel_fit","level","ll")

para_g <- para %>% pivot_longer(cols = "mu_phage":"rel_fit")
ggplot(para_g, aes(x=value, y = ll, group = level)) + geom_point(aes(col = factor(level))) + facet_wrap(~name, scales = "free")

plot(para[which(para$ll > -1080),c(2:9)])
cor(para[which(para$ll > -1070),c(2:9)])

ggplot(para_g[which(para_g$ll > -1070),], aes(x=value, group = level)) + geom_histogram(aes(fill = factor(level))) + facet_wrap(~name, scales = "free")

plot(para_g$ll)
ll_inv <- run_sim_logPosterior(Initial.Values)
ggplot(para_g, aes(x=ll, group = level)) + geom_histogram(binwidth = 10, aes(fill = factor(level))) + geom_vline(xintercept = ll_inv) + 
  geom_density(alpha=.2,aes(fill = factor(level)))

ggplot(para_g, aes(x=ll, group = level)) + geom_density(alpha=.2,aes(fill = factor(level))) + geom_vline(xintercept = ll_inv) 

# Find top 100 
para_gn <- para_g %>% filter(!ll == -Inf)
top_ll <- para_gn[order(para_gn$ll,decreasing = TRUE)[1:1000],]

ggplot(top_ll, aes(x=value,group = name)) + geom_histogram() + facet_wrap(~name, scales = "free")

ggplot(para_g, aes(x=ll, group = level)) + 
#  geom_vline(xintercept = ll_inv) + 
  geom_density(alpha=.2,aes(fill = factor(level))) + 
  geom_density(data = top_ll, aes(fill = factor(level)))

top_ll_wide <- top_ll %>% pivot_wider(names_from = name)
cor(top_ll_wide[,4:8])
plot(top_ll_wide[,4:8])


