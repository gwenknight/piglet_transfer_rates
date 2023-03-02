#### Figures 

###### ****** load ******* ######################################################
# Libraries 
library(zoo)
library(patchwork)
library(tidyverse)
theme_set(theme_bw(base_size = 11))

# Model code
source("code/functions.R")

# Data
data <- read.csv("data/data_to_fit.csv")[,-1] 
data_6 <- data %>% filter(name %in% c("phi6","phi2","p1","p2","p3","p4"))

### Likelihood 
dist_like <- read.csv("data/seen_predicted.csv")[,-1]

### Total number of bugs  time series
totals <- read.csv("data/totals_bug.csv")[,-1] %>% select(time,value, parent, lim) %>% 
  pivot_wider(names_from = lim) %>% mutate(weight = 1/(max - min))


###### ****** RUN ******* ######################################################

### Run with all rates the same for all elements
## All same
gain_all <- 0.08
gamma_all <- 0.009
fitness_all <- 0.01

theta = c(# gain
  mu2 = gain_all, mu5 = gain_all,
  mu6 = gain_all,mu7 = gain_all,mu8 = gain_all,mu10 = gain_all,
  # loss
  gamma2 = gamma_all,gamma5 = gamma_all,
  gamma6 = gamma_all,gamma7 = gamma_all,gamma8 = gamma_all,gamma10 = gamma_all,
  # fitness
  f2 = fitness_all,f5 = fitness_all,
  f6 = fitness_all,f7 = fitness_all,f8 = fitness_all,f10 = fitness_all,
  grow = 0.1)

tsteps = 384

## Run for these parameters
out <- run_sim(tsteps, theta)
max(out$P_all$time)
out$error

out$P_all %>% filter(time == 49) %>% arrange(desc(freq))
out$Q_all %>% filter(time == 49) %>% arrange(desc(freq))

if(max(out$P_all$time)==tsteps){ # if get to end 
  #### element prevalence from model 
  model_outputp <- out$prev_predict %>% filter(name %in% c("phi6","phi2","p1","p2","p3","p4"))
  model_outputp$prev <- round(model_outputp$prev,2)
  
  # Check what distribution of n_colonies at this prevalence in the model pig
  distributs <- left_join(model_outputp, dist_like, by = "prev") %>% select(parent, time, name, prob_all, n_colonies_prev)
  # e.g. to check 
  #distributs %>% filter(name == "p1", parent == 1, n_colonies_prev == 0.9900) %>% summarise(sum(prob_all))
  
  # lookup the probability from this distribution for the data
  likelihood_lookup_elements <- left_join(data_6, distributs, by = c("parent","time","name","n_colonies_prev")) %>% summarise(sum(log(prob_all)))
  
  #### total bugs output from model 
  model_outputt <- out$totl_predict
  
  likelihood_lookup_totals <- left_join(model_outputt, totals, by = c("time", "parent")) %>% rowwise() %>% mutate(val_in = as.numeric(between(total,min,max))) %>%
    as.data.frame() %>% mutate(likelihood = weight * val_in) %>% summarise(sum(log(likelihood)))
  
  #### Compare to data 
  compare_dat <- likelihood_lookup_elements + likelihood_lookup_totals
  
  compare_dat
  # compare_dat <- as.numeric(left_join(model_outputp, data_6,  by = c("time", "name", "prev", "parent")) %>% ungroup() %>% 
  #                             summarise(sum(log(weighted_prob_all))))  + # log likelihood for elements
  #   as.numeric(left_join(model_outputt, totals, by = c("time", "parent")) %>% rowwise() %>% mutate(val_in = as.numeric(between(total,min,max))) %>%
  #                as.data.frame() %>% mutate(likelihood = weight * val_in) %>% summarise(sum(log(likelihood)))) # log likelhood for totals
}else{compare_dat <- as.numeric(-Inf)}
# return log likelihood
compare_dat
likelihood_lookup_elements
likelihood_lookup_totals


model_outputp <- rename(model_outputp, parent_strain = parent)
pigg_elements <- read.csv("data/pigg_elements.csv")[,-1]

g1 <- ggplot(pigg_elements %>% filter(name %in% c("phi6","phi2","p1","p2","p3","p4")), 
             aes(x=time, y = sum_prop, group = interaction(name, pig))) + 
  geom_line(aes(col = name, linetype = factor(pig)),size = 1.5, alpha = 0.4) + 
  geom_line(data = model_outputp, aes(x = time, y = prev, group = interaction(parent_strain,name))) + 
  geom_point(aes(col = name),size = 1.5) + 
  facet_wrap(name~parent_strain, ncol = 2)  


totalsp <- read.csv("data/totals_bug.csv")[,-1]
g2 <- ggplot(totalsp, aes(x=time, y = value, group = interaction(time,name,parent))) +
  geom_point(aes(col = factor(parent))) + 
  geom_line(aes(group = interaction(lim, parent), col = factor(parent)), lty = "dashed") + 
  scale_y_log10() + 
  geom_line(data = model_outputt, aes(x=time, y = total, group = parent, col = factor(parent)))

g1 / g2


### Run with all rates the same for plasmids / phage elements
## All same
gain_plasmids <- 0.008
loss_plasmids <- 0.2
fitness_plasmids <- 0.1

gain_phage <- 0.1
loss_phage <- 0.009
fitness_phage <- 0.0001

theta_phage_plasmid = c(# gain
  mu2 = gain_phage, mu5 = gain_phage,
  mu6 = gain_plasmids,mu7 = gain_plasmids,mu8 = gain_plasmids,mu10 = gain_plasmids,
  # loss
  gamma2 = loss_phage,gamma5 = loss_phage,
  gamma6 = loss_plasmids,gamma7 = loss_plasmids,gamma8 = loss_plasmids,gamma10 = loss_plasmids,
  # fitness
  f2 = fitness_phage,f5 = fitness_phage,
  f6 = fitness_plasmids,f7 = fitness_plasmids,f8 = fitness_plasmids,f10 = fitness_plasmids,
  grow = 0.1)

tsteps = 384

## Run for these parameters
out <- run_sim(tsteps, theta_phage_plasmid)
max(out$P_all$time)
out$error

#ggplot(out$P_all, aes(x=time, y = freq, group = label)) + geom_line(aes(col = factor(label)))
#ggplot(out$Q_all, aes(x=time, y = freq, group = label)) + geom_line(aes(col = factor(label)))

out$P_all %>% filter(time == 49) %>% arrange(desc(freq))
out$Q_all %>% filter(time == 49) %>% arrange(desc(freq))

if(max(out$P_all$time)==tsteps){ # if get to end 
  #### element prevalence from model 
  model_outputp <- out$prev_predict %>% filter(name %in% c("phi6","phi2","p1","p2","p3","p4"))
  model_outputp$prev <- round(model_outputp$prev,2)
  
  # Check what distribution of n_colonies at this prevalence in the model pig
  distributs <- left_join(model_outputp, dist_like, by = "prev") %>% select(parent, time, name, prob_all, n_colonies_prev)
  # e.g. to check 
  #distributs %>% filter(name == "p1", parent == 1, n_colonies_prev == 0.9900) %>% summarise(sum(prob_all))
  
  # lookup the probability from this distribution for the data
  likelihood_lookup_elements <- left_join(data_6, distributs, by = c("parent","time","name","n_colonies_prev")) %>% summarise(sum(log(prob_all)))
  
  #### total bugs output from model 
  model_outputt <- out$totl_predict
  
  likelihood_lookup_totals <- left_join(model_outputt, totals, by = c("time", "parent")) %>% rowwise() %>% mutate(val_in = as.numeric(between(total,min,max))) %>%
    as.data.frame() %>% mutate(likelihood = weight * val_in) %>% summarise(sum(log(likelihood)))
  
  #### Compare to data 
  compare_dat <- likelihood_lookup_elements + likelihood_lookup_totals
  
  compare_dat
  # compare_dat <- as.numeric(left_join(model_outputp, data_6,  by = c("time", "name", "prev", "parent")) %>% ungroup() %>% 
  #                             summarise(sum(log(weighted_prob_all))))  + # log likelihood for elements
  #   as.numeric(left_join(model_outputt, totals, by = c("time", "parent")) %>% rowwise() %>% mutate(val_in = as.numeric(between(total,min,max))) %>%
  #                as.data.frame() %>% mutate(likelihood = weight * val_in) %>% summarise(sum(log(likelihood)))) # log likelhood for totals
}else{compare_dat <- as.numeric(-Inf)}
# return log likelihood
compare_dat
likelihood_lookup_elements
likelihood_lookup_totals

model_outputp <- rename(model_outputp, parent_strain = parent)
pigg_elements <- read.csv("data/pigg_elements.csv")[,-1]

g1 <- ggplot(pigg_elements %>% filter(name %in% c("phi6","phi2","p1","p2","p3","p4")), 
             aes(x=time, y = sum_prop, group = interaction(name, pig))) + 
  geom_line(aes(col = name, linetype = factor(pig)),size = 1.5, alpha = 0.4) + 
  geom_line(data = model_outputp, aes(x = time, y = prev, group = interaction(parent_strain,name))) + 
  geom_point(aes(col = name),size = 1.5) + 
  facet_wrap(name~parent_strain, ncol = 2)  


totalsp <- read.csv("data/totals_bug.csv")[,-1]
g2 <- ggplot(totalsp, aes(x=time, y = value, group = interaction(time,name,parent))) +
  geom_point(aes(col = factor(parent))) + 
  geom_line(aes(group = interaction(lim, parent), col = factor(parent)), lty = "dashed") + 
  scale_y_log10() + 
  geom_line(data = model_outputt, aes(x=time, y = total, group = parent, col = factor(parent)))

g1 / g2
