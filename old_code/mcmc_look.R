### Look at output

library(tmvtnorm)
library(tidyverse)
library(here)
library(coda)
library(ggplot2)
library(scales)
library(reshape2)
library(dplyr)
library(Rfast)
library(patchwork)
source("code/piglet_mrsa_functions.R")
source("code/mcmcmh.r") # had to change line172

### CHANGE
filename <- "scn4_a_2022-02-09_15-50-20_GMT"
mcmc.epi_every <- c()
mcmc.epi_every$trace <- read.csv(here::here("fits",paste0(filename,"_","trace",".csv")))[,-1]
#mcmc.epi_every$acceptance.rate <- read.csv(here::here("fits/phage_plasmid",paste0(filename,"_","acceptance_rates",".csv")))
#mcmc.epi_every$covmat.empirical <- read.csv(here::here("fits/phage_plasmid",paste0(filename,"_","covmat",".csv")))



# # # Look at output
mcmc.trace <- mcmc(mcmc.epi_every$trace)
summary(mcmc.trace)
acceptanceRate <- 1 - rejectionRate(mcmc.trace)
acceptanceRate
plot(mcmc.trace)
effectiveSize(mcmc.trace)


mcmc.trace.burned <- burnAndThin(mcmc.trace, burn = 100)
plot(mcmc.trace.burned)

autocorr.plot(mcmc.trace.burned)
#
mcmc.trace.burned.thinned <- burnAndThin(mcmc.trace.burned, thin = 5)
autocorr.plot(mcmc.trace.burned.thinned)

plotESSBurn(mcmc.trace)

tail(mcmc.trace)

#parameters_every <- mcmc.trace[1,1:(ncol(mcmc.trace))]
parameters_every <- mcmc.trace[nrow(mcmc.trace),1:(ncol(mcmc.trace))]

## Run for these parameters
out <- piglet_mrsa_movement(tsteps, parameters_every, ini$bacteria, ini$difference_list)

### Data
data <- read.csv("data/data_to_fit.csv")[,-1]
data_6 <- data %>% filter(variable %in% c("V2","V5","V6","V7","V8","V10"))
### Likelihood
dist_like <- read.csv("data/seen_predicted.csv")[,-1]

### Total number of bugs  time series
totals <- read.csv("data/totals_bug.csv")[,-1] %>% select(time,value, parent, lim) %>%
  pivot_wider(names_from = lim) %>% mutate(weight = 1/(max - min))

pigg_elements <- read.csv("data/pigg_elements.csv")[,-1]

### Likelihood
if(!is.null(out$prev_predict)){
  
  #### element prevalence from model (a)
  #c("SCCmec","phi6","SaPI","Tn916","phi2","p1","p2","p3","phi3","p4")
  model_outputp <- out$prev_predict %>% filter(variable %in% c("V2","V5","V6","V7","V8","V10"))
  model_outputp$prev <- round(model_outputp$value,2)
  
  # Check what distribution of n_colonies at this prevalence in the model pig
  distributs <- left_join(model_outputp, dist_like, by = "prev") %>% select(parent, time, variable, prob_all, n_colonies_prev)
  # e.g. to check 
  #distributs %>% filter(name == "p1", parent == 1, n_colonies_prev == 0.9900) %>% summarise(sum(prob_all))
  
  # lookup the probability from this distribution for the data
  likelihood_lookup_elements <- left_join(data_6, distributs, by = c("parent","time","variable","n_colonies_prev")) %>% 
    mutate(exp_miss = dnorm(prob_all, 1, 0.3)) %>% # experimental measure - 10% around a prob of 1
    summarise(sum(log(exp_miss))) #summarise(sum(log(prob_all)))
  
  #### total bugs output from model (b)
  model_outputt <- out$totl_predict
  
  likelihood_lookup_totals <- left_join(model_outputt, totals, by = c("time", "parent")) %>% rowwise() %>% 
    #mutate(val_in = as.numeric(between(total,min,max))) %>% # instead of 1 / 0 make distance 
    ungroup() %>% 
    mutate(val_in = dnorm(log10(total), mean = (log10(max) - log10(min))/2 + log10(min), sd = (log10(max) - log10(min))/20)) %>% #mean = (max - min)/2 + min, sd = (max - min)/10000)) %>% # instead of 1 / 0 make distance 
    as.data.frame() %>% mutate(likelihood = val_in) %>% summarise(sum(log(likelihood)))
  
  #### ensure certain profiles present (c)
  total_end <- unlist(out$all_results %>% filter(time == tsteps) %>% group_by(parent) %>% summarise(total = sum(value)) %>%select(total))
  ## Need to be present at > 80% and > 20% for parent 1 and parent 2 respectively 
  profile_end <- out$all_results %>% filter(time == tsteps, variable %in% profiles_needed_end) 
  profile_end$variable <- as.numeric(profile_end$variable)
  if(nrow(profile_end) == 3){prof_end <- profile_end$value}else{
    prf_end <- as.data.frame(cbind(profiles_needed_end,c(0,0,0)));colnames(prf_end) <- c("variable","end") 
    p <- left_join(prf_end, profile_end, by = "variable")
    p[which(is.na(p$value)), "value"] <- 0
    prof_end <- p$value
  }
  # If total_end = 0 then prof_end will be 0 too 
  if(total_end[1]>0){ prof_end[1] <-  prof_end[1]/total_end[1]}
  if(total_end[2]>0){ prof_end[2:3] <-  prof_end[2:3]/total_end[2]}
  # Could change sd etc if not close enough 
  likelihood_profile_end <- sum(log(dnorm(prof_end,mean = c(0.8,0.2,0.2), sd = 0.5)))
  
  # Don't make -Inf possible
  # if(nrow(profile_end) == 3 && total_end[2] > 0){
  #   likelihood_profile_end <- profile_end %>% 
  #     mutate(prop = value / c(total_end[1],total_end[2],total_end[2]),
  #            cutoff = pmax(0,prop - 0.05),#c(0.8,0.2,0.2)), # more the better # not using -- too strict
  #            likelihood = log(prop)) %>% summarise(sum(likelihood))} else{likelihood_profile_end <- -Inf}
  
  #### Compare to data 
  compare_dat <- likelihood_lookup_elements + likelihood_lookup_totals + 10 * likelihood_profile_end # add in a 10* weight for profile_end as otherwise only a small contribution relatively
}else{compare_dat <- -Inf}

# return log likelihood
as.numeric(compare_dat) #### - 718 for mock data
likelihood_lookup_elements
likelihood_lookup_totals
10 * likelihood_profile_end

model_outputp <- rename(model_outputp, parent_strain = parent)
model_outputp$name <- recode(model_outputp$variable, V2 = "phi6",V5 = "phi2",V6 = "p1",V7 = "p2",V8 = "p3",V10 ="p4")

g1 <- ggplot(pigg_elements %>% filter(name %in% c("phi6","phi2","p1","p2","p3","p4")),
             aes(x=time, y = sum_prop, group = interaction(name, pig))) +
  geom_line(aes(col = name, linetype = factor(pig)),size = 1.5, alpha = 0.4) +
  geom_line(data = model_outputp, aes(x = time, y = prev, group = interaction(parent_strain,name))) +
  geom_point(aes(col = name),size = 1.5) +
  facet_wrap(name~parent_strain, ncol = 2) + 
  ggtitle(paste0("likelihood = ", round(compare_dat,3)))

g2 <- ggplot(totalsp, aes(x=time, y = value, group = interaction(time,name,parent))) +
  geom_point(aes(col = factor(parent))) +
  geom_line(aes(group = interaction(lim, parent), col = factor(parent)), lty = "dashed") +
  scale_y_log10() +
  geom_line(data = model_outputt, aes(x=time, y = total, group = parent, col = factor(parent)))

g3 <- ggplot(out$all_results %>% filter(variable %in% profiles_needed_end), aes(x=time, y = value)) + geom_line(aes(col = variable)) + geom_vline(xintercept = tsteps) + 
  geom_hline(yintercept = c(0.05) * max(out$all_results$value)) + scale_color_manual(values = c("red","green","blue"), breaks = c(profiles_needed_end))

g1 / (g2 + g3)

ggsave(here::here("fits",paste0(filename,"_","fit",".pdf")))

dput(mcmc.trace[nrow(mcmc.trace),1:(ncol(mcmc.trace)-1)])

