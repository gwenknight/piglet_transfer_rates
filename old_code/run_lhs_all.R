### Run lhs on big range for all lhs

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
library(here)
library(sn) # skewed normal
theme_set(theme_bw(base_size = 11))

setwd(here::here())
source("code/iterative_lhs.R")
source("code/piglet_mrsa_functions.R")

# Scen1
Initial.Values = c(mu = 0.01,
                   gamma = 0.00000001,
                   f = 0.000001, 
                   grow = 0.17, 
                   rel_fit = 0.99)
limit1 = cbind(c(rep(0,2),rep(-0.5,1), rep(0,2)), # lower limits
               c(rep(0.5,2),rep(0.5,1),3,1)) # upper limits 
m1 <- lhs_build_run(Initial.Values, limit = limit1, "lhs_all/sc1", nsamples = 100)
m1$worked <- as.data.frame(m1$worked)
colnames(m1$worked) <- c("ll", names(Initial.Values))
setwd(here::here())
write.csv(m1$worked, "fits/lhs_all/sc1/worked.csv")

# Scen2
Initial.Values = c(mu_phage = 0.01, mu_plasmid = 0.01, 
                   gamma_phage = 0.00000001,gamma_plasmid = 0.00000001,
                   f_phage = 0.000001, f_plasmid = 0.00000001,
                   grow = 0.17, 
                   rel_fit = 0.99)
limit2 = cbind(c(rep(0,4),rep(-1,2), rep(0,2)),c(rep(1,4),rep(1,2),3,1.7))
m2 <- lhs_build_run(Initial.Values, limit = limit2, "lhs_all/sc2", nsamples = 1e4)
m2$worked <- as.data.frame(m2$worked)
colnames(m2$worked) <- c("ll", names(Initial.Values))
setwd(here::here())
write.csv(m2$worked, "fits/lhs_all/sc2/worked.csv")

# Scen3
Initial.Values = c(mu = 0.01,
                   gamma = 0.00000001,
                   f2 = 0.000001, f5 = 0.000001, f6 = 0.000001, 
                   f7 = 0.000001, f8 = 0.000001, f10 = 0.000001, 
                   grow = 0.17, 
                   rel_fit = 0.99)
limit3 = cbind(c(rep(0,2),rep(-0.7,6), rep(0,2)), c(rep(0.7,2),rep(0.7,6),3,1.5))
m3 <- lhs_build_run(Initial.Values, limit = limit3, "lhs_all/sc3", nsamples = 1e4)
m3$worked <- as.data.frame(m3$worked)
colnames(m3$worked) <- c("ll", names(Initial.Values))
setwd(here::here())
write.csv(m3$worked, "fits/lhs_all/sc3/worked.csv")


# Scen4
Initial.Values =  c(mu2 = 0, mu5 = 0, mu6 = 0, 
                    mu7 = 0.65, mu8 = 0, mu10 = 0.1, 
                    gamma2 = 0.01, gamma5 = 0.01, 
                    gamma6 = 0.000000000001, gamma7 = 0.0005, 
                    gamma8 = 0.000000000001, gamma10 = 0.00000001, 
                    f2 = -0.4, f5 = -0.4, f6 = 0.5, 
                    f7 = -0.5, f8 = 0.4, f10 = 0.4, 
                    grow = 0.08, rel_fit = 0.95)
limit4 = cbind(c(rep(0,12),rep(-1,6), rep(0,2)), c(rep(1,12),rep(1,6),3,1.5))
m4 <- lhs_build_run(Initial.Values, limit = limit4, "lhs_all/sc4", nsamples = 100)
m4$worked <- as.data.frame(m4$worked)
colnames(m4$worked) <- c("ll", names(Initial.Values))
setwd(here::here())
write.csv(m4$worked, "fits/lhs_all/sc4/worked.csv")

## No loss 
# Scen4
Initial.Values = c(mu2 = 0, mu5 = 0, mu6 = 0,
                   mu7 = 0, mu8 = 0, mu10 = 0,
                   gamma2 = 0.00005, gamma5 = 0.00005,
                   gamma6 = 0.00005, gamma7 = 0.00005,
                   gamma8 = 0.0000508, gamma10 = 0.00005,
                   f2 = -0.00005, f5 = -0.00005, f6 = 0.00005,
                   f7 = -0.00005, f8 = 0.00005, f10 = 0.00005,
                   grow = 0.12, rel_fit = 1.08)
limit4_nl = cbind(c(rep(0,12),rep(-1,6), rep(0,2)), c(rep(0,6),rep(1,6),rep(1,6),3,1.5))
m4_nl <- lhs_build_run(Initial.Values, limit = limit4_nl, "lhs_all/sc4_nl", nsamples = 1e4)
m4_nl$worked <- as.data.frame(m4_nl$worked)
colnames(m4_nl$worked) <- c("ll", names(Initial.Values))
setwd(here::here())
write.csv(m4_nl$worked, "fits/lhs_all/sc4_nl/worked.csv")


### Look at answers 
### Output
ggplot(m4$worked, aes(x=ll)) + geom_density()
cormat <- round(cor(m4$worked[,1:ncol(m4$worked)]),2)
melted_cormat <- melt(cormat)
ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()
library("PerformanceAnalytics")
my_data <- m4$worked[,2:ncol(m4$worked)]
chart.Correlation(my_data, histogram=TRUE, pch=19)

best <- m4$worked %>% filter(ll == max(m4$worked$ll))
m4$worked %>% filter(ll > round(best$ll,-2))
best <- m4$worked %>% filter(round(ll,2) == round(-595.00,2))
out <- piglet_mrsa_movement(tsteps, best[-1], ini$bacteria, ini$difference_list)


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
  likelihood_profile_end <- sum(log(dnorm(prof_end,mean = c(0.8,0.2,0.2), sd = 0.1)))
  
  #### Add in priors
  # Set up for all - if para not there then 0 
  if(length(best[-1]) < 32 && length(best[-1]) > 10){
    prior.mu = dunif(as.numeric(best[-1]["mu2"]), min = 0, max = 1, log = TRUE) + 
      dunif(as.numeric(best[-1]["mu5"]), min = 0, max = 1, log = TRUE) + 
      dunif(as.numeric(best[-1]["mu6"]), min = 0, max = 1, log = TRUE) + 
      dunif(as.numeric(best[-1]["mu7"]), min = 0, max = 1, log = TRUE) + 
      dunif(as.numeric(best[-1]["mu8"]), min = 0, max = 1, log = TRUE) + 
      dunif(as.numeric(best[-1]["mu10"]), min = 0, max = 1, log = TRUE)
    prior.gamma = dunif(as.numeric(best[-1]["gamma2"]), min = 0, max = 1, log = TRUE) + 
      dunif(as.numeric(best[-1]["gamma5"]), min = 0, max = 1, log = TRUE) + 
      dunif(as.numeric(best[-1]["gamma6"]), min = 0, max = 1, log = TRUE) + 
      dunif(as.numeric(best[-1]["gamma7"]), min = 0, max = 1, log = TRUE) + 
      dunif(as.numeric(best[-1]["gamma8"]), min = 0, max = 1, log = TRUE) + 
      dunif(as.numeric(best[-1]["gamma10"]), min = 0, max = 1, log = TRUE)
    prior.f = dnorm(as.numeric(best[-1]["f2"]), mean = 0, sd = 0.1, log = TRUE) + 
      dnorm(as.numeric(best[-1]["f5"]), mean = 0, sd = 0.1, log = TRUE) + 
      dnorm(as.numeric(best[-1]["f6"]), mean = 0, sd = 0.1, log = TRUE) + 
      dnorm(as.numeric(best[-1]["f7"]), mean = 0, sd = 0.1, log = TRUE) + 
      dnorm(as.numeric(best[-1]["f8"]), mean = 0, sd = 0.1, log = TRUE) + 
      dnorm(as.numeric(best[-1]["f10"]), mean = 0, sd = 0.1, log = TRUE) 
    prior.grow = dunif(as.numeric(best[-1][["grow"]]),0,3,log = TRUE)
    prior.relfit = dnorm(as.numeric(best[-1]["rel_fit"]), mean = 1, sd = 0.1, log = TRUE) 
    log.prior = prior.mu + prior.gamma + prior.f + prior.grow + prior.relfit
  }
  
  # If fixed input - same rates for all elements
  if(length(best[-1]) == 5){
    prior.mu = dunif(as.numeric(best[-1]["mu"]), min = 0, max = 1, log = TRUE) 
    prior.gamma = dunif(as.numeric(best[-1]["gamma"]), min = 0, max = 1, log = TRUE) 
    prior.f = dnorm(as.numeric(best[-1]["f"]), mean = 0, sd = 0.1, log = TRUE) 
    prior.grow = dunif(as.numeric(best[-1][["grow"]]),0,3,log = TRUE)
    prior.relfit = dnorm(as.numeric(best[-1]["rel_fit"]), mean = 1, sd = 0.1, log = TRUE) 
    log.prior = prior.mu + prior.gamma + prior.f + prior.grow + prior.relfit
  }
  
  # If fixed input - same rates for all elements and no fitness cost 
  if(length(best[-1]) == 4){
    prior.mu = dunif(as.numeric(best[-1]["mu"]), min = 0, max = 1, log = TRUE) 
    prior.gamma = dunif(as.numeric(best[-1]["gamma"]), min = 0, max = 1, log = TRUE) 
    prior.grow = dunif(as.numeric(best[-1][["grow"]]),0,3,log = TRUE)
    prior.relfit = dnorm(as.numeric(best[-1]["rel_fit"]), mean = 1, sd = 0.1, log = TRUE) 
    log.prior = prior.mu + prior.gamma + prior.grow + prior.relfit
  }
  
  
  # If fixed input - same rates for phage vs plasmids
  if(length(best[-1]) == 8){
    prior.mu = dunif(as.numeric(best[-1]["mu_phage"]), min = 0, max = 1, log = TRUE)  + dunif(as.numeric(best[-1]["mu_plasmid"]), min = 0, max = 1, log = TRUE) 
    prior.gamma = dunif(as.numeric(best[-1]["gamma_phage"]), min = 0, max = 1, log = TRUE) + dunif(as.numeric(best[-1]["gamma_plasmid"]), min = 0, max = 1, log = TRUE) 
    prior.f = dnorm(as.numeric(best[-1]["f_phage"]), mean = 0, sd = 0.1, log = TRUE) + dnorm(as.numeric(best[-1]["f_plasmid"]), mean = 0, sd = 0.1, log = TRUE)
    prior.grow = dunif(as.numeric(best[-1][["grow"]]),0,3,log = TRUE)
    prior.relfit = dnorm(as.numeric(best[-1]["rel_fit"]), mean = 1, sd = 0.1, log = TRUE) 
    log.prior = prior.mu + prior.gamma + prior.f + prior.grow + prior.relfit
  }
  
  
  # If fixed input - same loss/gain rates different fitness
  if(length(best[-1]) == 10){
    prior.mu = dunif(as.numeric(best[-1]["mu"]), min = 0, max = 1, log = TRUE) 
    prior.gamma = dunif(as.numeric(best[-1]["gamma"]), min = 0, max = 1, log = TRUE) 
    prior.f = dnorm(as.numeric(best[-1]["f2"]), mean = 0, sd = 0.1, log = TRUE) + 
      dnorm(as.numeric(best[-1]["f5"]), mean = 0, sd = 0.1, log = TRUE) + 
      dnorm(as.numeric(best[-1]["f6"]), mean = 0, sd = 0.1, log = TRUE) + 
      dnorm(as.numeric(best[-1]["f7"]), mean = 0, sd = 0.1, log = TRUE) + 
      dnorm(as.numeric(best[-1]["f8"]), mean = 0, sd = 0.1, log = TRUE) + 
      dnorm(as.numeric(best[-1]["f10"]), mean = 0, sd = 0.1, log = TRUE)
    prior.grow = dunif(as.numeric(best[-1][["grow"]]),0,3,log = TRUE)
    prior.relfit = dnorm(as.numeric(best[-1]["rel_fit"]), mean = 1, sd = 0.1, log = TRUE) 
    log.prior = prior.mu + prior.gamma + prior.f + prior.grow + prior.relfit
  }
  
  #### Compare to data 
  compare_dat <- log.prior + likelihood_lookup_elements + likelihood_lookup_totals + 10 * likelihood_profile_end # add in a 10* weight for profile_end as otherwise only a small contribution relatively
}else{compare_dat <- -Inf}

# return log likelihood
compare_dat #### - 718 for mock data
log.prior
likelihood_lookup_elements
likelihood_lookup_totals
10*likelihood_profile_end


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

