###### Find odd behaviour

### Look at output all same

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
theme_set(theme_bw(base_size = 11))

### Data
data <- read.csv("data/data_to_fit.csv")[,-1]
data_6 <- data %>% filter(variable %in% c("V2","V5","V6","V7","V8","V10"))
### Likelihood
dist_like <- read.csv("data/seen_predicted.csv")[,-1]
pigg_elements <- read.csv("data/pigg_elements.csv")[,-1]
totalsp <- read.csv("data/totals_bug.csv")[,-1]

### Total number of bugs  time series
totals <- read.csv("data/totals_bug.csv")[,-1] %>% select(time,value, parent, lim) %>%
  pivot_wider(names_from = lim) %>% mutate(weight = 1/(max - min))

ini <- initial_piglet_setup(384)

### parameters
parameters_every <- c(mu = 10, gamma = 10000,
                      f = 0, grow = 2.93)

parameters_every <- c(mu2 = 0.265289370464156, mu5 = 0.130932829264206, mu6 = 2.58577755844879, 
                      mu7 = 0.481693463586813, mu8 = 3.65375125367999, mu10 = 1.18123827559817, 
                      gamma2 = 0.592694899220884, gamma5 = 1.11892206012185, gamma6 = 3.75164109872491, 
                      gamma7 = 0.71398191791436, gamma8 = 0.240601848692034, gamma10 = 2.51037647414252, 
                      f2 = 0.0164854283546256, f5 = 0.213784972241758, f6 = 0.0382536099278368, 
                      f7 = 0.0400716377754492, f8 = 0.203513266049325, f10 = 0.466694316467684, 
                      grow = 1.97708717931669, rel_fit = 0.3)

parameters_every <- c(mu2 = 0.205586807324394, mu5 = 0.200943634306506, mu6 = 1.43214227025678, 
  mu7 = 0.509097403589266, mu8 = 1.19672115106187, mu10 = 0.3, 
  gamma2 = 0.2193937320621, gamma5 = 0.09100294518291, gamma6 = 0.49881135307676, 
  gamma7 = 0.0924492968656, gamma8 = 0.561307180279803, gamma10 = 0.000001, 
  f2 = 0.0853074236054307, f5 = 0.0626762521572307, f6 = 0.190218563660229, 
  f7 = -0.1, f8 = 0.0240887795738617, f10 = 0.555270549669397, 
  grow = 2.1013935939726, rel_fit = 1)

Initial.Values = c(mu2 = 0.265289370464156, mu5 = 0.130932829264206, mu6 = 2.58577755844879, 
                   mu7 = 0.481693463586813, mu8 = 3.65375125367999, mu10 = 10000, 
                   gamma2 = 2, gamma5 = 1.11892206012185, gamma6 = 3.75164109872491, 
                   gamma7 = 2, gamma8 = 0.240601848692034, gamma10 = 2.51037647414252, 
                   f2 = 0.010, f5 = 0.01, f6 = 0.03, 
                   f7 = 0.0402, f8 = -0.3, f10 = -0.32, 
                   grow = 3, rel_fit = 0.9)

## Run for these parameters
#out <- piglet_mrsa_movement(tsteps, parameters_every, ini$bacteria, ini$difference_list)

out <- piglet_mrsa_movement(tsteps, Initial.Values, ini$bacteria, ini$difference_list)


### fg = 0.0002
### f = 0.1 grow = 5    4.5
### f = 0.0001 grow = 3 NO
### f = 0.1 grow = 3 OK


### Likelihood
if(!is.null(out$prev_predict)){
  
  #### element prevalence from model 
  #c("SCCmec","phi6","SaPI","Tn916","phi2","p1","p2","p3","phi3","p4")
  model_outputp <- out$prev_predict %>% filter(variable %in% c("V2","V5","V6","V7","V8","V10"))
  model_outputp$prev <- round(model_outputp$value,2)
  
  # Check what distribution of n_colonies at this prevalence in the model pig
  distributs <- left_join(model_outputp, dist_like, by = "prev") %>% select(parent, time, variable, prob_all, n_colonies_prev)
  # e.g. to check 
  #distributs %>% filter(name == "p1", parent == 1, n_colonies_prev == 0.9900) %>% summarise(sum(prob_all))
  
  # lookup the probability from this distribution for the data
  likelihood_lookup_elements <- left_join(data_6, distributs, by = c("parent","time","variable","n_colonies_prev")) %>% summarise(sum(log(prob_all)))
  
  #### total bugs output from model 
  model_outputt <- out$totl_predict
  
  likelihood_lookup_totals <- left_join(model_outputt, totals, by = c("time", "parent")) %>% rowwise() %>% mutate(val_in = as.numeric(between(total,min,max))) %>%
    as.data.frame() %>% mutate(likelihood = weight * val_in) %>% summarise(sum(log(likelihood)))
  
  #### Compare to data 
  compare_dat <- likelihood_lookup_elements + likelihood_lookup_totals
}else{compare_dat <- -Inf}

# return log likelihood
compare_dat #### - 718 for mock data
likelihood_lookup_elements
likelihood_lookup_totals


model_outputp <- rename(model_outputp, parent_strain = parent)
model_outputp$name <- recode(model_outputp$variable, V2 = "phi6",V5 = "phi2",V6 = "p1",V7 = "p2",V8 = "p3",V10 ="p4")

g1 <- ggplot(pigg_elements %>% filter(name %in% c("phi6","phi2","p1","p2","p3","p4")),
             aes(x=time, y = sum_prop, group = interaction(name, pig))) +
  geom_line(aes(col = name, linetype = factor(pig)),size = 1.5, alpha = 0.4) +
  geom_line(data = model_outputp, aes(x = time, y = prev, group = interaction(parent_strain,name))) +
  geom_point(aes(col = name),size = 1.5) +
  facet_wrap(name~parent_strain, ncol = 2)

g2 <- ggplot(totalsp, aes(x=time, y = value, group = interaction(time,name,parent))) +
  geom_point(aes(col = factor(parent))) +
  geom_line(aes(group = interaction(lim, parent), col = factor(parent)), lty = "dashed") +
  scale_y_log10() +
  geom_line(data = model_outputt, aes(x=time, y = total, group = parent, col = factor(parent)))

### Everything
ev <- as.data.frame(out$everything)
colnames(ev) <- seq(1:2048)
ev$time <- seq(1:tsteps)
ev_m <- ev %>% pivot_longer(cols = seq(1:2048)) %>% mutate(parent = rep(c(rep(1, ncol(ev)/2),c(rep(2, ncol(ev)/2))),tsteps)) %>% filter(value > 0) %>%
  group_by(time, parent) %>% mutate(total_t = sum(value), prop_t = value/total_t)
g3a <- ggplot(ev_m %>% filter(prop_t > 0.05), aes(x=time, y = value, group = name)) + geom_line(aes(col = name)) + theme(legend.position = "none") + facet_wrap(~parent)

g3 <- ggplot(ev_m %>% filter(prop_t > 0.05), aes(x=time, y = value, group = name))  + geom_bar(stat = "identity", position = "fill", aes(fill = factor(name))) + theme(legend.position = "none") + facet_wrap(~parent)

g1 / (g2 + g3)
