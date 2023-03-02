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
library(prodlim)
library(ggforce)
library(ggridges)
library(truncnorm)
library(sn) # skewed normal
library(here)
source("code/piglet_mrsa_functions.R")
theme_set(theme_bw(base_size = 11))

setwd(here::here())

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
#Pinit = t(c(1,1,1,1,1,1,1,1,0,0,1)) # pig adapted
#Qinit = t(c(0,0,0,0,0,0,0,0,1,1,2)) # human 
profile_end1.1 <- c(1,1,1,1,1,1,1,1,0,0,1) # pig - same as at start
profile_end2.1 <- c(0,1,0,0,1,0,0,0,1,0,2) # human - gains phi6 / phi2 / loses p4 (2/5/10)
profile_end2.2 <- c(0,1,0,0,1,0,1,0,1,0,2) # human
profiles_needed_end <- c(row.match(profile_end1.1, as.data.frame(ini$bacteria)%>%select(-freq), nomatch = NA),
  row.match(profile_end2.1, as.data.frame(ini$bacteria)%>%select(-freq), nomatch = NA),
  row.match(profile_end2.2, as.data.frame(ini$bacteria)%>%select(-freq), nomatch = NA))
#c("SCCmec","phi6","SaPI","Tn916","phi2","p1","p2","p3","phi3","p4")
             #2                     5      6    7   8           10

### parameters
parameters_every <- c(mu = 10, gamma = 10000,
                      f = 0, grow = 2.93)

Initial.Values = c(mu2 = 0.0000265289370464156, mu5 = 0.000000130932829264206, mu6 = 0.01, 
                   mu7 = 0.00000481693463586813, mu8 = 0.01, mu10 = 0.005, 
                   gamma2 = 0.00001, gamma5 = 0.000001, gamma6 = 0.00000000001, 
                   gamma7 = 0.000001, gamma8 = 0.00000001, gamma10 = 0.000000000001, 
                   f2 = -0.3, f5 = -0.3, f6 = 0.9, 
                   f7 = -0.2, f8 = 0.3, f10 = 0.8, 
                   grow = 0.15, rel_fit = 1)

Initial.Values = c(mu_phage = 0.01, mu_plasmid = 0.01, 
                   gamma_phage = 0.00000001,gamma_plasmid = 0.00000001,
                   f_phage = 0.000001, f_plasmid = 0.00000001,
                   grow = 0.17, 
                   rel_fit = 0.99)
Initial.Values = c(mu2 = 0.023631663994558, mu5 = 7.96232764582261e-07, mu6 = 0.0362005899163667, 
                   mu7 = 3.84139033003573e-05, mu8 = 0.0516027884447939, mu10 = 0.0202872958267813, 
                   gamma2 = 2.36264639681026e-05, gamma5 = 1.14560454890954e-06, 
                   gamma6 = 3.33363775439206e-10, gamma7 = 6.92065123035721e-07, 
                   gamma8 = 5.05865587238981e-06, gamma10 = 5.34202975235121e-12, 
                   f2 = -2.10459304661928, f5 = -1.86778149978006, f6 = 0.673853238158728, 
                   f7 = -0.334461607733397, f8 = 2.00273114386821, f10 = 0.706988436254512, 
                   grow = 0.0436026095794833, rel_fit = 1.34916656665852)
#c("SCCmec","phi6","SaPI","Tn916","phi2","p1","p2","p3","phi3","p4")
             #2                     5      6    7   8           10
# LIKELIHOOD OF -11723
# Initial.Values = c(mu2 = 0, mu5 = 0, mu6 = 0, 
#                    mu7 = 0.65, mu8 = 0, mu10 = 0.1, 
#                    gamma2 = 0.01, gamma5 = 0.01, 
#                    gamma6 = 0.000000000001, gamma7 = 0.0005, 
#                    gamma8 = 0.000000000001, gamma10 = 0.00000001, 
#                    f2 = -0.4, f5 = -0.4, f6 = 0.5, 
#                    f7 = -0.5, f8 = 0.4, f10 = 0.4, 
#                    grow = 0.08, rel_fit = 0.95)

Initial.Values4 = c(mu2 = 0, mu5 = 0, mu6 = 0, 
                   mu7 = 0.65, mu8 = 0, mu10 = 0.1, 
                   gamma2 = 0.01, gamma5 = 0.01, 
                   gamma6 = 0.000000000001, gamma7 = 0.0005, 
                   gamma8 = 0.000000000001, gamma10 = 0.00000001, 
                   f2 = -0.4, f5 = -0.4, f6 = 0.5, 
                   f7 = -0.5, f8 = 0.4, f10 = 0.4, 
                   grow = 0.08, rel_fit = 0.95)
run_sim_logPosterior(Initial.Values4)


Initial.Values3 = c(mu = 0.01,
                   gamma = 0.00000001,
                   f2 = 0.000001, f5 = 0.000001, f6 = 0.000001, 
                   f7 = 0.000001, f8 = 0.000001, f10 = 0.000001, 
                   grow = 0.17, 
                   rel_fit = 0.99)
run_sim_logPosterior(Initial.Values3)

Initial.Values2 = c(mu_phage = 0.01, mu_plasmid = 0.01, 
                   gamma_phage = 0.00000001,gamma_plasmid = 0.00000001,
                   f_phage = 0.000001, f_plasmid = 0.00000001,
                   grow = 0.17, 
                   rel_fit = 0.99)
run_sim_logPosterior(Initial.Values2)

Initial.Values1 = c(mu = 0.1,
                   gamma = 0.000000000001,
                   f = 0.1, 
                   grow = 0.08, 
                   rel_fit = 0.95)
run_sim_logPosterior(Initial.Values1)

#### All give non infinite outputs


# Initial.Values = c(mu_phage = 0.0000000000001, mu_plasmid = 0.0000000000000001, 
#                    gamma_phage = 0.00001,gamma_plasmid = 0.000000001,
#                    f_phage = -0.01, f_plasmid = -0.00000001,
#                    grow = 0.17, 
#                    rel_fit = 0.99)
#Initial.Values = samp.coda[100,]

## Run for these parameters
#out <- piglet_mrsa_movement(tsteps, parameters_every, ini$bacteria, ini$difference_list)

out <- piglet_mrsa_movement(tsteps, Initial.Values, ini$bacteria, ini$difference_list)

#out <- piglet_mrsa_movement(tsteps, parameters, ini$bacteria, ini$difference_list)

### fg = 0.0002
### f = 0.1 grow = 5    4.5
### f = 0.0001 grow = 3 NO
### f = 0.1 grow = 3 OK


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
  #likelihood_profile_end <- sum(log(dnorm(prof_end,mean = c(0.8,0.2,0.2), sd = 0.5)))
  ## Strain 1 has to be mostly the parent strain! 
  ## Strain 2: two profiles are highly dominant: sum of these has to be > 0.5
  likelihood_profile_end <- log(dnorm(prof_end[1], mean = 0.975, sd = 0.0288)) + log(dtruncnorm(sum(prof_end[2:3]),a = 0, b = 1, mean = 0.75, sd = 0.17))

  #likelihood_profile_end <- sum(log(dsn(prof_end,xi = c(0.735,0.14,0.14), omega = c(0.08,0.2,0.2), alpha = 4))) # Made into a skewed normal - to bias to going over these limits
  # dsn(1,xi = c(0.735,0.14,0.14), omega = c(0.08,0.3,0.3), alpha = 4) # This skewed normal has density < 5% at 1
  # dsn(1,xi = c(0.735,0.14,0.14), omega = c(0.08,0.3,0.3), alpha = 4) # and density < 1% (but > 0 to avoid -Inf) at 0
  # mean(rsn(100000,xi = 0.735, omega = 0.08, alpha = 4)) # and 80%
  # mean(rsn(100000,xi = 0.14, omega = 0.08, alpha = 4)) # or 20% as the mean value.
  # The skew is fixed to be the same for both. 
  
  # Don't make -Inf possible
  # if(nrow(profile_end) == 3 && total_end[2] > 0){
  #   likelihood_profile_end <- profile_end %>% 
  #     mutate(prop = value / c(total_end[1],total_end[2],total_end[2]),
  #            cutoff = pmax(0,prop - 0.05),#c(0.8,0.2,0.2)), # more the better # not using -- too strict
  #            likelihood = log(prop)) %>% summarise(sum(likelihood))} else{likelihood_profile_end <- -Inf}
  
  #### Compare to data 
  compare_dat <- likelihood_lookup_elements + likelihood_lookup_totals + 100 * likelihood_profile_end # add in a 100* weight for profile_end as otherwise only a small contribution relatively
}else{compare_dat <- -Inf}

# return log likelihood: bigger better (-10 better than -100)
compare_dat #### - 718 for mock data
likelihood_lookup_elements
likelihood_lookup_totals
likelihood_profile_end


model_outputp <- rename(model_outputp, parent_strain = parent)
model_outputp$name <- recode(model_outputp$variable, V2 = "phi6",V5 = "phi2",V6 = "p1",V7 = "p2",V8 = "p3",V10 ="p4")

g1 <- ggplot(pigg_elements %>% filter(name %in% c("phi6","phi2","p1","p2","p3","p4")),
             aes(x=time, y = sum_prop, group = interaction(name, pig))) +
  geom_line(aes(col = name, linetype = factor(pig)),size = 1.5, alpha = 0.6) +
  geom_line(data = model_outputp, aes(x = time, y = prev, group = interaction(parent_strain,name)),size = 1) +
  geom_point(aes(col = name),size = 1.5) +
  facet_wrap(name~parent_strain, ncol = 2) +
  scale_color_discrete("MGE") + 
  scale_linetype_discrete("Piglet") + 
  scale_y_continuous("Proportion of population", lim = c(0,1.1)) + 
  ggtitle(paste0("likelihood = ", round(compare_dat,3)))

g2 <- ggplot(totalsp, aes(x=time, y = value, group = interaction(time,name,parent))) +
  geom_point(aes(col = factor(parent))) +
  geom_line(aes(group = interaction(lim, parent), col = factor(parent)), lty = "dashed") +
  scale_y_log10("Total bacterial count") +
  scale_colour_discrete("Parent") + 
  geom_line(data = model_outputt, aes(x=time, y = total, group = parent, col = factor(parent)))

g3 <- ggplot(out$all_results %>% filter(variable %in% profiles_needed_end), aes(x=time, y = value)) + geom_line(aes(col = variable)) + geom_vline(xintercept = tsteps) + 
   geom_hline(yintercept = c(0.8, 0.2) * total_end) + scale_color_manual(values = c("red","green","blue"), breaks = c(profiles_needed_end), "Profile") + 
  scale_y_continuous("Number of bacteria")

#geom_hline(yintercept = c(0.05, 0.2) * max(out$all_results$value))
### Everything
#ev_m <- out$all_results%>% filter(value > 0) %>%
#  group_by(time, parent) %>% mutate(total_t = sum(value), prop_t = value/total_t)
#g3a <- ggplot(ev_m %>% filter(prop_t > 0.05), aes(x=time, y = value, group = variable)) + geom_line(aes(col = variable)) + theme(legend.position = "none") + facet_wrap(~parent)
#
#g3 <- ggplot(ev_m %>% filter(prop_t > 0.05), aes(x=time, y = value, group = variable))  + geom_bar(stat = "identity", position = "fill", aes(fill = interaction(factor(variable), factor(parent)))) + 
#  theme(legend.position = "none") + facet_wrap(~parent)

#c("SCCmec","phi6","SaPI","Tn916","phi2","p1","p2","p3","phi3","p4")
              #2                     5      6    7   8           10

g1 / (g2 + g3) + 
  plot_layout(heights = c(3,1))
#ggsave("plots/example.jpeg")

#out$all_results %>% filter(time == tsteps, value > 0) %>% group_by(parent) %>% arrange("value")
oo <- out$all_results %>% filter(time == tsteps, value > 0) %>% group_by(parent) %>% arrange("value") %>% filter(value == max(value))
oo
ini$bacteria[c(256,1299, 1363, oo$variable),]

out$all_results %>% filter(time == tsteps, variable %in% c(256,1299,1363,768, 2035))
prof_end

#plot_circles(ini$bacteria, out$all_results,"play")
  
# For model fitting 
#left_join(data_6, distributs, by = c("parent","time","variable","n_colonies_prev")) %>% filter(prob_all == 0)
#model_outputp %>% filter(variable == "V10")

## Which profiles at end? 
# ini$bacteria[as.numeric(unlist(out$all_results %>% filter(time == tsteps, value > 0) %>% summarise(unique(variable)))),]
# bugs <- as.data.frame(ini$bacteria)
# bugs$variable = seq(1:nrow(bugs))
# bugs$parent = c(rep(1, nrow(bugs)/2), rep(2, nrow(bugs)/2))
# results <- out$all_results %>% filter(time == tsteps, value > 0) %>% select(variable, value, parent)
# results$variable <- as.numeric(results$variable)
# left_join(results, bugs %>% select(-freq), by = c("variable", "parent")) %>% arrange(value)
# # and which need...
# profiles_needed_end
# profile_end1.1 <- c(1,1,1,1,1,1,1,1,0,0,1) # pig - same as at start
# profile_end2.1 <- c(0,1,0,0,1,0,0,0,1,0,2) # human - gains phi6 / phi2 / loses p4 (2/5/10)
# profile_end2.2 <- c(0,1,0,0,1,0,1,0,1,0,2)
