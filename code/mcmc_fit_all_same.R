##### Fitting
library(tmvtnorm)
library(tidyverse)
library(here)
library(coda)


### Data 
data <- read.csv("data/data_to_fit.csv")[,-1] 
data_6 <- data %>% filter(variable %in% c("V2","V5","V6","V7","V8","V10"))

### Likelihood 
dist_like <- read.csv("data/seen_predicted.csv")[,-1]

### Total number of bugs  time series
totals <- read.csv("data/totals_bug.csv")[,-1] %>% select(time,value, parent, lim) %>% 
  pivot_wider(names_from = lim) %>% mutate(weight = 1/(max - min))

### Functions
source("code/piglet_mrsa_functions.R")
source("code/mcmcmh.r") # had to change line172

## All same
init.theta = c(mu = 0.15, gamma = 0.06, f = 0.03, grow = 0.0978)
lower.p <- init.theta
lower.p[] <- 0

mcmc.epi3_79 <- mcmcMH(target = run_sim_logPosterior,
                    limits=list(lower = lower.p),
                    init.theta = c(mu = 0.15, gamma = 0.06, f = 0.03, grow = 0.0978),
                    proposal.sd = c(rep(0.005,3), 0.005),
                    n.iterations = 2000,
                    adapt.size.start = 100,
                    adapt.shape.start = 500,
                    adapt.size.cooling=0.999, 
                    verbose = TRUE)

# Save output
filename = gsub(c(" "), "_", format(as.POSIXct(Sys.time()), tz = "Europe/London", usetz = TRUE))
filename = gsub(":", "-", filename)

write.csv(mcmc.epi3_79$trace, here::here("fits/all_same/",paste0(filename,"_","trace",".csv")))
write.csv(mcmc.epi3_79$acceptance.rate, here::here("fits/all_same/",paste0(filename,"_","acceptance_rates",".csv")))
write.csv(mcmc.epi3_79$covmat.empirical, here::here("fits/all_same/",paste0(filename,"_","covmat",".csv")))


# # # Look at output
# mcmc.trace <- mcmc(mcmc.epi3_79$trace)
# summary(mcmc.trace)
# acceptanceRate <- 1 - rejectionRate(mcmc.trace)
# acceptanceRate
# plot(mcmc.trace)
# effectiveSize(mcmc.trace)
# 
# 
# mcmc.trace.burned <- burnAndThin(mcmc.trace, burn = 1000)
# plot(mcmc.trace.burned)
# 
# autocorr.plot(mcmc.trace.burned)
# 
# mcmc.trace.burned.thinned <- burnAndThin(mcmc.trace.burned, thin = 5)
# autocorr.plot(mcmc.trace.burned.thinned)
# 
# plotESSBurn(mcmc.trace)


# #### Fit 13th Dec
# gain_all = 0.01894225 
# gamma_all = 0.006881459 
# fitness_all = 0.07604289 
# 
# theta = c(# gain
#   mu2 = gain_all, mu5 = gain_all,
#   mu6 = gain_all,mu7 = gain_all,mu8 = gain_all,mu10 = gain_all,
#   # loss
#   gamma2 = gamma_all,gamma5 = gamma_all,
#   gamma6 = gamma_all,gamma7 = gamma_all,gamma8 = gamma_all,gamma10 = gamma_all,
#   # fitness
#   f2 = fitness_all,f5 = fitness_all,
#   f6 = fitness_all,f7 = fitness_all,f8 = fitness_all,f10 = fitness_all,
#   grow = 0.1492127)
# 
# ## Run for these parameters
# out <- run_sim(tsteps, theta)
# max(out$P_all$time)
# out$error
# 
# ### Data 
# data <- read.csv("data/data_to_fit.csv")[,-1] 
# data_6 <- data %>% filter(name %in% c("phi6","phi2","p1","p2","p3","p4"))
# 
# ### Likelihood 
# dist_like <- read.csv("data/seen_predicted.csv")[,-1]
# 
# ### Total number of bugs  time series
# totals <- read.csv("data/totals_bug.csv")[,-1] %>% select(time,value, parent, lim) %>% 
#   pivot_wider(names_from = lim) %>% mutate(weight = 1/(max - min))
# 
# 
# if(max(out$P_all$time)==tsteps){ # if get to end 
#   #### element prevalence from model 
#   model_outputp <- out$prev_predict %>% filter(name %in% c("phi6","phi2","p1","p2","p3","p4"))
#   model_outputp$prev <- round(model_outputp$prev,2)
#   
#   ## check likelihood works
#   #mock_data <- read_csv("data/mock_perfect_data.csv")[,-1]
#   #model_outputp <- mock_data
#   
#   # Check what distribution of n_colonies at this prevalence in the model pig
#   distributs <- left_join(model_outputp, dist_like, by = "prev") %>% select(parent, time, name, prob_all, n_colonies_prev)
#   # e.g. to check 
#   #distributs %>% filter(name == "p1", parent == 1, n_colonies_prev == 0.9900) %>% summarise(sum(prob_all))
#   
#   # lookup the probability from this distribution for the data
#   likelihood_lookup_elements <- left_join(data_6, distributs, by = c("parent","time","name","n_colonies_prev")) %>% summarise(sum(log(prob_all)))
#   
#   #### total bugs output from model 
#   model_outputt <- out$totl_predict
#   
#   likelihood_lookup_totals <- left_join(model_outputt, totals, by = c("time", "parent")) %>% rowwise() %>% mutate(val_in = as.numeric(between(total,min,max))) %>%
#     as.data.frame() %>% mutate(likelihood = weight * val_in) %>% summarise(sum(log(likelihood)))
#   
#   #### Compare to data 
#   compare_dat <- likelihood_lookup_elements + likelihood_lookup_totals
#   
#   compare_dat
#   # compare_dat <- as.numeric(left_join(model_outputp, data_6,  by = c("time", "name", "prev", "parent")) %>% ungroup() %>% 
#   #                             summarise(sum(log(weighted_prob_all))))  + # log likelihood for elements
#   #   as.numeric(left_join(model_outputt, totals, by = c("time", "parent")) %>% rowwise() %>% mutate(val_in = as.numeric(between(total,min,max))) %>%
#   #                as.data.frame() %>% mutate(likelihood = weight * val_in) %>% summarise(sum(log(likelihood)))) # log likelhood for totals
# }else{compare_dat <- as.numeric(-Inf)}
# # return log likelihood
# compare_dat #### - 718 for mock data
# likelihood_lookup_elements
# likelihood_lookup_totals
# 
# 
# model_outputp <- rename(model_outputp, parent_strain = parent)
# pigg_elements <- read.csv("data/pigg_elements.csv")[,-1]
# 
# g1 <- ggplot(pigg_elements %>% filter(name %in% c("phi6","phi2","p1","p2","p3","p4")), 
#              aes(x=time, y = sum_prop, group = interaction(name, pig))) + 
#   geom_line(aes(col = name, linetype = factor(pig)),size = 1.5, alpha = 0.4) + 
#   geom_line(data = model_outputp, aes(x = time, y = prev, group = interaction(parent_strain,name))) + 
#   geom_point(aes(col = name),size = 1.5) + 
#   facet_wrap(name~parent_strain, ncol = 2)  
# 
# 
# totalsp <- read.csv("data/totals_bug.csv")[,-1]
# g2 <- ggplot(totalsp, aes(x=time, y = value, group = interaction(time,name,parent))) +
#   geom_point(aes(col = factor(parent))) + 
#   geom_line(aes(group = interaction(lim, parent), col = factor(parent)), lty = "dashed") + 
#   scale_y_log10() + 
#   geom_line(data = model_outputt, aes(x=time, y = total, group = parent, col = factor(parent)))
# 
# g1 / g2
# ggsave(here::here("fits/all_same/",paste0(filename,"_","fit",".pdf")))
