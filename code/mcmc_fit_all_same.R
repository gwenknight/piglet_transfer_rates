##### Fitting
library(tmvtnorm)
library(tidyverse)
library(here)
library(coda)


### Data 
data <- read.csv("data/data_to_fit.csv")[,-1] 
data_6 <- data %>% filter(name %in% c("phi6","phi2","p1","p2","p3","p4"))

### Likelihood 
dist_like <- read.csv("data/seen_predicted.csv")[,-1]

### Total number of bugs  time series
totals <- read.csv("data/totals_bug.csv")[,-1] %>% select(time,value, parent, lim) %>% 
  pivot_wider(names_from = lim) %>% mutate(weight = 1/(max - min))

### Functions
source("code/functions.R")

### Calculate log posterior from model at theta parameter values
run_sim_logPosterior <- function(theta){
  tsteps = 384
  
  ## Run for these parameters
  out <- run_sim(tsteps, theta)
  max(out$P_all$time)
  
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
    
  }else{compare_dat <- as.numeric(-Inf)}
  # return log likelihood
  as.numeric(compare_dat)
}

### FROM https://rdrr.io/github/sbfnk/fitR/src/R/mcmc.r
#library(fitR)
source("code/mcmcmh.r") # had to change line172
## All same
init.theta = c(mu = 0.15, gamma = 0.06, f = 0.03, grow = 0.0978)
lower.p <- init.theta
lower.p[] <- 0

mcmc.epi3_79 <- mcmcMH(target = run_sim_logPosterior,
                    limits=list(lower = lower.p),
                    init.theta = c(mu = 0.15, gamma = 0.06, f = 0.03, grow = 0.0978),
                    proposal.sd = c(rep(0.005,3), 0.005),
                    n.iterations = 2,
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


# # Look at output
# mcmc.trace <- mcmc(mcmc.epi3_79$trace)
# summary(mcmc.trace)
# acceptanceRate <- 1 - rejectionRate(mcmc.trace)
# acceptanceRate
# plot(mcmc.trace)
