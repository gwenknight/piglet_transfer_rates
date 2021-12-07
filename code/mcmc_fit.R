##### Fitting
library(tmvtnorm)
library(tidyverse)

### Data 
data <- read.csv("data_to_fit.csv")[,-1] 
data_6 <- data %>% filter(name %in% c("phi6","phi2","p1","p2","p3","p4"))

### Likelihood 
dist_like <- read.csv("seen_predicted.csv")[,-1]

### Total number of bugs  time series
totals <- read.csv("totals_bug.csv")[,-1] %>% select(time,value, parent, lim) %>% 
  pivot_wider(names_from = lim) %>% mutate(weight = 1/(max - min))

# Just focus on 6 elements that actually move 
#c("phi6","phi2","p1","p2","p3","p4")
# theta = c(# gain
#   mu2 = 0.008,mu5 = 0.008,
#   mu6 = 0.0079,mu7 = 0.043,mu8 = 0.000195,mu10 = 0.0008,
#   # loss
#   gamma2 = 0.004,gamma5 = 0.0009,
#   gamma6 = 0.049,gamma7 = 0.089,gamma8 = 0.65,gamma10 = 0.4,
#   # fitness
#   f2 = 0.2,f5 = 0.06,
#   f6 = 0.48,f7 = 0.06,f8 = 0.03,f10 = 0.38,
#   grow = 0.102)

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
library(fitR)
source("code/mcmcmh.r") # had to change line172
## All same
gain_all <- 0.08
gamma_all <- 0.009
fitness_all <- 0.01



init.theta = c(# gain
  mu2 = gain_all, mu5 = gain_all,
  mu6 = gain_all,mu7 = gain_all,mu8 = gain_all,mu10 = gain_all,
  # loss
  gamma2 = gamma_all,gamma5 = gamma_all,
  gamma6 = gamma_all,gamma7 = gamma_all,gamma8 = gamma_all,gamma10 = gamma_all,
  # fitness
  f2 = fitness_all,f5 = fitness_all,
  f6 = fitness_all,f7 = fitness_all,f8 = fitness_all,f10 = fitness_all,
  grow = 0.1)
lower.p <- init.theta
lower.p[] <- 0

mcmc.epi3_79 <- mcmcMH(target = run_sim_logPosterior,
                    limits=list(lower = lower.p),
                    init.theta = c(# gain
                      mu2 = 0.08, mu5 = 0.08,
                      mu6 = 0.08,mu7 = 0.08,mu8 = 0.08,mu10 = 0.08,
                      # loss
                      gamma2 = 0.009,gamma5 = 0.009,
                      gamma6 = 0.009,gamma7 = 0.009,gamma8 = 0.009,gamma10 = 0.009,
                      # fitness
                      f2 = 0.01,f5 = 0.01,
                      f6 = 0.01,f7 = 0.01,f8 = 0.01,f10 = 0.01,
                      grow = 0.1),
                    proposal.sd = c(rep(0.005,12), rep(0.001,6),0.005),
                    n.iterations = 5000,
                    adapt.size.start = 100,
                    adapt.shape.start = 500,
                    adapt.size.cooling=0.999, 
                    verbose = TRUE)

parameter


my_dLogPosterior <- function(fitmodel, theta, init.state, data) {
  
  # calculate the fitmodel prior for parameter vector theta using
  # fitmodel$dprior, and assign to variable log.prior
  log.prior <- fitmodel$dprior(theta, log = TRUE)
  
  # calculate the log-likelihood of `theta`
  # and `init.state` with respect to the data using `dTrajObs`
  # and assign to a variable `log.likelihood`    
  log.likelihood <- dTrajObs(fitmodel, theta, init.state, data, log = TRUE)
  
  # calulate the log-posterior using the log-prior and log-likelihood
  log.posterior <- log.prior + log.likelihood
  
  return(log.posterior)
  
}

my_dLogPosterior_R0_epi1 <- function(R0) {
  
  return(my_dLogPosterior(fitmodel = SIR,
                          theta = c(R0 = R0, D_inf = 2),
                          init.state = c(S = 999, I = 1, R = 0),
                          data = epi1))
}

mcmc.epi3 <- mcmcMH(target = dLogPosterior_epi3,
                    init.theta = c(R0 = 1),
                    proposal.sd = c(0.01),
                    n.iterations = 1000)
dLogPosterior_epi3 <- function(R0) {
  
  return(my_dLogPosterior(fitmodel = SIR,
                          theta = c(R0 = R0, D_inf = 2),
                          init.state = c(S = 999, I = 1, R = 0),
                          data = epi3))
}  

my_dLogPosterior_R0_epi1 <- function(R0) {
  
  return(my_dLogPosterior(fitmodel = SIR,
                          theta = c(R0 = R0, D_inf = 2),
                          init.state = c(S = 999, I = 1, R = 0),
                          data = epi1))
}
my_dLogPosterior_R0_epi1(R0 = 3)  
my_dLogPosterior_epi3 <- function(theta) {
  
  return(my_dLogPosterior(fitmodel = SIR,
                          theta = theta,
                          init.state = c(S = 999, I = 1, R = 0),
                          data = epi3))
  
}  
mcmc.epi3 <- mcmcMH(target = my_dLogPosterior_epi3,
                    init.theta = c(R0 = 1, D_inf = 2),
                    proposal.sd = c(0.01, 0.1),
                    n.iterations = 1000)
