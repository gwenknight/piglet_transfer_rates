#### Scenario 1: All same loss / gain / fitness costs

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
source("code/piglet_mrsa_functions.R")
theme_set(theme_bw(base_size = 11))
set.seed(42) # to get reproducible results

library(fmcmc) #https://uscbiostats.github.io/fmcmc/
library(coda)
library(adaptMCMC)
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

# parameters
Initial.Values = c(mu = 0.01,
                   gamma = 0.00000001,
                   f = 0.000001, 
                   grow = 0.17, 
                   rel_fit = 0.99)

# Check non-zero likelihood for start
#out <- piglet_mrsa_movement(tsteps, Initial.Values, ini$bacteria, ini$difference_list)
run_sim_logPosterior(Initial.Values)



samp <- MCMC(run_sim_logPosterior, 
             n=1000, 
             init=Initial.Values, 
             scale=c(0.001, 0.000001, 0.000001, 0.0001, 0.000001),
             adapt=TRUE, 
             acc.rate=0.234)

samp.coda <- convert.to.coda(samp)
class(samp.coda)
write.csv(samp.coda, "fits/220125_scn1_second_1000.csv")
## ----------------------
## use functions of package 'coda'
require(coda)
plot(samp.coda)
cumuplot(samp.coda)

### Try 
library(fmcmc)
khaario <- kernel_adapt(freq = 1, warmup = 500)
set.seed(12) 
out_haario_1 <- MCMC(
  init   = Initial.Values,                       
  fun       = run_sim_logPosterior, 
  nsteps    = 1000,    # We will only run the chain for 100 steps                    
  kernel    = khaario, # We passed the predefined kernel
  thin      = 1,       # No thining here
  nchains   = 1L,      # A single chain
  multicore = FALSE    # Running in serial
)

out <- fmcmc::MCMC(
  initial = Initial.Values,
  fun     = run_sim_logPosterior,
  nsteps  = 1e3,
  kernel  = kernel_normal(scale = .00000000000001) 
)
plot(out[,1:3])

out <- fmcmc::MCMC(
  initial = out,
  fun     = run_sim_logPosterior,
  nsteps  = 1e4,
  kernel  = kernel_normal(scale = .000000000000001),
  conv_checker = convergence_geweke(200)
)

set.seed(112) 

out_final <- fmcmc::MCMC(
  initial   = out,                       # Automatically takes the last 2 points
  fun       = run_sim_logPosterior, 
  nsteps    = 5e4,                       # Increasing the sample size
  kernel    = kernel_normal(scale = .000000000000001),
  thin      = 10
)

#### ADAPT
khaario <- kernel_adapt(Sd = 0.000000001,freq = 1, warmup = 500, ub = c(0.2,0.2,0.5,3,1.5),
                        lb = c(0,0,-0.5,rep(0,2)))

set.seed(12) 

out_haario_1 <- fmcmc::MCMC(
  initial   = Initial.Values,                       
  fun       = run_sim_logPosterior, 
  nsteps    = 5000,    
  kernel    = khaario, # We passed the predefined kernel
  thin      = 1,       # No thining here
  nchains   = 1L,      # A single chain
  multicore = FALSE    # Running in serial
)
plot(out_haario_1[,1:5])
