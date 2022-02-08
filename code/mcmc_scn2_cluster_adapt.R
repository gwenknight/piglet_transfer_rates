#### Scenario 2: phage / plasmid split

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

library(fmcmc) #https://uscbiostats.github.io/fmcmc/
library(adaptMCMC)

theme_set(theme_bw(base_size = 11))

#### CODE
source("code/piglet_mrsa_functions.R")

set.seed(42) # to get reproducible results

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
Initial.Values = c(mu_phage = 0.01, mu_plasmid = 0.01, 
                   gamma_phage = 0.00000001,gamma_plasmid = 0.00000001,
                   f_phage = 0.000001, f_plasmid = 0.00000001,
                   grow = 0.17, 
                   rel_fit = 0.99)

# Check non-zero likelihood for start
#out <- piglet_mrsa_movement(tsteps, Initial.Values, ini$bacteria, ini$difference_list)
#run_sim_logPosterior(Initial.Values)

### Try 
out_final <- fmcmc::MCMC(
  initial   = Initial.Values,                      
  fun       = run_sim_logPosterior, 
  nsteps    = 2.2e3,                       # Increasing the sample size
  kernel    = kernel_adapt(freq = 1, warmup = 500, ub = c(rep(0.3,4),rep(0.5,2),3,1.2),
                           lb = c(rep(0,4),rep(-0.5,2), rep(0,2))), 
  thin      = 1
)

# Save output
filename = gsub(c(" "), "_", format(as.POSIXct(Sys.time()), tz = "Europe/London", usetz = TRUE))
filename = gsub(":", "-", filename)

write.csv(out_final, here::here("fits/",paste0("scn2_a_",filename,"_","trace",".csv")))




