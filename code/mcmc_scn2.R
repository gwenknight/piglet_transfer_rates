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
Initial.Values = c(mu_phage = 0.01, mu_plasmid = 0.01, 
                   gamma_phage = 0.00000001,gamma_plasmid = 0.00000001,
                   f_phage = 0.000001, f_plasmid = 0.00000001,
                   grow = 0.17, 
                   rel_fit = 0.99)

# Check non-zero likelihood for start
#out <- piglet_mrsa_movement(tsteps, Initial.Values, ini$bacteria, ini$difference_list)
run_sim_logPosterior(Initial.Values)


out <- fmcmc::MCMC(
  initial = Initial.Values,
  fun     = run_sim_logPosterior,
  nsteps  = 1e3,
  kernel  = kernel_normal(scale = .0000001) 
)
plot(out[,1:3])


out <- fmcmc::MCMC(
  initial = out,
  fun     = run_sim_logPosterior,
  nsteps  = 1e4,
  kernel  = kernel_normal(scale = .0000001),
  conv_checker = convergence_geweke(200)
)
