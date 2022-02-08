#### Scenario 4: all diff

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
Initial.Values = c(mu2 = 2.58492135046923e-05, mu5 = 1.23694298645124e-07, mu6 = 0.0091543430464811, 
                   mu7 = 5.63000577873516e-06, mu8 = 0.0100363784007989, mu10 = 0.00644058325472071, 
                   gamma2 = 1.03260990391265e-05, gamma5 = 1.29510851107398e-06, 
                   gamma6 = 7.43707130531284e-12, gamma7 = 9.90849379073332e-07, 
                   gamma8 = 1.18084441913925e-08, gamma10 = 1.25138645913259e-12, 
                   f2 = -0.222792921010004, f5 = -0.311241161192935, f6 = 0.635171096823821, 
                   f7 = -0.172511628201815, f8 = 0.281506548980894, f10 = 0.781415083339214, 
                   grow = 0.125465546487709, rel_fit = 1.08138114848714)

# Check non-zero likelihood for start
#out <- piglet_mrsa_movement(tsteps, Initial.Values, ini$bacteria, ini$difference_list)
#run_sim_logPosterior(Initial.Values)

### Try 
sd <- read.csv("fits/scn4_sd_lhs3.csv")[,-1]

out_final <- fmcmc::MCMC(
  initial   = Initial.Values,                      
  fun       = run_sim_logPosterior, 
  nsteps    = 2.2e3,                       # Increasing the sample size
  kernel    = kernel_adapt(Sd = sd, freq = 1, warmup = 500, ub = c(rep(0.2,12),rep(0.2,6),3,1.2),
                           lb = c(rep(0,12),rep(-0.2,6), rep(0,2))), 
  thin      = 1
)

# Save output
filename = gsub(c(" "), "_", format(as.POSIXct(Sys.time()), tz = "Europe/London", usetz = TRUE))
filename = gsub(":", "-", filename)

write.csv(out_final, here::here("fits/",paste0("scn4_a_",filename,"_","trace",".csv")))




