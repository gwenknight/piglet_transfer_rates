##### Fitting
library(tmvtnorm)
library(tidyverse)
library(here)
library(coda)
library(ggplot2)
library(scales)
library(reshape2)
library(dplyr)
library(Rfast)


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

## Phage same, plasmids same 
# init.theta = c(mu_phage = 0.15, gamma_phage = 0.06, f_phage = 0.03, 
#                mu_plasmid = 0.15, gamma_plasmid = 0.06, f_plasmid = 0.03,
#                grow = 0.0978)
init.theta = c(mu_phage = 0.250183096038211, gamma_phage = 1.40119453558188, 
  f_phage = 0.00834685332776346, mu_plasmid = 1.97624042546612, 
  gamma_plasmid = 0.348806989403094, f_plasmid = 0.143797591886829, 
  grow = 1.67736788787884)

lower.p <- init.theta
lower.p[] <- 0

mcmc.epi3_79 <- mcmcMH(target = run_sim_logPosterior,
                       limits=list(lower = lower.p),
                       init.theta = c(mu_phage = 0.250183096038211, gamma_phage = 1.40119453558188, 
                                      f_phage = 0.00834685332776346, mu_plasmid = 1.97624042546612, 
                                      gamma_plasmid = 0.348806989403094, f_plasmid = 0.143797591886829, 
                                      grow = 1.67736788787884),
                       proposal.sd = c(rep(0.005,8)),
                       n.iterations = 2000,
                       adapt.size.start = 100,
                       adapt.shape.start = 500,
                       adapt.size.cooling=0.999, 
                       verbose = TRUE)

# Save output
filename = gsub(c(" "), "_", format(as.POSIXct(Sys.time()), tz = "Europe/London", usetz = TRUE))
filename = gsub(":", "-", filename)

write.csv(mcmc.epi3_79$trace, here::here("fits/phage_plasmid/",paste0(filename,"_","trace",".csv")))
write.csv(mcmc.epi3_79$acceptance.rate, here::here("fits/phage_plasmid/",paste0(filename,"_","acceptance_rates",".csv")))
write.csv(mcmc.epi3_79$covmat.empirical, here::here("fits/phage_plasmid/",paste0(filename,"_","covmat",".csv")))


# # Look at output
# mcmc.trace <- mcmc(mcmc.epi3_79$trace)
# summary(mcmc.trace)
# acceptanceRate <- 1 - rejectionRate(mcmc.trace)
# acceptanceRate
# plot(mcmc.trace)
