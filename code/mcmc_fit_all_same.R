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

## All same
init.theta = c(mu = 0.7953883, gamma = 1.963830, f = 0.003424697, grow = 1.436710)
#init.theta = c(mu = 4.58960662274629, gamma = 0.774851268517529, f = 0.000899333599447348, 
#               grow = 1.39733915676241)
# v different to above! 
#next.theta = c(mu = 0.673698428160987, gamma = 2.92459069176649, f = 0.0444864304993408, 
#               grow = 1.4699165045315)
# home computer 16/12
c(mu = 0.751744625832429, gamma = 2.0345341340617, f = 0.0215217218291819, 
  grow = 1.47016534147811)
# git 15/12 - v different! 
c(mu = 3.83421041070614, gamma = 2.49493900039174, f = 0.00256323172862037, 
  grow = 1.39538895651162)

lower.p <- init.theta
lower.p[] <- 0

mcmc.epi3_79 <- mcmcMH(target = run_sim_logPosterior,
                    limits=list(lower = lower.p),
                    init.theta =  c(mu = 3.83421041070614, gamma = 2.49493900039174, f = 0.00256323172862037, 
                                    grow = 1.39538895651162),
                    proposal.sd = c(rep(0.005,3), 0.005),
                    n.iterations = 10000,
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



