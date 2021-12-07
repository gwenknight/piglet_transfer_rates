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
source("code/mcmcmh.r") # had to change line172

## All same
init.theta = c(mu = 0.15, gamma = 0.06, f = 0.03, grow = 0.0978)
lower.p <- init.theta
lower.p[] <- 0

mcmc.epi3_79 <- mcmcMH(target = run_sim_logPosterior,
                    limits=list(lower = lower.p),
                    init.theta = c(mu = 0.15, gamma = 0.06, f = 0.03, grow = 0.0978),
                    proposal.sd = c(rep(0.005,3), 0.005),
                    n.iterations = 5000,
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
